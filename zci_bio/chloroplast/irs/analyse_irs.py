from collections import defaultdict
from datetime import datetime
import statistics
import numpy as np
from itertools import chain
from step_project.common.table.steps import TableStep
from zci_bio.sequences.steps import SequencesStep
from zci_bio.sequences.fetch import do_fetch_sequences
from zci_bio.utils.extract_data import ExtractData
from zci_bio.utils.stat_by_taxonomy import GroupByTaxonomy
from common_utils.value_data_types import sheets_2_excel
from common_utils.properties_db import PropertiesDB


METHOD_NAMES = ('ncbi', 'airpg', 'ge_seq', 'small_d', 'small_d_P', 'small_d_D', 'small_d_all',
                'chloe', 'chloroplot', 'pga', 'pga_sb', 'plann', 'plann_sb', 'org_annotate', 'zci')
METHODS_USE_SEQUENCES = ('ncbi', 'airpg', 'small_d', 'small_d_P', 'small_d_D', 'small_d_all', 'chloroplot',
                         'pga', 'pga_sb', 'plann', 'plann_sb', 'org_annotate', 'zci')
METHODS_SEPARATE_PATH = ('ge_seq', 'chloe')
METHOD_NAMES_RESEARCH = ('chloe', 'chloroplot', 'ge_seq', 'org_annotate', 'pga', 'plann', 'airpg')
_column_types_acc = [
    ('Organism', 'str'), ('Accession', 'seq_ident'),
    ('Created', 'date'), ('Published', 'date'),
    ('Length', 'int'), ('not_dna', 'int')]
_column_types_method = [
    ('Method', 'str'),
    ('IRa_start', 'int'), ('IRa_end', 'int'), ('IRb_start', 'int'), ('IRb_end', 'int'),
    ('IRa_len', 'int'), ('IRb_len', 'int'), ('diff_len', 'int'), ('diff_type', 'str'),
    ('IR_type', 'str'), ('IR_wraps', 'int'),
    ('SSC_len', 'int'), ('LSC_len', 'int'),
    ('replace_num', 'int'), ('replace_sum', 'int'), ('indel_num', 'int'), ('indel_sum', 'int'),
    ('not_dna_ira', 'int'), ('not_dna_irb', 'int'), ('not_dna_irs', 'int')]
_irs_null_row = [None] * (len(_column_types_method) - 1)
_irs_null_row[_column_types_method.index(('IR_type', 'str')) - 1] = 'no'


class _ByYear:
    def __init__(self, methods, data):
        self.num_by_year = defaultdict(int)  # year -> num
        self.by_year = defaultdict(list)     # (year, method) -> list of data
        self.all_years = set()
        self.methods = methods
        for date_str, m_data in data:
            assert date_str
            if isinstance(date_str, str):
                current_year = datetime.fromisoformat(date_str).year
            else:
                current_year = date_str.year
            self.all_years.add(current_year)
            self.num_by_year[current_year] += 1
            for method, m_d in zip(methods, m_data):
                if 'ira' in m_d:
                    self.by_year[(current_year, method)].append(m_d)

    def get_rows(self):
        by = self.by_year.get
        rows = []
        sum_all = 0
        sum_m = dict((m, 0) for m in self.methods)
        min_ir_l_method = dict()
        max_ir_l_method = dict()
        num_not_eq_method = defaultdict(int)
        num_q_method = defaultdict(int)
        for y in sorted(self.all_years):
            n_all = self.num_by_year.get(y, 0)
            sum_all += n_all
            row = [y, n_all]
            d_row = [''] * len(row)
            # Min, max IR length
            # num with diff
            for idx, m in enumerate(self.methods):
                nm = min_ir_l = max_ir_l = num_not_eq = num_q = 0
                if m_ds := by((y, m)):
                    nm = len(m_ds)
                    if not_q := [d for d in m_ds if d['type'][0] != '?']:
                        min_ir_l = min(min(d['ir_lengths']) for d in not_q)
                        max_ir_l = max(max(d['ir_lengths']) for d in not_q)
                        min_ir_l_method[m] = min(min_ir_l_method.get(m, min_ir_l), min_ir_l)
                        max_ir_l_method[m] = max(max_ir_l_method.get(m, max_ir_l), max_ir_l)
                        #
                        num_not_eq = sum(1 for d in not_q if d['type'] != '+')
                        num_not_eq_method[m] += num_not_eq
                    #
                    num_q = sum(1 for d in m_ds if d['type'][0] == '?')
                    num_q_method[m] += num_q
                #
                sum_m[m] += nm
                rows.append((row if idx == 0 else d_row) +
                            [m, nm, round(100 * (nm / n_all), 2) if n_all else 0, min_ir_l, max_ir_l, num_not_eq, num_q])
        #
        row = ['all', sum_all]
        rows.extend((row if idx == 0 else d_row) +
                    [m, sum_m[m], round(100 * (sum_m[m] / sum_all), 2) if sum_all else 0,
                     min_ir_l_method.get(m, 0), max_ir_l_method.get(m, 0), num_not_eq_method[m], num_q_method[m]]
                    for idx, m in enumerate(self.methods))
        return rows

    def get_columns(self):
        return ['year', 'num_sequences', 'method', 'num annotated', '% annotated',
                'Min IR length', 'Max IR length', 'Num not equal', 'Num questionable']


def analyse_irs_collect_needed_data(step_data, table_step, method, seqs_methods, common_db):
    # Finds sequences needed for the analysis. Depends on cached data (PropertiesDB)
    # Creates step with these sequences.
    assert method in ['seqs', 'ge_seq', 'chloe'], method

    if method == 'chloe':
        step = TableStep(table_step.project, step_data, remove_data=True)
    else:
        step = SequencesStep(table_step.project, step_data, remove_data=True)
    table_step.propagate_step_name_prefix(step)

    seq_idents = table_step.get_column_values('ncbi_ident')
    properties_db = PropertiesDB()

    if method == 'seqs':
        to_fetch = properties_db.not_stored_keys1(seq_idents, 'NCBI GenBank data')  # Set
        for m in seqs_methods:
            to_fetch.update(properties_db.not_stored_keys1(seq_idents, f'annotation {m}'))
        do_fetch_sequences(step, to_fetch, common_db)
    elif method == 'ge_seq':
        do_fetch_sequences(step, properties_db.not_stored_keys1(seq_idents, 'annotation ge_seq'), common_db)
    elif method == 'chloe':
        not_in = properties_db.not_stored_keys1(seq_idents, 'annotation chloe')
        step.set_table_data([[s] for s in not_in], [('seq_ident', 'seq_ident')])

    #
    step.save()
    return step


def analyse_irs(step_data, table_step, seqs_step, ge_seq_step, chloe_step, methods, taxa_ranks, taxa_names):
    step = TableStep(table_step.project, step_data, remove_data=True)
    table_step.propagate_step_name_prefix(step)

    seq_ident_2_length = table_step.mapping_between_columns('ncbi_ident', 'length')
    seq_idents = sorted(seq_ident_2_length.keys())
    nc_2_taxid = table_step.mapping_between_columns('ncbi_ident', 'tax_id')
    extract_data = ExtractData(properties_db=PropertiesDB(), sequences_step=seqs_step)
    # seq_ident -> dict with data from NCBI genbak file
    gb_data = extract_data.cache_keys1_genbank_data(seq_idents, seqs_step)

    # Collect sequence data. seq_ident -> dict(seq_values+, annotations+)
    # Dict's store all needed data for latter outputs.
    acc_data = dict(
        (seq_ident, dict(
            seq_ident=seq_ident,
            taxid=nc_2_taxid[seq_ident],
            length=(length := gb['length']),
            organism=(organism := gb.get('organism', '?')),
            first_date=(first_date := gb.get('first_date')),
            update_date=(update_date := gb.get('update_date')),
            not_dna=(not_dna := len(gb.get('not_dna', []))),
            _seq_row=[organism, seq_ident, first_date, update_date, length, len(gb.get('not_dna', []))]))
        for seq_ident, gb in gb_data.items())
    assert all(x in acc_data for x in seq_idents), [x for x in seq_idents if x not in acc_data]

    # Add annotation data. Add values into a dicts. <annotation_method> -> data
    for method in methods:
        # ExtractData method to use.
        m_call = getattr(extract_data, f'cache_keys1_annotation_{method}')
        # find input step
        seq_step = seqs_step
        if method == 'ge_seq':
            seq_step = ge_seq_step
        elif method == 'chloe':
            seq_step = chloe_step
        # Call annotation method and store data
        m_data = m_call(seq_idents, seq_step=seq_step)
        for seq_ident, irs in m_data.items():
            irs['_method_row'] = _irs_2_row(irs)
            acc_data[seq_ident][method] = irs

    # Export data in Excel file, in more sheets
    assert taxa_ranks or taxa_names
    group_bt = GroupByTaxonomy(list(acc_data.values()), ranks=taxa_ranks, names=taxa_names,
                               taxid_attr='taxid')
    g_data = list(group_bt.sorted_nodes_objects(objects_sort=(lambda d: d['organism']),
                                                return_names=True, compact_names=True))

    grouped_columns = group_bt.grouped_columns()
    empty_row_part = [''] * len(grouped_columns)

    # All data
    _excel_columns = grouped_columns + [c for c, _ in (_column_types_acc + _column_types_method)] + ['Same', 'Similar']
    rows = []
    table_rows = []
    not_compact_names = [''] * len(grouped_columns)
    for node_names, objects in g_data:
        nn = node_names
        not_compact_names = [a or b for a, b in zip(nn, not_compact_names)]
        for o in objects:
            seq_row = o['_seq_row']
            d_sr = [''] * len(seq_row)
            same_m = _group_same_irs(o, methods)
            similar_m = _group_similar_irs(o, methods)
            rows.extend(((nn + seq_row) if i == 0 else (empty_row_part + d_sr)) + [m] + o[m]['_method_row'] + [sm, sim]
                        for i, (m, sm, sim) in enumerate(zip(methods, same_m, similar_m)))
            table_rows.extend(not_compact_names + seq_row + [m] + o[m]['_method_row'] + [sm, sim]
                              for i, (m, sm, sim) in enumerate(zip(methods, same_m, similar_m)))
            nn = empty_row_part
    sheets = [('All methods', _excel_columns, rows)]

    # One sheet per method
    _excel_columns = grouped_columns + [c for c, _ in (_column_types_acc + _column_types_method[1:])]
    sheets.extend(
        (m, _excel_columns,
         list(chain(*([(node_names if i == 0 else empty_row_part) + o['_seq_row'] + o[m]['_method_row']
                       for i, o in enumerate(objects)]
                      for node_names, objects in g_data)))) for m in methods)

    # By year
    by_year_fd = _ByYear(methods, ((d['first_date'], (d[m] for m in methods)) for d in acc_data.values()))
    by_year_ud = _ByYear(methods, ((d['update_date'], (d[m] for m in methods)) for d in acc_data.values()))
    columns = by_year_fd.get_columns()
    sheets.append(('Year created', columns, by_year_fd.get_rows()))
    sheets.append(('Year published', columns, by_year_ud.get_rows()))

    # By taxonomy
    for group_idx, rank in [(-1, 'all')] + list(enumerate(grouped_columns)):
        _g_data = g_data if rank == taxa_ranks[-1] else list(
            group_bt.sorted_nodes_objects(objects_sort=(lambda d: d['organism']),
                                          return_names=True, compact_names=True,
                                          lowest_rank=rank if rank in taxa_ranks else group_idx))

        rows = []
        for node_names, objects in _g_data:
            # Taxonomy and lengths
            row_p = node_names + _lengths_row(objects)
            dummy_p = [''] * len(row_p)

            # Methods
            num_seqs = len(objects)
            for idx, m in enumerate(methods):
                if ir_lengths := [max(o[m]['ir_lengths']) for o in objects if 'ira' in o[m]]:
                    num = len(ir_lengths)
                    perc = round(100 * num / num_seqs, 2)
                    _min = min(ir_lengths)
                    _max = max(ir_lengths)
                    avg = round(statistics.mean(ir_lengths), 1)
                    bs = boxplot_stats(ir_lengths)
                else:
                    num = perc = 0
                    _min = _max = avg = None
                    bs = _boxplot_empty
                #
                rows.append((row_p if idx == 0 else dummy_p) + [m, num, perc, _min, _max, avg] + bs)

        columns = grouped_columns[:(group_idx + 1)] + \
            ['Num sequences', 'Min length', 'Max length', 'Length span', 'Avg length', 'Std length',
             'Method', 'Num annotated', '% annotated', 'Min IR length', 'Max IR length', 'Avg IR length'] + \
            _boxplot_columns
        sheets.append((('All summed' if group_idx == -1 else f'By {rank}'), columns, rows))

    # # Comparison
    # sheets.append((
    #   'Comparisons',
    #   [''] + methods[1:],
    #   [[m] + [''] * idx +
    #    [_cm(acc_data, m, x) for x in methods[idx + 1:]] for idx, m in enumerate(methods[:-1])]))

    #
    sheets.append(('Overall stats', ['Stat'] + methods, _stat_rows(methods, acc_data)))

    # Figure data
    rank = 'family'
    if rank in grouped_columns:
        n_methods = [m for m in methods if m != 'airpg']
        group_idx = grouped_columns.index(rank)
        num_species_below = group_bt.ncbi_taxonomy.num_species_below
        fp = extract_data.properties_db.fetch_property
        #
        rows = []
        for nodes, objects in group_bt.sorted_nodes_objects(objects_sort=(lambda d: d['organism']), lowest_rank=rank):
            _min_d = min(o['first_date'] for o in objects)
            not_dna = sum(1 for o in objects if o['not_dna'])
            similar = [0] * len(n_methods)
            for o in objects:
                similar[max(_group_similar_irs(o, n_methods)) - 1] += 1

            rows.append([n.name for n in nodes] +
                        [fp(nodes[-1].name, 'NCBI num species', num_species_below, nodes[-1].taxid)] +
                        _lengths_row(objects) + [_min_d, not_dna] + similar)

        columns = grouped_columns[:(group_idx + 1)] + \
            ['Num species', 'Num sequences', 'Min length', 'Max length', 'Length diff', 'Avg length', 'Std length',
             'Min date', 'Num non DNA'] + \
            [f'With {i}' for i in range(1, len(n_methods) + 1)]
        sheets.append(('Figure data', columns, rows))

    # Sequence data, no clade grouping
    rows = []
    for node_names, objects in group_bt.sorted_nodes_objects(objects_sort=(lambda d: d['organism']),
                                                             return_names=True, compact_names=False):
        rows.extend((node_names + o['_seq_row']) for o in objects)
    sheets.append(('Sequences ', grouped_columns + [c for c, _ in _column_types_acc], rows))

    # Excel: method sheets
    sheets_2_excel('chloroplast_irs_analysis.xlsx', sheets)

    # Store step data. This finishes step
    columns = [(c, 'str') for c in grouped_columns] + _column_types_acc + _column_types_method + \
        [('Same', 'str'), ('Similar', 'int')]
    step.set_table_data(table_rows, columns)
    step.save()

    return step


def _lengths_row(objects):
    lengths = [o['length'] for o in objects]
    _min = min(lengths)
    _max = max(lengths)
    avg = round(statistics.mean(lengths), 1)
    std = round(statistics.stdev(lengths), 1) if len(lengths) > 1 else 0
    return [len(lengths), _min, _max, _max - _min, avg, std]


def _stat_rows(methods, acc_data):
    num_seqs = len(acc_data)
    rows = []

    def _add_stat(title, add_diffs, filter_data):
        per_method = [dict(annot=0, lengths=[], wraps=0, on_start=0,
                           not_dna=0, max_not_dna=0,
                           diff_lens=0, max_diff_len=0, max_indel_len=0) for _ in methods]
        for seq_ident, ir_data in acc_data.items():
            for m, pm in zip(methods, per_method):
                m_data = ir_data[m]
                # assert all(x > 0 for x in m_data.get('ir_lengths', [])), (seq_ident, m, m_data['ir_lengths'])
                if filter_data(m_data):
                    pm['annot'] += 1
                    ir_lengths = m_data['ir_lengths']
                    pm['lengths'].extend(ir_lengths)
                    if any(e < s for s, e in (m_data['ira'], m_data['irb'])):
                        pm['wraps'] += 1
                    if any(s == 0 for s, e in (m_data['ira'], m_data['irb'])):
                        pm['on_start'] += 1
                    if (nd := m_data.get('not_dna')):
                        pm['not_dna'] += 1
                        pm['max_not_dna'] = max(pm['max_not_dna'], *nd)
                    if add_diffs:
                        if dl := abs(ir_lengths[0] - ir_lengths[1]):
                            pm['diff_lens'] += 1
                            pm['max_diff_len'] = max(pm['max_diff_len'], dl)
                        pm['max_indel_len'] = max(pm['max_indel_len'], m_data.get('max_indel_length', 0))

        rows.extend([
            [title] + [''] * len(methods),
            ['Num annotated'] + [d['annot'] for d in per_method],
            ['% annotated'] + [round(100 * d['annot'] / num_seqs, 2) for d in per_method],
            ['Min IR length'] + [min(d['lengths'], default=0)for d in per_method],
            ['Max IR length'] + [max(d['lengths'], default=0) for d in per_method],
            ['Avg IR length'] + [round(statistics.mean(d['lengths']), 1) if d['lengths'] else 0 for d in per_method],
            ['Num wraps'] + [d['wraps'] for d in per_method],
            ['Num on start'] + [d['on_start'] for d in per_method],
            ['Num with not DNA'] + [d['not_dna'] for d in per_method],
            ['Max not DNA length'] + [d['max_not_dna'] for d in per_method],
        ])
        if add_diffs:
            rows.extend([
                ['Num different lengths'] + [d['diff_lens'] for d in per_method],
                ['Max different length'] + [d['max_diff_len'] for d in per_method],
                ['Max indel length'] + [d['max_indel_len'] for d in per_method],
            ])
        #
        box_s = [boxplot_stats(d['lengths']) for d in per_method]
        for idx, bl in enumerate(_boxplot_columns):
            rows.append([bl] + [d[idx] for d in box_s])
    #
    e_row = [''] * (len(methods) + 1)
    _add_stat('Annotated sequences', True, lambda ir_data: 'ira' in ir_data and ir_data['type'][0] != '?')
    rows.append(e_row)
    _add_stat('Exact IRs', False, lambda ir_data: ir_data.get('type') == '+')
    rows.append(e_row)
    _add_stat('IRs differ', True, lambda ir_data: 'ira' in ir_data and ir_data['type'] != '+' and ir_data['type'][0] != '?')
    rows.append(e_row)
    rows.append(['Num sequences', num_seqs] + e_row[2:])
    return rows


def _group_same_irs(obj, methods):
    same = []
    by_method = []
    for m in methods:
        if ira := obj[m].get('ira'):
            d = ira, obj[m].get('irb')
            try:
                ind = same.index(d) + 1
            except ValueError:
                same.append(d)
                ind = len(same)
            by_method.append(ind)
        else:
            by_method.append(None)
    return by_method


def _group_similar_irs(obj, methods, diff_length=20):
    # Not perfect clustering, but works. Note: indexes not located IRs also, as max index
    same = []
    by_method = []
    for m in methods:
        if ira := obj[m].get('ira'):
            s, e = ira
            for ind, (o_s, o_e) in enumerate(same):
                if abs(s - o_s) + abs(e - o_e) <= diff_length:
                    by_method.append(ind + 1)
                    break
            else:
                same.append(ira)
                by_method.append(len(same))  # Indexed form 1!
        else:
            by_method.append(None)
    return [(b or (len(same) + 1)) for b in by_method]  # Set max index to None's


def _cm(acc_data, m1, m2):
    c = _compare_methods(acc_data, m1, m2)
    return ','.join(str(c[k]) for k in ('same', 'in_1', 'in_2', 'longer_1', 'longer_2', 'not_same'))


def _compare_methods(acc_data, m1, m2):
    comp = dict(same=0, in_1=0, in_2=0, longer_1=0, longer_2=0, not_same=0)
    for data in acc_data.values():
        d1 = data[m1]
        d2 = data[m2]
        if d1 == d2:
            comp['same'] += 1
        elif 'ira' in d1:
            if 'ira' in d2:
                # larger_1=0, larger_2=0, not_same=0
                if _longer(d1, d2):
                    comp['longer_1'] += 1
                elif _longer(d2, d1):
                    comp['longer_2'] += 1
                else:
                    comp['not_same'] += 1
            else:
                comp['in_1'] += 1
        elif 'ira' in d2:
            comp['in_2'] += 1
    return comp


def _longer(d1, d2):
    # Returns True if d1 IRs have larger span than d2
    return _longer_ir(d1['ira'], d2['ira']) and _longer_ir(d1['irb'], d2['irb'])


def _longer_ir(ir1, ir2):
    s1, e1 = ir1
    s2, e2 = ir2
    if s1 < e1:
        return s1 <= s2 < e2 <= e1
    # s1 > e1 -> wraps
    if s2 < e2:
        return s1 <= s2
    return s1 <= s2 and e2 <= e1


def _irs_2_row(irs_data):
    if irs_data is None:    # No data, method was not called?
        return ['?'] * (len(_column_types_method) - 1)
    if 'ira' not in irs_data:  # No annotation
        return _irs_null_row
    #
    dt = irs_data.get('type', '')
    num_sum = [0] * 4
    for f in dt.split(';'):
        if f[0] == 'R':
            num_sum[0], num_sum[1] = map(int, f[2:].split(','))
        elif f[0] == 'I':
            num_sum[2], num_sum[3] = map(int, f[2:].split(','))
    ir_type = ('exact' if dt == '+' else 'differs')
    ir_wraps = (irs_data['ira'][1] < irs_data['ira'][0]) or (irs_data['irb'][1] < irs_data['irb'][0])
    not_dna = irs_data.get('not_dna', [0, 0])
    return irs_data['ira'] + irs_data['irb'] + irs_data['ir_lengths'] + \
        [irs_data['diff_len'], dt, ir_type, int(ir_wraps)] + \
        [irs_data['ssc_length'], irs_data['lsc_length']] + \
        num_sum + not_dna + [sum(not_dna)]


#
# check matplotlib method boxplot_stats(), implemented in matplotlib/cbook/__init__.py
_boxplot_columns = ['Median', 'IRQ', 'Outliers', 'Q1', 'Q3', 'LoVal', 'HiVal', 'LoOutliers', 'HiOutliers']
_boxplot_empty = [None] * 9


def boxplot_stats(x, whis=1.5):
    x = np.asarray(x)

    # medians and quartiles
    q1, med, q3 = np.percentile(x, [25, 50, 75])
    iqr = q3 - q1

    # get high extreme
    hival = q3 + whis * iqr
    wiskhi = x[x <= hival]
    whishi = q3 if (len(wiskhi) == 0 or np.max(wiskhi) < q3) else np.max(wiskhi)

    # get low extreme
    loval = q1 - whis * iqr
    wisklo = x[x >= loval]
    whislo = q1 if (len(wisklo) == 0 or np.min(wisklo) > q1) else np.min(wisklo)

    # compute number of outliers
    num_lo = np.count_nonzero(x < whislo)
    num_hi = np.count_nonzero(x > whishi)

    # np.mean(x)
    return [med, iqr, num_lo + num_hi, q1, q3, loval, hival, num_lo, num_hi]
