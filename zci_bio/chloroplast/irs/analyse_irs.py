from collections import defaultdict
from datetime import datetime
import statistics
from step_project.common.table.steps import TableStep
from zci_bio.sequences.steps import SequencesStep
from zci_bio.sequences.fetch import do_fetch_sequences
from zci_bio.utils.extract_data import ExtractData
from zci_bio.utils.stat_by_taxonomy import StatByTaxonomy
from common_utils.value_data_types import sheets_2_excel
from common_utils.properties_db import PropertiesDB
from ..utils import cycle_distance_lt

METHOD_NAMES = ('ncbi', 'ge_seq', 'small_d', 'small_d_P', 'small_d_D', 'small_d_all',
                'chloe', 'chloroplot', 'pga', 'pga_sb', 'plann', 'plann_sb', 'org_annotate')
METHODS_USE_SEQUENCES = ('ncbi', 'small_d', 'small_d_P', 'small_d_D', 'small_d_all', 'chloroplot',
                         'pga', 'pga_sb', 'plann', 'plann_sb', 'org_annotate')
METHODS_SEPARATE_PATH = ('ge_seq', 'chloe')
_column_types_acc = [
    ('accession', 'seq_ident'), ('organism', 'str'), ('first_date', 'date'), ('length', 'int'), ('not_dna', 'int')]
_column_types_method = [
    ('method', 'str'),
    ('IRa_start', 'int'), ('IRa_end', 'int'), ('IRb_start', 'int'), ('IRb_end', 'int'),
    ('IRa_len', 'int'), ('IRb_len', 'int'), ('diff_len', 'int'), ('type', 'str')]


class _ByTaxonomy(StatByTaxonomy):
    def __init__(self, methods, taxids, ranks=None, names=None):
        super().__init__(taxids, ranks=ranks, names=names)
        self.methods = methods

    def stat_attrs(self):
        return self.methods

    def excel_columns(self):
        return ['Min length', 'Max length', 'Avg. length', 'Length stdev'] + self.methods

    def _add(self, node, methods_data, seq_length):
        node['lengths'].append(seq_length)
        for method, data in methods_data.items():
            if 'ira' in data:
                node[method] += 1

    def _init_new_node(self, node):
        super()._init_new_node(node)
        node['lengths'] = []

    def _sum(self, stats):
        ss = super()._sum(stats)
        for s in stats:
            ss['lengths'].extend(s['lengths'])
        return ss

    def _to_row_label(self, s, max_depth, depth, label):
        lengths = s['lengths']
        avg = statistics.mean(lengths)
        std = round(statistics.stdev(lengths), 1) if len(lengths) > 1 else 0
        row = [''] * max_depth + \
            [s['num'], min(lengths), max(lengths), int(round(avg)), std] + \
            [s[a] for a in self.stat_attrs()]
        row[depth] = label
        return row


class _ByYear:
    def __init__(self):
        self.by_year = defaultdict(int)  # (year, method|None) -> num
        self.all_years = set()
        self.current_year = None

    def set_year(self, date_str):
        if not date_str:
            self.current_year = None
        elif isinstance(date_str, str):
            self.current_year = datetime.fromisoformat(date_str).year
        else:
            self.current_year = date_str.year
        self.all_years.add(self.current_year)
        self.by_year[(self.current_year, None)] += 1

    def add_method(self, method):
        self.by_year[(self.current_year, method)] += 1

    def get_rows(self, methods):
        by = self.by_year.get
        rows = []
        sum_all = 0
        sum_m = dict((m, 0) for m in methods)
        for y in sorted(self.all_years):
            n_all = by((y, None), 0)
            sum_all += n_all
            row = [y, n_all]
            for m in methods:
                nm = by((y, m), 0)
                sum_m[m] += nm
                row.extend([nm, round(100 * (nm / n_all), 2) if n_all else 0])
            rows.append(row)
        #
        sum_row = ['all', sum_all]
        for m in methods:
            sum_row.extend([sum_m[m], round(100 * (sum_m[m] / sum_all), 2) if sum_all else 0])
        rows.append(sum_row)
        return rows

    def get_columns(self, methods):
        cs = ['year', 'num_sequences']
        for m in methods:
            cs.extend([m, f'{m} %'])
        return cs


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
        print('chloe', not_in)
        step.set_table_data([[s] for s in not_in], [('seq_ident', 'seq_ident')])

    #
    step.save()
    return step


def analyse_irs(step_data, table_step, seqs_step, ge_seq_step, chloe_step, methods, taxa_ranks, taxa_names):
    step = TableStep(table_step.project, step_data, remove_data=True)
    table_step.propagate_step_name_prefix(step)

    seq_ident_2_length = table_step.mapping_between_columns('ncbi_ident', 'length')
    seq_idents = sorted(seq_ident_2_length.keys())
    extract_data = ExtractData(properties_db=PropertiesDB(), sequences_step=seqs_step)

    # Collect annotation data
    # seq_ident -> dict([annotation_method->data]+)
    acc_data = dict((s, dict()) for s in seq_idents)
    for method in methods:
        m_call = getattr(extract_data, f'cache_keys1_annotation_{method}')
        seq_step = seqs_step
        if method == 'ge_seq':
            seq_step = ge_seq_step
        elif method == 'chloe':
            seq_step = chloe_step
        m_data = m_call(seq_idents, seq_step=seq_step)
        for seq_ident, irs in m_data.items():
            acc_data[seq_ident][method] = irs

    # Format data for storing in table step and Excel file(s)
    # seq_ident -> dict with data from NCBI genbak file
    gb_data = extract_data.cache_keys1_genbank_data(seq_idents, seqs_step)
    table_rows = []
    per_method = [[] for _ in methods]
    by_first_date_year = _ByYear()
    if taxa_ranks or taxa_names:
        nc_2_taxid = table_step.mapping_between_columns('ncbi_ident', 'tax_id')
        by_taxonomy = _ByTaxonomy(methods, nc_2_taxid.values(), ranks=taxa_ranks, names=taxa_names)
    else:
        by_taxonomy = None
    n_acc_c = len(_column_types_acc)
    for seq_ident, data in acc_data.items():
        gb = gb_data[seq_ident]
        first_date = gb.get('first_date')
        by_first_date_year.set_year(first_date)
        row = [seq_ident, gb.get('organism', '?'), first_date, gb['length'], len(gb.get('not_dna', []))]
        for method, pm in zip(methods, per_method):  # Garanties column order!
            m_data = data[method]
            ir_row = _irs_2_row(m_data)
            row.extend([method] + ir_row)
            pm.append(row[:n_acc_c] + ir_row)
            if m_data and 'ira' in m_data:
                by_first_date_year.add_method(method)
        table_rows.append(row)
        #
        if by_taxonomy:
            by_taxonomy.add(nc_2_taxid[seq_ident], dict((m, data[m]) for m in methods), seq_ident_2_length[seq_ident])

    # Table data
    columns = []
    for method in methods:
        columns.extend([(f'{method}_{n}', t) for n, t in _column_types_method])
    step.set_table_data(table_rows, _column_types_acc + columns)
    step.save()

    # Excel: method sheets
    _excel_columns = [c for c, _ in (_column_types_acc + _column_types_method[1:])]
    sheets = [(m, _excel_columns, rows) for m, rows in zip(methods, per_method)]
    sheets.append((
      'Comparisons',
      [''] + methods[1:],
      [[m] + [''] * idx +
       [_cm(acc_data, m, x) for x in methods[idx + 1:]] for idx, m in enumerate(methods[:-1])]))
    sheets.append(('By year', by_first_date_year.get_columns(methods), by_first_date_year.get_rows(methods)))
    if by_taxonomy:
        sheets.append(by_taxonomy.export_sheet('By taxonomy', 0))
    sheets_2_excel('chloroplast_irs_analysis.xls', sheets)

    return step


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
        return ['?'] * 8
    if len(irs_data) == 1:  # Only length, means method didn't find IRs
        return ['-'] * 8
    #
    seq_length = irs_data['length']
    ira = irs_data['ira']
    irb = irs_data['irb']
    ira_l = cycle_distance_lt(*ira, seq_length)
    irb_l = cycle_distance_lt(*irb, seq_length)
    diff_len = abs(ira_l - irb_l)
    return [ira[0], ira[1], irb[0], irb[1], ira_l, irb_l, diff_len, irs_data['type']]
