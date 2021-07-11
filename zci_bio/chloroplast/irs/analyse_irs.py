from collections import namedtuple, Counter
from step_project.common.table.steps import TableStep
from zci_bio.sequences.steps import SequencesStep
from common_utils.exceptions import ZCItoolsValueError
from zci_bio.chloroplast.utils import find_chloroplast_irs
from zci_bio.sequences.fetch import do_fetch_sequences
# from zci_bio.utils.features import Feature
from zci_bio.utils.extract_data import ExtractData
from common_utils.value_data_types import sheets_2_excel
from common_utils.properties_db import PropertiesDB
from ..utils import cycle_distance_lt

METHOD_NAMES = ['ncbi', 'ge_seq', 'small_d']
_column_types_acc = [('accession', 'seq_ident'), ('organism', 'str'), ('first_date', 'date'), ('length', 'int')]
_column_types_method = [
    ('method', 'str'),
    ('IRa_start', 'int'), ('IRa_end', 'int'), ('IRb_start', 'int'), ('IRb_end', 'int'),
    ('IRa_len', 'int'), ('IRb_len', 'int'), ('diff_len', 'int'), ('type', 'str')]


def analyse_irs_collect_needed_data(step_data, table_step, method, seqs_methods, common_db):
    # Finds sequences needed for the analysis. Depends on cached data (PropertiesDB)
    # Creates step with these sequences.
    assert method in ['seqs', 'ge_seq'], method

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

    #
    step.save()
    return step


def analyse_irs(step_data, table_step, seqs_step, ge_seq_step, methods):
    step = TableStep(table_step.project, step_data, remove_data=True)
    table_step.propagate_step_name_prefix(step)

    table_step.mapping_column_2_columns('ncbi_ident', 'length')  # Mapiranje seq_ident -> sto?

    seq_idents = sorted(table_step.get_column_values('ncbi_ident'))
    extract_data = ExtractData(properties_db=PropertiesDB(), sequences_step=seqs_step)

    # seq_ident -> dict([annotation_method->data]+)
    acc_data = dict((s, dict()) for s in seq_idents)
    gb_data = extract_data.cache_keys1_genbank_data(seq_idents, seqs_step)

    for method in methods:
        m_call = getattr(extract_data, f'cache_keys1_annotation_{method}')
        m_data = m_call(seq_idents, seq_step=(ge_seq_step if method == 'ge_seq' else seqs_step))
        for seq_ident, irs in m_data.items():
            acc_data[seq_ident][method] = irs

    # Store step data
    table_rows = []
    per_method = [[] for _ in methods]
    for seq_ident, data in acc_data.items():
        gb = gb_data[seq_ident]
        row = [seq_ident, gb.get('organism', '?'), gb.get('first_date'), gb['length']]
        for method, pm in zip(methods, per_method):  # Garanties column order!
            ir_row = _irs_2_row(data[method])
            row.extend([method] + ir_row)
            pm.append(row[:4] + ir_row)
        table_rows.append(row)

    columns = []
    for method in methods:
        columns.extend([(f'{method}_{n}', t) for n, t in _column_types_method])
    step.set_table_data(table_rows, _column_types_acc + columns)
    step.save()

    # Excel: sheet-ovi metode
    _excel_columns = [c for c, _ in (_column_types_acc + _column_types_method[1:])]
    sheets = [(m, _excel_columns, rows) for m, rows in zip(methods, per_method)]
    sheets_2_excel('chloroplast_irs_analysis.xls', sheets)

    # Excel: Accession -> koje metode su uspjele

    return step


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


# # _calc methods return dict with attrs: length, ira, irb, desc.
# # ira, irb and are pairs (start, end) of integers.
# # Indices are indexed from 0, and end index means that range is < end.
# # It is possible to have start > end.
# def _calc_ncbi(step, seq_ident):
#     return _find_in_seq_record(seq_ident, step.get_sequence_record(seq_ident))


# def _calc_ge_seq(step, seq_ident):
#     return _find_in_seq_record(seq_ident, step.get_sequence_record(seq_ident))


# def _calc_small_d(step, seq_ident):
#     seq_rec = step.get_sequence_record(seq_ident)
#     return dict(length=len(seq_rec.seq))  # ToDo:


# def _find_in_seq_record(seq_ident, seq_rec):
#     if not seq_rec:
#         raise ZCItoolsValueError(f'No sequence record for {seq_ident}!')

#     if irs := find_chloroplast_irs(seq_rec, check_length=False):
#         ira, irb = irs
#         # To be sure!
#         ira_p = ira.location.parts
#         irb_p = irb.location.parts
#         d = dict(length=len(seq_rec.seq),
#                  ira=[int(ira_p[0].start) - 1, int(ira_p[-1].end)],
#                  irb=[int(irb_p[0].start) - 1, int(irb_p[-1].end)])
#         #
#         ira_s = ira.extract(seq_rec)
#         irb_s = irb.extract(seq_rec)
#         if ira.strand == irb.strand:
#             irb_s = irb_s.reverse_complement()
#         d.update(_irs_desc(ira_s, irb_s))
#         return d
#     return dict(length=len(seq_rec.seq))


# def _irs_desc(ira, irb):
#     ira = str(ira.seq)
#     irb = str(irb.seq)
#     # print(len(ira), len(irb))
#     # print('  ira', ira[:20], ira[-20:])
#     # print('  irb', irb[:20], irb[-20:])
#     if ira == irb:
#         return dict(type='+')

#     # Quite slow method. Not unbearable, but slow.
#     # ToDo: possible speedup is to remove equal ends of sequences.
#     print('eva', len(ira), len(irb))
#     diff = difflib.SequenceMatcher(a=ira, b=irb, autojunk=False)
#     opcodes = diff.get_opcodes()
#     # R, I -> [num blocks, num bps]
#     RM = [0, 0]
#     IN = [0, 0]
#     for x in opcodes:
#         print(x)
#         if x[0] == 'equal':
#             continue
#         if x[0] == 'replace':
#             c, bp = RM, (x[2] - x[1])
#         elif x[0] == 'delete':
#             c, bp = IN, (x[2] - x[1])
#         elif x[0] == 'insert':
#             c, bp = IN, (x[4] - x[3])
#         c[0] += 1
#         c[1] += bp

#     return dict(type=';'.join(f'{label}:{d[0]},{d[1]}' for d, label in ((RM, 'R'), (IN, 'I')) if d[0]),
#                 diff=opcodes)


# ############################################################################################
# class _Characterization(namedtuple('_Characterization', 'len_ira, len_irb, len_diff')):
#     def __new__(cls, irs):
#         ira, irb = irs
#         len_ira, len_irb = len(ira), len(irb)
#         return super(_Characterization, cls).__new__(cls, len_ira, len_irb, abs(len_ira - len_irb))

#     len_irb_x = property(lambda self: self.len_irb if self.len_ira != self.len_irb else None)

#     def is_problematic(self):
#         return self.len_diff > 0


# def get_methods():
#     return AnalyseIRs.get_methods()


# class AnalyseIRs:
#     def __init__(self, step, annotations_step):
#         project = step.project
#         self.step = step
#         self.annotations_step = annotations_step  # GeSeq annotations
#         self.table_step = project.find_previous_step_of_type(self.annotations_step, 'table')
#         self.sequences_step = project.find_previous_step_of_type(self.annotations_step, 'sequences')  # NCBI annotations

#     def get_method_callables(self, methods):
#         methods_fs = []
#         for m in methods:
#             method_name = f'_method_{m.lower()}'
#             if not hasattr(self, method_name):
#                 m_names = ", ".join(sorted(m[8:] for m in dir(cls) if m.startswith('_method_')))
#                 raise ZCItoolsValueError(f'No function for method {m}! Available methods: {m_names}')
#             methods_fs.append(getattr(self, method_name))
#         return methods_fs

#     def run(self, methods, export_all):
#         seq_idents = sorted(self.annotations_step.all_sequences())
#         # ToDo: Caching?
#         rows = []
#         sheets = []
#         for name, m_callable in zip(methods, self.get_method_callables(methods)):
#             sheet_rows = []
#             sheets.append((name, sheet_rows))
#             found_irs = 0
#             problematic = 0
#             for seq_ident in seq_idents:
#                 if irs := m_callable(seq_ident):
#                     found_irs += 1
#                     # ToDo: moooooore
#                     c = _Characterization(irs)
#                     rows.append([name, seq_ident] + [c.len_ira, c.len_irb_x, c.len_diff])
#                     if c.is_problematic():
#                         problematic += 1
#                         sheet_rows.append(rows[-1])
#                     elif export_all:
#                         sheet_rows.append(rows[-1])
#                 else:
#                     rows.append([name, seq_ident] + [None] * 3)
#                     if export_all:
#                         sheet_rows.append(rows[-1])

#             # Add summary
#             sheet_rows.append([])
#             sheet_rows.append([None, 'Num sequences', len(seq_idents)])
#             sheet_rows.append([None, 'Num IRs annotated', found_irs])
#             sheet_rows.append([None, 'Num problematic', problematic])

#         self.step.set_table_data(rows, _column_types)
#         return sheets

#     # Methods return IRs (tuple (ira, irb)) or None
#     def _method_ge_seq(self, seq_ident):
#         return self._find_in_seq_record(self.annotations_step.get_sequence_record(seq_ident))

#     def _method_ncbi(self, seq_ident):
#         return self._find_in_seq_record(self.sequences_step.get_sequence_record(seq_ident))

#     def _find_in_seq_record(self, seq_rec):
#         if irs := find_chloroplast_irs(seq_rec, check_length=False):
#             ira, irb = irs
#             s_l = len(seq_rec)
#             return Feature(s_l, name='ira', feature=ira), Feature(s_l, name='irb', feature=irb)
