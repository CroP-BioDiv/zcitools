from collections import namedtuple
from step_project.common.table.steps import TableStep
from common_utils.exceptions import ZCItoolsValueError
from zci_bio.chloroplast.utils import find_chloroplast_irs
from zci_bio.utils.features import Feature


class _Characterization(namedtuple('_Characterization', 'len_ira, len_irb, len_diff')):
    def __new__(cls, len_ira, len_irb):
        return super(_Characterization, cls).__new__(cls, len_ira, len_irb, abs(len_ira - len_irb))

    len_irb_x = property(lambda self: self.len_irb if self.len_ira != self.len_irb else None)

    def is_problematic(self):
        return self.len_diff > 0


def get_methods():
    return AnalyseIRs.get_methods()


def analyse_irs(step_data, annotations_step, methods, only_problematic):
    step = TableStep(annotations_step.project, step_data, remove_data=True)
    annotations_step.propagate_step_name_prefix(step)
    sheets = AnalyseIRs(step, annotations_step).run(methods, only_problematic)
    step.save()
    print(sheets)
    # step.to_excel('chloroplast_irs_analysis.xls')  # Store data for a check
    return step


class AnalyseIRs:
    def __init__(self, step, annotations_step):
        project = step.project
        self.step = step
        self.annotations_step = annotations_step  # GeSeq annotations
        self.table_step = project.find_previous_step_of_type(self.annotations_step, 'table')
        self.sequences_step = project.find_previous_step_of_type(self.annotations_step, 'sequences')  # NCBI annotations

    @classmethod
    def get_method_names(cls):
        return sorted(m[8:] for m in dir(cls) if m.startswith('_method_'))

    def get_method_callables(self, methods):
        methods_fs = []
        for m in methods:
            method_name = f'_method_{m.lower()}'
            if not hasattr(self, method_name):
                m_names = ", ".join(self.get_method_names())
                raise ZCItoolsValueError(f'No function for method {m}! Available methods: {m_names}')
            methods_fs.append(getattr(self, method_name))
        return methods_fs

    def run(self, methods, only_problematic):
        seq_idents = sorted(self.annotations_step.all_sequences())
        # ToDo: Caching?
        rows = []
        sheets = []
        for name, m_callable in zip(methods, self.get_method_callables(methods)):
            sheets.append((name, []))
            sheet_rows = sheets[-1][1]
            for seq_ident in seq_idents:
                if irs := m_callable(seq_ident):
                    ira, irb = irs
                    # ToDo: moooooore
                    c = _Characterization(len(ira), len(irb))
                    rows.append([name, seq_ident] + [c.len_ira, c.len_irb_x, c.len_diff])
                    if not only_problematic or c.is_problematic():
                        sheet_rows.append(rows[-1])
                else:
                    rows.append([name, seq_ident] + [None] * 3)
                    if not only_problematic:
                        sheet_rows.append(rows[-1])

        self.step.set_table_data(
            rows,
            [('Method', 'str'), ('seq_ident', 'seq_ident'),
             ('IRa_len', 'int'), ('IRb_len', 'int'), ('diff_len', 'int')])

        return sheets

    # Methods return IRs (tuple (ira, irb)) or None
    def _method_ge_seq(self, seq_ident):
        return self._find_in_seq_record(self.annotations_step.get_sequence_record(seq_ident))

    def _method_ncbi(self, seq_ident):
        return self._find_in_seq_record(self.sequences_step.get_sequence_record(seq_ident))

    def _find_in_seq_record(self, seq_rec):
        if irs := find_chloroplast_irs(seq_rec, check_length=False):
            ira, irb = irs
            s_l = len(seq_rec)
            return Feature(s_l, name='ira', feature=ira), Feature(s_l, name='irb', feature=irb)
