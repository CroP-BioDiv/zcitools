import datetime
from collections import defaultdict
from step_project.base_workflow import BaseWorkflow
from common_utils.exceptions import ZCItoolsValueError


class IRsStatistics(BaseWorkflow):
    _WORKFLOW = 'irs_statistics'

    @staticmethod
    def required_parameters():
        # , 'plastids'
        return ('taxons', 'methods')

    @staticmethod
    def format_parameters(params):
        # Methods
        from ..chloroplast.irs.analyse_irs import METHOD_NAMES
        methods = set()
        not_known_methods = set()
        for m in params['methods'].lower().split(','):
            if m == 'all':
                methods.update(METHOD_NAMES)
            elif m in METHOD_NAMES:
                methods.add(m)
            else:
                not_known_methods.add(m)
        if not_known_methods:
            raise ZCItoolsValueError(f'Not know method(s): {", ".join(not_known_methods)}!')
        params['methods'] = methods
        #
        params['plastids'] = int(params.get('plastids', 0))
        params['remove_irl'] = int(params.get('remove_irl', 0))
        if taxa_ranks := params.get('taxa_ranks'):
            params['taxa_ranks'] = taxa_ranks.split(',')
        if 'max_update_date' in params:
            params['max_update_date'] = datetime.date.fromisoformat(params['max_update_date'])
        return params

    def _actions(self):
        from ..chloroplast.irs.analyse_irs import METHODS_USE_SEQUENCES, METHODS_SEPARATE_PATH

        taxons = ' '.join(f'-t {t}' for t in self.parameters['taxons'].split(','))
        plastids = '-P' if self.parameters['plastids'] else ''
        if max_update_date := self.parameters.get('max_update_date', ''):
            max_update_date = f'--max-update-date {max_update_date}'
        remove_irl = '--remove-irl' if int(self.parameters.get('remove_irl', '0')) else ''
        methods = self.parameters['methods']

        actions = [('01_chloroplast_list', f"ncbi_chloroplast_list {taxons} {plastids} {max_update_date} {remove_irl}")]

        # Collect data
        # Methods that use NCBI sequences from common step
        stats = [f'-m {m}' for m in methods]
        seqs_methods = [m for m in methods if m in METHODS_USE_SEQUENCES]
        seqs_methods = ' '.join(f'-s {m}' for m in seqs_methods)
        actions.append(('02_seqs', f"analyse_irs_collect_needed_data 01_chloroplast_list seqs {seqs_methods}"))

        # Methods that use separate path to collect data
        for m in methods:
            if m in METHODS_SEPARATE_PATH:
                stats.append(f'-{m[0]} 03_{m}')
                actions.append((f'02_{m}', f"analyse_irs_collect_needed_data 01_chloroplast_list {m}"))
                actions.append((f'03_{m}', f"{m} 02_{m}"))

        # Analysis
        cmd = f"analyse_irs 01_chloroplast_list 02_seqs {' '.join(stats)}"
        if ranks := self.parameters.get('taxa_ranks'):
            cmd += ' ' + ' '.join(f'-r {r}' for r in ranks)
        if names := self.parameters.get('taxons'):
            cmd += ' ' + ' '.join(f'-n {n}' for n in names.split(','))
        actions.append(('04_stats', cmd))

        # Summary, result, ...
        return actions

    def get_summary(self):
        # ---------------------------------------------------------------------
        # Collect sequence data
        # ---------------------------------------------------------------------
        if not (step_01 := self.project.read_step_if_in('01_chloroplast_list')):
            return dict(text='Project not started!')

        if not (results := self.project.read_step_if_in('04_stats')):
            return dict(text='Analysis not done!')

        #
        methods = self.parameters['methods']
        max_l = max(len(m) for m in methods)
        it_2_idx = dict(exact=0, differs=1, no=2)
        m_2_it = dict((m, [0, 0, 0]) for m in methods)
        m_2_wraps = dict((m, [0, 0]) for m in methods)
        #
        m_2_dl, dl_splits = _idx_data(methods, (0, 10, 100))
        m_2_blocks, block_splits = _idx_data(methods, (1, 10, 100))
        m_2_ddd, ddd_splits = _idx_data(methods, (1, 10, 100))
        for method, ir_type, ir_wraps, diff_len, replace_num, replace_sum, indel_num, indel_sum in results.select(
                ('Method', 'IR_type', 'IR_wraps', 'diff_len', 'replace_num', 'replace_sum', 'indel_num', 'indel_sum')):
            m_2_it[method][it_2_idx[ir_type]] += 1
            if ir_type != 'no':       # With IRs
                m_2_wraps[method][int(ir_wraps)] += 1
                m_2_dl[method][_idx(diff_len, dl_splits)] += 1
            if ir_type == 'differs':  # Not exact IRs
                m_2_blocks[method][_idx(replace_num + indel_num, block_splits)] += 1
                m_2_ddd[method][_idx(replace_sum + indel_sum, block_splits)] += 1

        #
        text = "Summary\n"
        text += self._data_table(m_2_it, 'Number of Sequences with Annotated IRs', ('Exact IRs', 'IRs differ', 'No IRs'))
        text += self._data_table(m_2_wraps, 'IR wraps', ('No', 'Yes'))
        text += self._data_table(m_2_dl, 'IR difference in length', _idx_labels(dl_splits, measure=' bp'))
        text += self._data_table(m_2_blocks, 'IR indel/replace number of blocks', _idx_labels(block_splits))
        text += self._data_table(m_2_ddd, 'IR indel/replace number of bps', _idx_labels(ddd_splits))
        return dict(text=text)

    def _data_table(self, data, title, labels):
        methods = self.parameters['methods']
        max_l = max(len(m) for m in methods)
        text = f"""\n{title}:
             {' '.join(m.rjust(max_l) for m in methods)}
"""
        for idx, lab in enumerate(labels):
            text += f"  {lab:<10} {' '.join(str(data[m][idx]).rjust(max_l) for m in methods)}\n"
        return text


def _idx_data(methods, lens):
    n = len(lens) + 1
    return dict((m, [0] * n) for m in methods), lens


def _idx(x, lens):
    for idx, l in enumerate(lens):
        if x <= l:
            return idx
    return len(lens)


def _idx_labels(lens, measure=''):
    lens = [str(l) for l in lens]
    max_l = len(lens[-1]) + 1
    labels = []
    for l in lens:
        if l == '0':
            labels.append(f'   {l.rjust(max_l)}')
        else:
            labels.append(f'<= {l.rjust(max_l)}')
    labels.append(f'>   {lens[-1]}')
    return [l + measure for l in labels] if measure else labels
