import datetime
from step_project.base_workflow import BaseWorkflow


class IRsStatistics(BaseWorkflow):
    _WORKFLOW = 'irs_statistics'

    @staticmethod
    def required_parameters():
        # , 'plastids'
        return ('taxons', 'methods')

    @staticmethod
    def format_parameters(params):
        from ..chloroplast.irs.analyse_irs import METHOD_NAMES
        params['methods'] = params['methods'].lower().split(',')
        assert all(m in METHOD_NAMES for m in params['methods']), (params['methods'], METHOD_NAMES)
        params['plastids'] = int(params.get('plastids', 0))
        if 'max_update_date' in params:
            params['max_update_date'] = datetime.date.fromisoformat(params['max_update_date'])
        return params

    def _actions(self):
        taxons = ' '.join(f'-t {t}' for t in self.parameters['taxons'].split(','))
        plastids = '-P' if self.parameters['plastids'] else ''
        if max_update_date := self.parameters.get('max_update_date', ''):
            max_update_date = f'--max-update-date {max_update_date}'
        methods = self.parameters['methods']

        actions = [('01_chloroplast_list', f"ncbi_chloroplast_list {taxons} {plastids} {max_update_date}")]

        # Collect data
        seqs_methods = [m for m in methods if m == 'ncbi' or m.startswith('small_d')]
        seqs_methods = ' '.join(f'-s {m}' for m in seqs_methods)
        actions.append((f'02_seqs', f"analyse_irs_collect_needed_data 01_chloroplast_list seqs {seqs_methods}"))

        if 'ge_seq' in methods:
            actions.append((f'02_ge_seq', f"analyse_irs_collect_needed_data 01_chloroplast_list ge_seq"))
            actions.append((f'03_ge_seq', f"ge_seq 02_ge_seq"))

        # Find IRS
        stats = ['-m ge_seq -g 03_ge_seq' if m == 'ge_seq' else f'-m {m}' for m in methods]
        actions.append((f'04_stats', f"analyse_irs 01_chloroplast_list 02_seqs {' '.join(stats)}"))
        # Summary, result, ...
        return actions

    def get_summary(self):
        # ---------------------------------------------------------------------
        # Collect sequence data
        # ---------------------------------------------------------------------
        if not (step := self.project.read_step_if_in('01_chloroplast_list')):
            return dict(text='Project not started!')

        text = ''
        return dict(text=text)
