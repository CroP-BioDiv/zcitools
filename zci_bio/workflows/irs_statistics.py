from step_project.base_workflow import BaseWorkflow


class IRsStatistics(BaseWorkflow):
    _WORKFLOW = 'irs_statistics'

    @staticmethod
    def required_parameters():
        return ('taxons', 'methods')

    def _actions(self):
        taxons = ' '.join(f'-t {t}' for t in self.parameters['taxons'].split(','))
        methods = self.parameters['methods']
        return [
            ('01_chloroplast_list', f"ncbi_chloroplast_list {taxons}"),
            ('02_seqs', 'fetch_seqs 01_chloroplast_list'),
            ('03_GeSeq', 'ge_seq 02_seqs'),  # ?
            ('04_AnalyseIRs', ('analyse_irs', '03_GeSeq', '-m', methods)),
            # ('04_AnalyseIRs', f'analyse_irs 03_GeSeq -m "{methods}"')]
        ]

    def get_summary(self):
        # ---------------------------------------------------------------------
        # Collect sequence data
        # ---------------------------------------------------------------------
        if not (step := self.project.read_step_if_in('01_chloroplast_list')):
            return dict(text='Project not started!')

        text = ''
        return dict(text=text)
