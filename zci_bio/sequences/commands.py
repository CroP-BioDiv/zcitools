from step_project.base_commands import CreateStepCommand


class FetchSequencesStep(CreateStepCommand):
    _COMMAND = 'fetch_seqs'
    _HELP = "Creates sequences step. Mandatory argument is a table step."
    _STEP_BASE_NAME = 'seqs'

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Input table step')

    def _prev_steps(self):
        return [self.args.step]

    def cache_identifier(self):
        return dict(static=True, data_identifier=['sequences'])

    def run(self, step_data):
        from .fetch import fetch_sequences
        step = self.project.read_step(self.args.step, check_data_type='table')
        return fetch_sequences(step_data, step, self.get_cache_object())
