from .base import _CreateStepCommand

# Check base.py for description


class FetchSequencesStep(_CreateStepCommand):
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
        from ..processing.sequence.fetch import fetch_sequences
        from ..steps import read_step
        step = read_step(self.args.step, check_data_type='table')
        return fetch_sequences(step_data, step, self.get_cache_object())
