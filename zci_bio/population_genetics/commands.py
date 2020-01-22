from step_project.base_commands import CreateStepCommand
from common_utils.exceptions import ZCItoolsValueError


class NewHybrids(CreateStepCommand):
    _COMMAND = 'new_hybrids'
    _HELP = "Run NewHybrids on given data"
    _STEP_BASE_NAME = 'NewHybrids'
    # Note: no global caching

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('data_file', help='Input data file')
        parser.add_argument('gtyp_cat_file', help='Genotype category probabilities')
        parser.add_argument('-t', '--theta-prior', default='Jeffreys', help='Theta prior. Jeffrey or Uniform')
        parser.add_argument('-p', '--pi-prior', default='Jeffreys', help='Pi prior. Jeffrey or Uniform')
        parser.add_argument('-b', '--burn-in', default=5, type=int, help='Number of burn-in sweeps')
        parser.add_argument('-s', '--num-sweeps', default=20, type=int, help='Number of sweeps')
        parser.add_argument('-n', '--num-runs', default=1, type=int, help='Number of runs')
        parser.add_argument('-r', '--run', action='store_true', help='Run NewHybrid locale')

    def run(self, step_data):
        from .new_hybrids import create_new_hybrids_data
        return create_new_hybrids_data(self.project, step_data, self.args)

    def finish(self, step_obj):
        from .new_hybrids import finish_new_hybrids
        finish_new_hybrids(step_obj)
