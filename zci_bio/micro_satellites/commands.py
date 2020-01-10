from step_project.base_commands import CreateStepCommand
from common_utils.exceptions import ZCItoolsValueError


class QTLCartPermutation(CreateStepCommand):
    _COMMAND = 'qtl_cart_perm'
    _HELP = "Run QTL cartographer's permutation"
    _STEP_BASE_NAME = 'QTLCart'
    # Note: no global caching

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('raw', help='Input MapMaker raw file')
        parser.add_argument('-p', '--permutations', default=1000, type=int, help='Number of permutations')
        parser.add_argument(
            '-t', '--num-traits', type=int, help='Number of traits. Calculated from input files if not set.')
        parser.add_argument('-r', '--run', action='store_true', help='Run QTLCart locale')

    def run(self, step_data):
        from .qtl_cartographer import create_permutations
        a = self.args
        return create_permutations(self.project, step_data, a.raw, a.permutations, num_traits=a.num_traits, run=a.run)

    def finish(self, step_obj):
        from .qtl_cartographer import finish_permutations
        finish_permutations(step_obj)
