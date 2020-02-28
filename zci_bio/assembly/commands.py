from step_project.base_commands import CreateStepCommand
from common_utils.exceptions import ZCItoolsValueError


class NOVOPlastyStep(CreateStepCommand):
    _COMMAND = 'novoplasty'
    _HELP = "Assembly with NOVOPlasty assembler"
    _STEP_BASE_NAME = 'NOVOPlasty'

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('organelle_type', help='Organelle type (mito, chloro). Sets default genome range')
        parser.add_argument('-r', '--sequence-reads', help='Specify sequence reads file')
        parser.add_argument('-s', '--seed', help='Seed Input')
        # ToDo: sequences other way
        parser.add_argument('-p', '--project-name', help='Project name')
        parser.add_argument('-g', '--genome-range-from', type=int, help='Genome range from')
        parser.add_argument('-G', '--genome-range-to', type=int, help='Genome range to')
        parser.add_argument('-k', '--k-mer', default=39, help='K-mer')

    def run(self, step_data):
        from .novoplasty import novoplasty_run
        return novoplasty_run(self.project, step_data, self.args)
