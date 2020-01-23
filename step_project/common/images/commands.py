import os.path
from step_project.base_commands import CreateStepCommand

# Check base.py for description


class CircosCorrelation(CreateStepCommand):
    _COMMAND = 'circos_correlation'
    _HELP = """Creates image step from correlation data with Circos."""

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('-s', '--show-image', action='store_true', help='Display result image.')
        parser.add_argument('-i', '--input-filename', help='Input filename')
        parser.add_argument('-c', '--group-color', action='append', help='Group color. Format <group>,<color>')

        parser.add_argument('--one-width', default=1000, help='Width of correlation 1. Base scale factor!')
        parser.add_argument('--gap-correlations', default=40, help='Width of gap between neighbouring correlations')
        # parser.add_argument('-s', '--input-step', help='Input step')
        # parser.add_argument('-f', '--format', help='Input data format (if needed). Values: excel, csv, txt')
        # parser.add_argument('-c', '--columns', help='Columns. Format name1,type1:name2,type2:...')

    def run(self, step_data):
        from .circos_correlation import create_circos_correlation
        return create_circos_correlation(self.project, step_data, self.args)
