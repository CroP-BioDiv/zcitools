import os.path
from step_project.base_commands import Command, CreateStepCommand

# Check base.py for description


class TableStep(CreateStepCommand):
    _COMMAND = 'table'
    _HELP = """
Creates table step. Mandatory argument is input data file.
Additional arguments specify how to interpret input data.
"""

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('filename', help='Input data')
        parser.add_argument('-f', '--format', help='Input data format (if needed). Values: excel, csv, txt')
        parser.add_argument('-c', '--columns', help='Columns. Format name1,type1:name2,type2:...')

    def run(self, step_data):
        from .input_file import create_table_step
        args = self.args
        return create_table_step(self.project, step_data, args.filename, data_format=args.format, columns=args.columns)


class SQLSelect(Command):
    _COMMAND = 'select'
    _HELP = "SQL select query on table steps"
    _PRESENTATION = False

    @staticmethod
    def set_arguments(parser):
        parser.add_argument(
            'steps', nargs='+', help='Table steps. From part. Tables are named a, b, c, ... in order of appearance.')
        parser.add_argument('-s', '--select', help='Select part. Format "[<table>.]column [[AS] name], *"')
        parser.add_argument('-w', '--where', help='Where part')
        parser.add_argument('-g', '--group-by', help='Group by part')
        parser.add_argument('-H', '--having', help='Having part')
        parser.add_argument('-o', '--order-by', help='Order part')
        #
        parser.add_argument('-r', '--result', default='print', help='Result into: print, excel, step')
        parser.add_argument('-f', '--output-filename', help='Export filename')

    def _is_step_cmd(self):
        return self.args.result[0].lower() == 's'

    def get_command_type(self):
        return 'new_step' if self._is_step_cmd() else None

    def prev_steps(self):
        return [os.path.normpath(p) for p in self.args.steps]

    def step_base_name(self):
        return self._COMMAND

    def run(self, step_data=None):
        from .select import select_data
        ps = self.args
        steps = [self.project.read_step(s, check_data_type=('table', 'table_grouped')) for s in ps.steps]
        select_data(ps.result, step_data, steps, ps.select,
                    where_part=ps.where, group_by_part=ps.group_by, having_part=ps.having, order_by_part=ps.order_by,
                    output_filename=ps.output_filename)
