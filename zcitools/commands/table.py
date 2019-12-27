from zcitools.base.commands import CreateStepCommand

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
        from ..processing.input_file import create_table_step
        args = self.args
        return create_table_step(self.zcit, step_data, args.filename, data_format=args.format, columns=args.columns)
