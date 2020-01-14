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


class TableSelect(Command):
    _COMMAND = 'table_select'
    _HELP = "SQL select like query"

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Main table step name')
        parser.add_argument('-o', '--output', help='Output filename')
        parser.add_argument('-s', '--select', help='Select part')
        # Format "column_1 table_step column_2". Means: relation.column_1 = table_step column_2
        parser.add_argument('-j', '--join', nargs=3, action='append', help='From part (join)')
        parser.add_argument('-w', '--where', help='Where part')

    def run(self):
        params = self.args
        step = self.project.read_step(params.step)
        relation = step.get_relation()
        #
        if params.join:
            for c1, step, c2 in join:
                s = self.project.read_step(step)
                relation = relation.join(c1, s.get_relation(), c2)
        #
        if params.where:
            relation = relation.where(params.where)
        #
        if params.select:
            relation = relation.select([s.strip() for s in params.select.split(',')])
        #
        if params.output:
            relation.to_excel(params.output)
        else:
            relation.print_data()
