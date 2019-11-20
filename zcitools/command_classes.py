"""
Command classes.
Class has to implement:
 - __init__(args)        : constructor
 - _COMMAND              : (attribute) command argument to use
 - _HELP                 : (attribute) command help text
 - set_arguments(parser) : (static method) sets command's arguments
 - run(step_name)        : runs command with given arguments.
                           Returns data that will be stored into step description
 - prev_steps()          : returns list of input steps
 - needs_editing()       : returns True if step data needs editing
"""


class _Command:
    def __init__(self, args):
        self.args = args

    def prev_steps(self):
        return []

    def needs_editing(self):
        return False


class _InitProject(_Command):
    _COMMAND = 'init'
    _HELP = "Initialize project in given direcotry name."

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('dirname', help='Directory name')
        parser.add_argument('-d', '--description', help='Project description text')

    def run(self, _):
        from .init_project import init_project
        init_project(self.args.dirname, self.args.description)


class _TableStep(_Command):
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
        # parser.add_argument('-m', '--convert-method', help='Set convert method')

    def run(self, step_data):
        from .create_step.input_file import create_table_step
        args = self.args
        return create_table_step(step_data, args.filename, data_format=args.format, columns=args.columns)


commands_map = dict((cls._COMMAND, cls) for cls in locals().values() if hasattr(cls, '_COMMAND'))
