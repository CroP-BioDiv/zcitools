"""
Command classes.
Class has to implement:
 - _COMMAND              : (class attribute) command argument to use
 - _HELP                 : (class attribute) command help text
 - set_arguments(parser) : (static method) sets command's arguments
 - run(args, step_name)  : (static method) runs command with given arguments.
                           Returns data that will be stored into step description
 - prev_steps(args)      : (static method) returns list of input steps
 - needs_editing(args)   : (static method) returns True if step data needs editing
"""


class _Command:
    @staticmethod
    def prev_steps(args):
        return []

    @staticmethod
    def needs_editing(args):
        return False


class _InitProject(_Command):
    _COMMAND = 'init'
    _HELP = "Initialize project in given direcotry name."

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('dirname', help='Directory name')
        parser.add_argument('-d', '--description', help='Project description text')

    @staticmethod
    def run(args, _):
        from .init_project import init_project
        init_project(args.dirname, args.description)


class _TableStep(_Command):
    _COMMAND = 'table'
    _HELP = """
Creates table step. Mandatory argument is input data file.
Additional arguments specify how to interpret input data.
"""

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('data', help='Input data')
        parser.add_argument('-f', '--format', help='Input data format (if needed). Values: excel, csv, txt')
        parser.add_argument('-c', '--convert-method', help='Set convert method')

    @staticmethod
    def run(args, step_name):
        from .steps.input_file import create_table_step
        return create_table_step(step_name, args.data, format=args.format, convert_method=args.convert_method)


commands_map = dict((cls._COMMAND, cls) for cls in locals().values() if hasattr(cls, '_COMMAND'))
