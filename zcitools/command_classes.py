"""
Command classes. Class has to implement:
 - __init__(args)        : constructor
 - _COMMAND              : (attribute) command argument to use
 - _HELP                 : (attribute) command help text
 - set_arguments(parser) : (static method) sets command's arguments
 - run(step_data)        : runs command with given arguments. For create step method returns step object.
 - prev_steps()          : returns list of input steps
"""

# Note: importing is done in run() methods to prevent crashes because of not used missing libraries!
from .steps import read_step
from .utils.exceptions import ZCItoolsValueError


class _Command:
    STEP_COMMAND = False

    def __init__(self, args):
        self.args = args

    def prev_steps(self):
        return []

    def run(self):
        raise NotImplementedError(f'Method {self.__class__.__name__}.run() is not implemented!')


class _StepCommand(_Command):
    STEP_COMMAND = True
    _STEP_BASE_NAME = None

    def run(self, step_data):
        raise NotImplementedError(f'Method {self.__class__.__name__}.run(step_data) is not implemented!')

    def step_base_name(self):
        return self._STEP_BASE_NAME or self._COMMAND


# --------------------------------------------------
# Not step commands
# --------------------------------------------------
class _InitProject(_Command):
    _COMMAND = 'init'
    _HELP = "Initialize project in given direcotry name."

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('dirname', help='Directory name')
        parser.add_argument('-d', '--description', help='Project description text')

    def run(self):
        from .init_project import init_project
        init_project(self.args.dirname, self.args.description)


# --------------------------------------------------
# Step commands
# --------------------------------------------------
class _TableStep(_StepCommand):
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
        from .create_step.input_file import create_table_step
        args = self.args
        return create_table_step(step_data, args.filename, data_format=args.format, columns=args.columns)


class _NCBIStep(_StepCommand):
    _COMMAND = 'ncbi'
    _HELP = "Creates sequences step. Mandatory argument is a table step."

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Input table step')
        parser.add_argument(
            '-f', '--force-download', action='store_true', help='Download even if step already contains data.')

    def prev_steps(self):
        return [self.args.step]

    def run(self, step_data):
        from .create_step.ncbi import download_ncbi
        step = read_step(self.args.step, check_data_type='table')
        if not step:
            raise ZCItoolsValueError(f"Step {self.args.step} is not of data type 'table'!")
        return download_ncbi(step_data, step, force_download=self.args.force_download)


#
commands_map = dict((cls._COMMAND, cls) for cls in locals().values() if hasattr(cls, '_COMMAND'))
