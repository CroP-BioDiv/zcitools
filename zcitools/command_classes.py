"""
Command classes, called from zcit script.
There are few different types of commands, depending on needed data to run and return value.
 - commands to create new step
   receive data of step to create, return None or step object,
 - commnads to make calculation on existing step:
   receive step object, return bool if calculation was successfully
 - commnads to make general work
   no input arguments, return None

Class has to implement:
 - __init__(args)        : constructor
 - _COMMAND              : (attribute) command argument to use
 - _HELP                 : (attribute) command help text
 - set_arguments(parser) : (static method) sets command's arguments
 - run(step_data/step/-) : runs command with given arguments
 - finish(step_object)   : finish previously created step. Used when editing is neeed.
 - prev_steps()          : returns list of input steps
"""

# Note: importing is done in run() methods to prevent crashes because of not used missing libraries!
import os.path
from .steps import read_step
from .utils.exceptions import ZCItoolsValueError


class _Command:
    _COMMAND_TYPE = None  # General work

    def __init__(self, args):
        self.args = args

    def prev_steps(self):
        return [os.path.normpath(p) for p in self._prev_steps()]

    def _prev_steps(self):
        return []

    def run(self):
        raise NotImplementedError(f'Method {self.__class__.__name__}.run() is not implemented!')

    def finish(self, step_obj):
        pass


class _CreateStepCommand(_Command):
    _COMMAND_TYPE = 'new_step'
    _STEP_BASE_NAME = None

    def run(self, step_data):
        raise NotImplementedError(f'Method {self.__class__.__name__}.run(step_data) is not implemented!')

    def step_base_name(self):
        return self._STEP_BASE_NAME or self._COMMAND


class _CalculateCommand(_Command):
    _COMMAND_TYPE = 'calculate'
    _STEP_DATA_TYPE = None         # Type of step commnad works on
    _CALCULATION_DIRECTORY = None  # Name of subdirectory

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Step name')

    def run(self, step):
        raise NotImplementedError(f'Method {self.__class__.__name__}.run(step) is not implemented!')


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


class _Finish(_Command):
    _COMMAND = 'finish'
    _HELP = "Finish step that needed editing."

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Step name')

    def run(self):
        from .steps import read_step
        step = read_step(self.args.step, update_mode=True)  # Set to be in update mode
        if not step.get_step_needs_editing():
            command_obj = commands_map[step.get_step_command()](None)
            command_obj.finish(step)
        else:
            print("Info: step {self.args.step} is already finished!")


class _CleanCache(_Command):
    _COMMAND = 'cache'
    _HELP = "Remove cache of given steps"

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', nargs='+', help='Step name')

    def run(self):
        from .steps import read_step
        for s in self.args.step:
            step = read_step(s)
            step.remove_cache_files()


class _Show(_Command):
    _COMMAND = 'show'
    _HELP = "Print step(s) data"

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Step name')
        parser.add_argument('params', nargs='*', help='Additional format option (free format, depends on step type)')

    def run(self):
        from .steps import read_step
        step = read_step(self.args.step)
        step.show_data(params=self.args.params)


# --------------------------------------------------
# Calculation commands
# --------------------------------------------------
class _OGDRAW(_CalculateCommand):
    pass


# --------------------------------------------------
# Create step commands
# --------------------------------------------------
class _TableStep(_CreateStepCommand):
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
        from .processing.input_file import create_table_step
        args = self.args
        return create_table_step(step_data, args.filename, data_format=args.format, columns=args.columns)


class _NCBIStep(_CreateStepCommand):
    _COMMAND = 'ncbi'
    _HELP = "Creates sequences step. Mandatory argument is a table step."
    _STEP_BASE_NAME = 'NCBI'

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Input table step')
        parser.add_argument(
            '-f', '--force-download', action='store_true', help='Download even if step already contains data.')

    def _prev_steps(self):
        return [self.args.step]

    def run(self, step_data):
        from .processing.ncbi import download_ncbi
        step = read_step(self.args.step, check_data_type='table')
        return download_ncbi(step_data, step, force_download=self.args.force_download)


# Annotations
class _GeSeqStep(_CreateStepCommand):
    _COMMAND = 'ge_seq'
    _HELP = "Annotates chloroplast sequences with GeSeq"
    _STEP_BASE_NAME = 'GeSeq'

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Input sequences step')

    def _prev_steps(self):
        return [self.args.step]

    def run(self, step_data):
        from .processing.annotation.ge_seq import create_ge_seq_data
        step = read_step(self.args.step, check_data_type='sequences')
        return create_ge_seq_data(step_data, step)

    def finish(self, step_obj):
        from .processing.annotation.ge_seq import finish_ge_seq_data
        finish_ge_seq_data(step_obj)


class _GeSeqStep(_GeSeqStep):
    _COMMAND = 'cpgavas'
    _HELP = "Annotates chloroplast sequences with CPGAVAS"
    _STEP_BASE_NAME = 'CPGAVAS'

    def run(self, step_data):
        from .processing.annotation.cpgavas import create_cpgavas_data
        step = read_step(self.args.step, check_data_type='sequences')
        return create_cpgavas_data(step_data, step)

    def finish(self, step_obj):
        from .processing.annotation.ge_seq import finish_cpgavas_data
        finish_cpgavas_data(step_obj)


#
commands_map = dict((cls._COMMAND, cls) for cls in locals().values() if hasattr(cls, '_COMMAND'))
