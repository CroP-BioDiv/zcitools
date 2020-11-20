import os.path
from common_utils.common_db import CommonDB

"""
Command classes, called from zcit script.
There are two different types of commands, depending on needed data to run and return value.
Commands:
 - to make general work.
   No input arguments, return None.

 - to create new step.
   Receive data of step to create, return None or step object.

   Step can be presentation step, which means it's data is addition to existing step data and probably will
   not be used for further calculations.
   Used only for naming that step. Presentation step is named same as base step with additional description.

Class has to implement:
 - __init__(args)        : constructor
 - _COMMAND              : (attribute) command argument to use
 - _HELP                 : (attribute) command help text
 - set_arguments(parser) : (static method) sets command's arguments
 - run(step_data/step/-) : runs command with given arguments
 - common_db_identifier(): returns None or data_identifier (tuple of strings)
                           Uniquelly specify how step data is created from project start.
                           Used to specify common DB data.

Create new step command also has to implement:
 - _PRESENTATION         : (attribute) bool describing is step presentation.
 - finish(step_object)   : finish previously created step. Used when editing is neeed.
 - prev_steps()          : returns list of input steps

Create step command's run() method can store additional data for finish() method.
It is good practice to store that data in file finish.yml.

"""


class _Command:
    _PROJECT_COMMAND = True
    _COMMAND_TYPE = None  # General work
    _COMMAND = None
    _COMMAND_GROUP = None
    _COMMON_DB_IDENT = None

    def __init__(self, project, args):
        self.project = project
        self.args = args

    def get_command_type(self):
        return self._COMMAND_TYPE

    @staticmethod
    def set_arguments(parser):
        pass

    def run(self):
        raise NotImplementedError(f'Method {self.__class__.__name__}.run() is not implemented!')

    # Common DB
    def common_db_identifier(self):
        return self._COMMON_DB_IDENT

    def get_common_db_object(self):
        return self._common_db_object(self)

    def get_step_db_object(self, step):
        return self._common_db_object(step)

    def get_idents_common_db_object(self, idents):
        assert (isinstance(idents, tuple) and all(isinstance(i, str) for i in idents)), idents
        return CommonDB.get_zci_db(idents)

    def _common_db_object(self, obj):
        idents = obj.common_db_identifier()
        assert idents is None or (isinstance(idents, tuple) and all(isinstance(i, str) for i in idents)), idents
        if idents:
            return CommonDB.get_zci_db(idents)


#
class ProjectCommand(_Command):
    pass


class NonProjectCommand(_Command):
    _PROJECT_COMMAND = False


# Create one step commands
class CreateStepCommand(ProjectCommand):
    _COMMAND_TYPE = 'new_step'
    _PRESENTATION = False
    _STEP_BASE_NAME = None

    def step_base_name(self):
        return self._STEP_BASE_NAME or self._COMMAND

    def prev_steps(self):
        return [os.path.normpath(p) for p in self._prev_steps()]

    def _prev_steps(self):
        return []

    def run(self, step_data):
        raise NotImplementedError(f'Method {self.__class__.__name__}.run(step_data) is not implemented!')

    def finish(self, step_obj):
        pass


class CreateStepFromStepCommand(CreateStepCommand):
    # Same as above, but assumes one step as input parameter
    _INPUT_STEP_DATA_TYPE = None

    def step_base_name(self):
        return self._format_step_name(super().step_base_name())

    def _format_step_name(self, name):
        step = self._input_step(no_data_check=True)
        name_prefix = step.get_step_name_prefix()
        if self.args.append_input_step_name:
            if sn := step.get_base_step_name():
                name = f'{sn}_{name}'  # Clearer naming with prefix, than sufix
        if name_prefix:
            name = f'{name_prefix}_{name}'
        return name

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Input sequences step')
        parser.add_argument('-A', '--append-input-step-name', action='store_true',
                            help='Append input sequence step name')

    def common_db_identifier(self):
        return self._COMMON_DB_IDENT or self._input_step(no_data_check=True).common_db_identifier()

    def _prev_steps(self):
        return [self.args.step]

    def _input_step(self, no_data_check=False):
        assert self._INPUT_STEP_DATA_TYPE
        a = self.args
        return self.project.read_step(
            a.step, check_data_type=self._INPUT_STEP_DATA_TYPE, no_check=(no_data_check or a.no_data_check))


# Create more steps commands
class CreateStepsCommand(ProjectCommand):
    _COMMAND_TYPE = 'new_steps'
    _PRESENTATION = False
    _STEP_BASE_NAME = None

    def step_base_name(self):
        return self._STEP_BASE_NAME or self._COMMAND

    def prev_steps(self):
        return [os.path.normpath(p) for p in self._prev_steps()]

    def _prev_steps(self):
        return []

    def run(self, step_data):
        raise NotImplementedError(f'Method {self.__class__.__name__}.run(step_data) is not implemented!')


class CreateStepsFromStepCommand(CreateStepsCommand):
    _INPUT_STEP_DATA_TYPE = None

    @staticmethod
    def set_arguments(parser):
        CreateStepsCommand.set_arguments(parser)
        parser.add_argument('step', help='Input sequences step')

    def common_db_identifier(self):
        return self._COMMON_DB_IDENT or self._input_step(no_data_check=True).common_db_identifier()

    def _prev_steps(self):
        return [self.args.step]

    def _input_step(self, no_data_check=False):
        assert self._INPUT_STEP_DATA_TYPE
        a = self.args
        return self.project.read_step(
            a.step, check_data_type=self._INPUT_STEP_DATA_TYPE, no_check=(no_data_check or a.no_data_check))
