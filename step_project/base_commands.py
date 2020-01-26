import os
from common_utils.cache import Cache

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
 - cache_identifier()    : returns None or dict(static=bool, data_identifier=list of strings)
                           Uniquelly specify how step data is created from project start.
                           Used for caching data.

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
    _CACHE_DIR_PROJECT = '_project_cache_'

    def __init__(self, project, args):
        self.project = project
        self.args = args

    def get_command_type(self):
        return self._COMMAND_TYPE

    @staticmethod
    def set_arguments(parser):
        raise NotImplementedError(f'Method {self.__class__.__name__}.set_arguments(parser) is not implemented!')

    def run(self):
        raise NotImplementedError(f'Method {self.__class__.__name__}.run() is not implemented!')

    # Caching
    def cache_identifier(self):
        return None

    def get_cache_object(self):
        d = self.cache_identifier()
        if d:
            if d['static']:
                # Default is not so global cache :-)
                _dir = os.environ['ZCI_GLOBAL_CACHE'] or os.path.join('..', '..', '_global_cache_')
            else:
                _dir = self._CACHE_DIR_PROJECT
            if _dir:
                return Cache(os.path.join(_dir, *d['data_identifier']))


class ProjectCommand(_Command):
    pass


class NonProjectCommand(_Command):
    _PROJECT_COMMAND = False


class CreateStepCommand(ProjectCommand):
    _COMMAND_TYPE = 'new_step'
    _PRESENTATION = False
    _STEP_BASE_NAME = None

    def prev_steps(self):
        return [os.path.normpath(p) for p in self._prev_steps()]

    def _prev_steps(self):
        return []

    def run(self, step_data):
        raise NotImplementedError(f'Method {self.__class__.__name__}.run(step_data) is not implemented!')

    def step_base_name(self):
        return self._STEP_BASE_NAME or self._COMMAND

    def finish(self, step_obj):
        pass


class CreateStepFromStepCommand(CreateStepCommand):
    # Same as above, but assumes one step as input parameter
    _INPUT_STEP_DATA_TYPE = None

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Input sequences step')

    def _prev_steps(self):
        return [self.args.step]

    def _input_step(self):
        assert self._INPUT_STEP_DATA_TYPE
        return self.project.read_step(self.args.step, check_data_type=self._INPUT_STEP_DATA_TYPE)
