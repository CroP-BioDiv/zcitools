import os
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
 - db_identifier()       : returns None or dict(static=bool, data_identifier=list of strings)
                           Uniquelly specify how step data is created from project start.
                           Used to specify common DB data.

Create new step command also has to implement:
 - _PRESENTATION         : (attribute) bool describing is step presentation.
 - finish(step_object)   : finish previously created step. Used when editing is neeed.
 - prev_steps()          : returns list of input steps

Create step command's run() method can store additional data for finish() method.
It is good practice to store that data in file finish.yml.

"""

ZCI_COMMON_DB_DIR = os.environ.get('ZCI_COMMON_DB')
_SEQUENCE_DBS_RELATIVE_DIR = 'sequence_dbs'
_SEQUENCE_DBS_DIR = os.path.join(ZCI_COMMON_DB_DIR, _SEQUENCE_DBS_RELATIVE_DIR) if ZCI_COMMON_DB_DIR else None


class _Command:
    _PROJECT_COMMAND = True
    _COMMAND_TYPE = None  # General work
    _COMMAND = None
    _COMMAND_GROUP = None

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

    # Common DB
    def db_identifier(self):
        return None

    def get_common_db_object(self):
        if ZCI_COMMON_DB_DIR:
            d = self.db_identifier()
            if d:
                return CommonDB(os.path.join(ZCI_COMMON_DB_DIR, *d['data_identifier']), base_dir=ZCI_COMMON_DB_DIR)

    @staticmethod
    def get_sequence_dbs():
        if _SEQUENCE_DBS_DIR:
            return [f for f in os.listdir(_SEQUENCE_DBS_DIR) if os.path.isdir(os.path.join(_SEQUENCE_DBS_DIR, f))]
        return []

    def get_sequence_db_ident(self, *idents, db=None):
        db = db or self.args.database
        dbs = self.get_sequence_dbs()
        if db not in dbs:
            raise ZCItoolsValueError(f"Database {db} doesn't exist. Possible values: {', '.join(dbs)}!")
        return [_SEQUENCE_DBS_RELATIVE_DIR, db, *idents]


#
class ProjectCommand(_Command):
    pass


class NonProjectCommand(_Command):
    _PROJECT_COMMAND = False


# Create step mixins (interfaces)
class _CreateStepsMixin:
    _PRESENTATION = False
    _STEP_BASE_NAME = None

    def prev_steps(self):
        return [os.path.normpath(p) for p in self._prev_steps()]

    def _prev_steps(self):
        return []

    def finish(self, step_obj):
        pass

    def step_base_name(self):
        return self._STEP_BASE_NAME or self._COMMAND


class _CreateStepsFromStepMixin(_CreateStepsMixin):
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


class _CreateSteps_CommonDB_Mixin(_CreateStepsFromStepMixin):
    _COMMON_DB_IDENT = None  # 'sequences'

    @classmethod
    def set_arguments(cls, parser):
        parser.add_argument('step', help='Input table step')
        dbs = cls.get_sequence_dbs()
        parser.add_argument('-d', '--database', default='base', help=f'Database to use: {", ".join(dbs)}')

    def db_identifier(self):
        assert self._COMMON_DB_IDENT
        return dict(static=True, data_identifier=self.get_sequence_db_ident(self._COMMON_DB_IDENT))


# Create one step commands
class CreateStepCommand(_CreateStepsMixin, ProjectCommand):
    _COMMAND_TYPE = 'new_step'

    def run(self, step_data):
        raise NotImplementedError(f'Method {self.__class__.__name__}.run(step_data) is not implemented!')


class CreateStepFromStepCommand(_CreateStepsFromStepMixin, CreateStepCommand):
    pass


class CreateStepFromStepCommand_CommonDB(_CreateSteps_CommonDB_Mixin, CreateStepFromStepCommand):
    pass


# Create more steps commands
class CreateStepsCommand(_CreateStepsMixin, ProjectCommand):
    _COMMAND_TYPE = 'new_steps'

    def run(self, step_data, args):
        raise NotImplementedError(f'Method {self.__class__.__name__}.run(step_data) is not implemented!')


class CreateStepsFromStepCommand(_CreateStepsFromStepMixin, CreateStepsCommand):
    pass


class CreateStepsFromStepCommand_CommonDB(_CreateSteps_CommonDB_Mixin, CreateStepsFromStepCommand):
    pass
