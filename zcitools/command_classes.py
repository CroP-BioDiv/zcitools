"""
Command classes, called from zcit script.
There are few different types (groups) of commands, depending on needed data to run and return value.
Commands:
 - to create new step.
   Receive data of step to create, return None or step object.
   Step can be presentation step, which means that it a
 - to make presentation data of existing step.
   Receive step object, return bool if calculation finished.
   Note: calculation can be called more times. It should test calculation state and continue from where it stopped.
 - to make general work.
   No input arguments, return None.

Class has to implement:
 - __init__(args)        : constructor
 - _COMMAND              : (attribute) command argument to use
 - _HELP                 : (attribute) command help text
 - set_arguments(parser) : (static method) sets command's arguments
 - run(step_data/step/-) : runs command with given arguments

Create new step command also has to implement:
 - finish(step_object)   : finish previously created step. Used when editing is neeed.
 - prev_steps()          : returns list of input steps
 - cache_identifier()    : returns None or dict(static=bool, data_identifier=list of strings)
"""

# Note: importing is done in run() methods to prevent crashes because of not used missing libraries!
import os.path
from .steps import read_step
from .utils.exceptions import ZCItoolsValueError
from .commands_create_step import *
from .commands_presentation import *


class _Command:
    _COMMAND_TYPE = None  # General work
    _COMMAND = None

    def __init__(self, args):
        self.args = args

    def run(self):
        raise NotImplementedError(f'Method {self.__class__.__name__}.run() is not implemented!')


class _InitProject(_Command):
    _COMMAND = 'init'
    _HELP = "Initialize project in given directory name."

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
        step = read_step(self.args.step, update_mode=True)  # Set to be in update mode
        if step.get_step_needs_editing():
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
        step = read_step(self.args.step)
        step.show_data(params=self.args.params)


#
commands_map = dict((cls._COMMAND, cls) for cls in locals().values() if getattr(cls, '_COMMAND', None))
