from types import SimpleNamespace
from .general import *
from .table import *
from .sequences import *
from .annotations import *
from .alignments import *
from .phylogenetics import *


# Needs commands_map!
class Finish:  # (Command):
    _COMMAND_TYPE = None  # General work
    _COMMAND = 'finish'
    _HELP = "Finish step that needed editing."

    def __init__(self, args):
        self.args = args

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Step name')

    def run(self):
        from ..steps import read_step
        step = read_step(self.args.step, update_mode=True)  # Set to be in update mode
        if step.is_completed():
            print(f"Info: step {self.args.step} is completed!")
        else:
            orig_args = SimpleNamespace(**step._step_data['command_args'])
            command_obj = commands_map[step.get_command()](orig_args)
            command_obj.finish(step)


commands_map = dict((cls._COMMAND, cls) for cls in locals().values() if getattr(cls, '_COMMAND', None))
