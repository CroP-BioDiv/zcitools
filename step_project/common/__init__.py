from .general_commands import *
from common_utils.misc import find_registered

general_commands = [c for c in locals().values() if getattr(c, '_COMMAND', None)]
registered_commands, registered_steps = find_registered(__file__, 'step_project.common')
registered_commands += general_commands
