from .general import *
from .table import *
from .phylogenetics import *

# commands_map = dict((cls._COMMAND, cls) for cls in locals().values() if getattr(cls, '_COMMAND', None))
registered_commands = [cls for cls in locals().values() if getattr(cls, '_COMMAND', None)]
