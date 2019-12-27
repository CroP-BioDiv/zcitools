from .general import *
from .table import *
from .sequences import *
from .annotations import *
from .alignments import *
from .phylogenetics import *

# commands_map = dict((cls._COMMAND, cls) for cls in locals().values() if getattr(cls, '_COMMAND', None))
registered_commands = [cls for cls in locals().values() if getattr(cls, '_COMMAND', None)]
