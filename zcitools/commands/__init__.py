from .general import *
from .annotations import *

commands_map = dict((cls._COMMAND, cls) for cls in locals().values() if getattr(cls, '_COMMAND', None))
