from .general import *
from .table import *

registered_commands = [cls for cls in locals().values() if getattr(cls, '_COMMAND', None)]
