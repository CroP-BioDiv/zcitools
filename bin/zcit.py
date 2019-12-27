#!/usr/bin/python3

from zcitools.zcit import ZCIT
from bio import registered_commands, registered_steps

# Note: this script is called from project main directory, all used filenames are relative to it!
zcit = ZCIT(registered_commands=registered_commands, registered_steps=registered_steps)
zcit.run()
