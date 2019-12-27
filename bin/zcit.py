#!/usr/bin/python3

from step_project.zcit import run
from zci_bio import registered_commands, registered_steps

# Note: this script is called from project main directory, all used filenames are relative to it!
run(registered_commands, registered_steps)
