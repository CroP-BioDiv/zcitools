#!/usr/bin/python3

from step_project import run
from zci_bio import registered_commands, registered_steps, registered_workflows

# Note: this script is called from project main directory, all used filenames are relative to it!
run(registered_commands, registered_steps, registered_workflows)
