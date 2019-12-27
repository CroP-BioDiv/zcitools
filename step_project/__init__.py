from .step_project import StepProject


# Note: this script is called from project main directory, all used filenames are relative to it!
def run(registered_commands, registered_steps):
    StepProject(registered_commands=registered_commands, registered_steps=registered_steps).run()
