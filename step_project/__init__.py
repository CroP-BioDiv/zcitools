from .run_command import RunCommand


# Note: this script is called from project main directory, all used filenames are relative to it!
def run(registered_commands, registered_steps, registered_workflows):
    RunCommand(registered_commands=registered_commands,
               registered_steps=registered_steps,
               registered_workflows=registered_workflows).run()
