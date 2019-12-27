import os.path
import sys
import argparse
from common_utils.file_utils import write_yaml, read_yaml
from common_utils.exceptions import ZCItoolsValueError
from .commands import registered_commands as common_commands
from .steps import registered_steps as common_steps


# Note: this script is called from project main directory, all used filenames are relative to it!
def run(registered_commands, registered_steps):
    ZCIT(registered_commands=registered_commands, registered_steps=registered_steps).run()


class ZCIT:
    def __init__(self, registered_commands=None, registered_steps=None):
        self.commands_map = dict()  # command -> command class that processes it
        self.steps_map = dict()
        # Add common commands and step types
        self.register_commands(common_commands)
        self.register_steps(common_steps)
        #
        if registered_commands:
            self.register_commands(registered_commands)
        if registered_steps:
            self.register_steps(registered_steps)

    def register_commands(self, command_classes):
        for command_cls in command_classes:
            cmd = command_cls._COMMAND
            if cmd in self.commands_map:
                raise ZCItoolsValueError(
                    f"Command {cmd} already registerd ({self.commands_map[cmd].__class__.__name__})!")
            self.commands_map[cmd] = command_cls

    def register_steps(self, step_classes):
        for step_cls in step_classes:
            s_type = step_cls._STEP_TYPE
            if s_type in self.steps_map:
                raise ZCItoolsValueError(
                    f"Step {s_type} already registerd ({self.steps_map[s_type].__class__.__name__})!")
            self.steps_map[s_type] = step_cls

    #
    def run(self):
        # Read command line parameters
        if len(sys.argv) == 1 or sys.argv[1].lower() == 'help':
            if len(sys.argv) > 2 and sys.argv[2].lower() in self.commands_map:
                command = sys.argv[2].lower()
                parser = self._get_parser(command, True)
                parser.print_help()
            else:
                print(f"""Usage: python {sys.argv[0]} <command> <arguments>
        Help: python {sys.argv[0]} help <command>

        Command is one of: {', '.join(sorted(self.commands_map.keys()))}""")
            return

        command = sys.argv[1].lower()
        if command not in self.commands_map:
            print(f'Command "{command}" is not supported!')
            sys.exit(0)

        parser = self._get_parser(command, False)
        args = parser.parse_args()
        command_obj = self.commands_map[command](self, args)

        # General work
        if not command_obj._COMMAND_TYPE:
            if command != 'init':
                if not self._check_is_project_valid():
                    return
            command_obj.run()

        # Create new step
        elif command_obj._COMMAND_TYPE == 'new_step':
            if not self._check_is_project_valid():
                return

            # Run command
            command_args = dict((k, v) for k, v in vars(args).items()
                                if k not in ('command', 'step_num', 'step_description'))
            step_data = dict(step_name=self._new_step_name(command_obj, args),
                             prev_steps=command_obj.prev_steps(),
                             command=command,
                             command_args=command_args,
                             cmd=' '.join(sys.argv[1:]))
            step_obj = command_obj.run(step_data)

            if step_obj:
                # Store log data into project_log.yml
                step_data = dict((k, v) for k, v in step_data.items() if k in ('cmd', 'step_name'))
                # Do not store if step_data is equal as from last command?
                log = read_yaml('project_log.yml')
                if not log or log[-1] != step_data:
                    write_yaml([step_data], 'project_log.yml', mode='a')  # Appends yml list
            else:
                print("Warning: create step command didn't return step object!")

    def _get_parser(self, command, for_help):
        command_cls = self.commands_map[command]

        parser = argparse.ArgumentParser(description=command_cls._HELP)
        if for_help:
            parser.add_argument(command, help=command)
        else:
            parser.add_argument('command', help=command)
        if command != 'init':
            parser.add_argument('-N', '--step-num', type=int, help='Step num prefix')
            parser.add_argument('-D', '--step-description', help='Step description')
        command_cls.set_arguments(parser)
        return parser

    def _check_is_project_valid(self):
        # Check is command run inside a project
        if os.path.isfile('project_log.yml'):
            return True
        print('Error: script is not called on valid project!')
        return False

    def _new_step_name(self, command_obj, args):
        # Find step name. Format <num>_<step_base_name>[_<description>]
        prev_steps = command_obj.prev_steps()

        if command_obj._PRESENTATION:
            # For presentation just add description
            assert prev_steps
            max_st = max(prev_steps)
            return f'{max_st}-{command_obj.step_base_name()}'

        if args.step_num:
            num = args.step_num
        elif prev_steps:
            num = max(int(s.split('_', 1)[0]) for s in prev_steps) + 1
        else:
            num = 1
        #
        desc = None
        if args.step_description:
            desc = args.step_description
        elif prev_steps:
            # Note: for now commands do not have '_'
            for s in sorted(prev_steps, reverse=True):
                d = s.split('_')[2:]
                if d:
                    desc = '_'.join(d)
                    break
        #
        sn = f'{num:02}_{command_obj.step_base_name()}'
        if desc:
            sn += f'_{desc}'

        return sn

    # Read step method
    def read_step(self, step_name, check_data_type=None, update_mode=False):
        if isinstance(step_name, str):
            desc_data = read_yaml(os.path.join(step_name, 'description.yml'))
        else:
            assert isinstance(step_name, list), type(step_name)
            desc_data = read_yaml(os.path.join(*step_name, 'description.yml'))
        if not desc_data:
            raise ZCItoolsValueError(f"'{step_name}' is not a step!")

        data_type = desc_data['data_type']

        if check_data_type:
            if isinstance(check_data_type, str):
                if check_data_type != data_type:
                    raise ZCItoolsValueError(f"Step {step_name} is not of data type '{check_data_type}'!")
            else:
                if data_type not in check_data_type:
                    raise ZCItoolsValueError(f"Step {step_name} is not of data types: {', '.join(check_data_type)}!")

        cls = self.steps_map.get(data_type)
        if not cls:
            raise ZCItoolsValueError(f"No step class for data type {data_type}!")

        return cls(self, desc_data['project'], update_mode=update_mode)
