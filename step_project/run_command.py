import os.path
import sys
import argparse
from collections import defaultdict
from common_utils.file_utils import write_yaml, read_yaml
from common_utils.exceptions import ZCItoolsValueError
from .common import registered_commands as common_commands, registered_steps as common_steps

_zci_general_help = """Usage: python {exe} <command> <arguments>
Help: python {exe} help <command>

Commands:
{commands}"""


def _format_commands(commands_map, group, filter_m):
    cmd_groups = defaultdict(dict)
    for cmd, cls in commands_map.items():
        if filter_m(cls):
            cmd_groups[cls._COMMAND_GROUP or ''][cmd] = cls  # Override None with ''
    if not cmd_groups:
        return ''
    bb = '\n'
    ret = [f'{group} commands:']
    for c_g, cls in sorted(cmd_groups.items()):
        ret.append(f'  {c_g or "General"}:')
        max_l = max(len(c) for c in cls.keys())
        ret.extend(f"    {cmd.ljust(max_l)} {cl._HELP.strip().split(bb)[0]}" for cmd, cl in sorted(cls.items()))
    return '\n' + '\n'.join(ret) + '\n'


class RunCommand:
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
                commands = _format_commands(
                    self.commands_map, 'Non-project', lambda cl: not cl._PROJECT_COMMAND)
                commands += _format_commands(
                    self.commands_map, 'Project', lambda cl: cl._PROJECT_COMMAND and not cl._COMMAND_TYPE)
                commands += _format_commands(
                    self.commands_map, 'Step', lambda cl: cl._PROJECT_COMMAND and cl._COMMAND_TYPE)
                # commands = ', '.join(sorted(self.commands_map.keys()))

                print(_zci_general_help.format(exe=sys.argv[0], commands=commands))
            return

        if len(sys.argv) == 2 and sys.argv[1].lower() == 'list_commands':
            cmds = set(self.commands_map.keys())
            cmds.add('help')
            print(' '.join(sorted(cmds)))
            return

        command = sys.argv[1].lower()
        if command not in self.commands_map:
            print(f'Command "{command}" is not supported!')
            sys.exit(0)

        parser = self._get_parser(command, False)
        return self._run_command(command, parser.parse_args())

    def _run_command(self, command, args):
        self._args = args  # Store commands args
        command_obj = self.commands_map[command](self, args)
        command_type = command_obj.get_command_type()

        # General work
        if not command_type:
            if command_obj._PROJECT_COMMAND and not self._check_is_project_valid():
                return
            command_obj.run()

        # Create new step
        elif command_type in ('new_step', 'new_steps'):
            if not self._check_is_project_valid():
                return

            # Run command
            command_args = dict((k, v) for k, v in vars(args).items()
                                if k not in ('command', 'step_num', 'step_description'))
            db_id = command_obj.common_db_identifier()
            step_data = dict(prev_steps=command_obj.prev_steps(),
                             common_db_identifier=list(db_id) if db_id else None,
                             command=command,
                             command_args=command_args,
                             cmd=' '.join(sys.argv[1:]))
            ret = None
            if command_type == 'new_step':
                step_data['step_name'] = self.new_step_name(command_obj, args)
                ret = command_obj.run(step_data)
                if ret:
                    if not ret.is_completed():
                        print(f'Step is not finished, check instruction ({ret.directory}/INSTRUCTIONS.txt)!')
                else:
                    print("Warning: create step command didn't return step object!")
            else:
                ret = command_obj.run(step_data)
                if ret is not None:
                    for s in ret:
                        if not s.is_completed():
                            print(f'Step is not finished, check instruction ({s.directory}/INSTRUCTIONS.txt)!')
                else:
                    print("Warning: create steps command didn't return any step object!")

            if ret:
                # Store log data into project_log.yml
                step_data = dict((k, v) for k, v in step_data.items() if k in ('cmd', 'step_name'))
                # Do not store if step_data is equal as from last command?
                log = read_yaml('project_log.yml')
                if not log or log[-1] != step_data:
                    write_yaml([step_data], 'project_log.yml', mode='a')  # Appends yml list

        else:
            print(f"Warning: not supported command_type {command_type}?!")

    def _get_parser(self, command, for_help):
        command_cls = self.commands_map[command]

        parser = argparse.ArgumentParser(
            description=command_cls._HELP,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        if for_help:
            parser.add_argument(command, help=command)
        else:
            parser.add_argument('command', help=command)

        if command_cls._COMMAND_TYPE in ('new_step', 'new_steps'):
            parser.add_argument('-N', '--step-num', type=int, help='Step num prefix')
            parser.add_argument('-D', '--step-description', help='Step description')
            parser.add_argument('--no-data-check', action='store_true', help=f'Do not check step data on loading.')
        command_cls.set_arguments(parser)
        return parser

    def _check_is_project_valid(self):
        # Check is command run inside a project
        if os.path.isfile('project_log.yml'):
            return True
        print('Error: script is not called on valid project!')
        return False

    def new_step_name(self, command_obj, args):
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
        # elif prev_steps:  # ??? desc by previous step?
        #     # Note: for now commands do not have '_'
        #     for s in sorted(prev_steps, reverse=True):
        #         d = s.split('_')[2:]
        #         if d:
        #             desc = '_'.join(d)
        #             break
        #
        sn = f'{num:02}_{command_obj.step_base_name()}'
        if desc:
            sn += f'_{desc}'

        return sn

    # Read step method
    def read_step(self, step_name, check_data_type=None, update_mode=False, no_check=False):
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

        return cls(self, desc_data['project'], update_mode=update_mode, no_check=no_check)

    def new_step(self, cls, step_data, remove_data=False, update_mode=False, no_check=False):
        return cls(self, step_data, remove_data=remove_data, update_mode=update_mode,
                   no_check=no_check or self._args.no_data_check)

    def new_step_by_type(self, data_type, step_data, remove_data=False, update_mode=False, no_check=False):
        cls = self.steps_map.get(data_type)
        if not cls:
            raise ZCItoolsValueError(f"No step class for data type {data_type}!")
        return cls(self, step_data, remove_data=remove_data, update_mode=update_mode,
                   no_check=no_check or self._args.no_data_check)

    def find_previous_step_of_type(self, step, prev_step_type):
        # ToDo:
        pass
