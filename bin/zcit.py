#!/usr/bin/python3

import os.path
import sys
import argparse
from zcitools.utils.file_utils import write_yaml, read_yaml, remove_directory, ensure_directory
from zcitools.commands import commands_map
from zcitools.steps import read_step

# Note: this script is called from project main directory, all used filenames are relative to it!


def _get_parser(command, for_help):
    cls = commands_map[command]

    parser = argparse.ArgumentParser(description=cls._HELP)
    if for_help:
        parser.add_argument(command, help=command)
    else:
        parser.add_argument('command', help=command)
    if command != 'init':
        parser.add_argument('-N', '--step-num', type=int, help='Step num prefix')
        parser.add_argument('-D', '--step-description', help='Step description')
    cls.set_arguments(parser)
    return parser


def _check_is_project_valid():
    # Check is command run inside a project
    if not os.path.isfile('project_log.yml'):
        print('Error: script is not called on valid project!')
        sys.exit(0)


def _new_step_name(command_obj, args):
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


# ---------------------------------------------------------
if len(sys.argv) == 1 or sys.argv[1].lower() == 'help':
    if len(sys.argv) > 2 and sys.argv[2] in commands_map:
        parser = _get_parser(sys.argv[2], True)
        parser.print_help()
    else:
        print(f"""Usage: python {sys.argv[0]} <command> <arguments>
Help: python {sys.argv[0]} help <command>

Command is one of: {', '.join(sorted(commands_map.keys()))}""")
    sys.exit(0)

command = sys.argv[1].lower()
if command not in commands_map:
    print(f'Command "{command}" is not supported!')
    sys.exit(0)

parser = _get_parser(command, False)
args = parser.parse_args()
command_obj = commands_map[command](args)


# General work
if not command_obj._COMMAND_TYPE:
    if command != 'init':
        _check_is_project_valid()
    command_obj.run()

# Create new step
elif command_obj._COMMAND_TYPE == 'new_step':
    _check_is_project_valid()

    # Run command
    step_data = dict(step_name=_new_step_name(command_obj, args),
                     prev_steps=command_obj.prev_steps(),
                     cache=command_obj.cache_identifier(),
                     command=command,
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
