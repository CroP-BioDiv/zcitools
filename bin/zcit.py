#!/usr/bin/python3

import os.path
import sys
import argparse
from zcitools.utils.file_utils import write_yaml
from zcitools.command_classes import commands_map


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

if not command_obj.STEP_COMMAND:
    command_obj.run()
else:
    # Check is command run inside a project
    if not os.path.isfile('project_log.yml'):
        print('Error: script is not called on valid project!')
        sys.exit(0)
    #
    # Find step name. Format <num>_<step_base_name>[_<description>]
    prev_steps = command_obj.prev_steps()
    #
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

    # Store log data into project_log.yml
    step_data = dict(step_name=sn,
                     prev_steps=prev_steps,
                     cmd=' '.join(sys.argv[1:]))
    write_yaml([step_data], 'project_log.yml', mode='a')  # Appends yml list

    # Run command
    step_obj = command_obj.run(step_data)
