#!/usr/bin/python3

import argparse


# Supported commnads
def _init_project(*args, **kwargs):
    if args:
        from zcitools.init_project import init_project
        init_project(args[0], *args[1:])
    else:
        print('Directory is not set!')


_commands = dict(
    init=_init_project,
)


# Handling of commnad line params
def _format_k(a):
    assert a[0].isalpha()
    a = a.lower()
    a = re.sub('[^0-9a-zA-Z]+', '_', a)
    return a


def _format_kwarg(a):
    if a.startswith('--'):
        return _format_k(a[2:])
    if a.startswith('-'):
        return _format_k(a[1:])
    return _format_k(a)


def _format_value(v):
    try:
        return int(v)
    except ValueError:
        pass
    try:
        return float(v)
    except ValueError:
        pass
    return v


#
parser = argparse.ArgumentParser(description='Runs zcitools command')

parser.add_argument('command', help='Command, one of: ' + ', '.join(sorted(_commands.keys())))
parser.add_argument('args', nargs='*', help='Arguments specific to command')

args, unknownargs = parser.parse_known_args()

if args.command not in _commands:
    import sys
    print(f'Command {args.command} is not supported!')
    sys.exit(0)

# Process unknownargs
all_args = list(args.args)
kwargs = dict()
in_kwarg = None
for a in unknownargs:
    if is_kwarg(a):
        if in_kwarg:
            kwargs[in_kwarg] = True
        in_kwarg = _format_kwarg(a)
    else:
        a = _format_value(a)
        if in_kwarg:
            kwargs[in_kwarg] = a
            in_kwarg = None
        else:
            all_args.append(a)
if in_kwarg:
    kwargs[in_kwarg] = True

_commands[args.command](*all_args, **kwargs)
