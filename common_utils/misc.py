import os
import re
import time
import datetime
import importlib
from functools import wraps
from .exceptions import ZCItoolsValueError


def find_registered(init_filename, package, dir_list=None):
    # Checks submodules __init__.py for attributes (lists) registered_commands and registered_steps
    # If attributes are not set, checks for classes in files commands.py and steps.py
    commands, steps = [], []
    _dir = os.path.dirname(init_filename)

    for f in (dir_list or os.listdir(_dir)):
        if not f.startswith('_') and os.path.isdir(os.path.join(_dir, f)):
            found_cmds = found_steps = False
            if os.path.isfile(os.path.join(_dir, f, '__init__.py')):
                r = importlib.import_module(f'.{f}', package=package)
                data = r.__dict__
                #
                found_cmds = 'registered_commands' in data
                if found_cmds:
                    commands.extend(data['registered_commands'])
                #
                found_steps = 'registered_steps' in data
                if found_steps:
                    steps.extend(data['registered_steps'])

            if not found_cmds:
                if os.path.isfile(os.path.join(_dir, f, 'commands.py')):
                    r = importlib.import_module(f'.{f}.commands', package=package)
                    commands.extend(c for c in r.__dict__.values() if getattr(c, '_COMMAND', None))
                else:
                    print(f"Warning: module {package}.{f} doesn't have registeded commands!")

            if not found_steps:
                if os.path.isfile(os.path.join(_dir, f, 'steps.py')):
                    r = importlib.import_module(f'.{f}.steps', package=package)
                    steps.extend(c for c in r.__dict__.values() if getattr(c, '_STEP_TYPE', None))
                else:
                    print(f"Warning: module {package}.{f} doesn't have registeded steps!")
    #
    return commands, steps


# Other
def split_list(data, num_items):
    assert isinstance(data, list), type(data)
    for i in range((len(data) // num_items) + 1):
        n = i * num_items
        yield data[n:(n + num_items)]


# Data checks
def sets_equal(have_to_exist, exist, description, step=None):
    # Is all data presented
    not_exist = have_to_exist - exist
    if not_exist:
        sd = f'({step}) ' if step else ''
        raise ZCItoolsValueError(f"{sd}Data for {description}(s) not presented: {', '.join(sorted(not_exist))}")

    # Is there more data than needed
    more_data = exist - have_to_exist
    if more_data:
        sd = f'({step}) ' if step else ''
        raise ZCItoolsValueError(f"{sd}Data exists for not listed {description}(s): {', '.join(sorted(more_data))}")


#
def YYYYMMDD_2_date(s):
    nums = re.findall(r'\d+', s)
    assert len(nums) == 3, s
    return datetime.date(*map(int, nums))


def time_it(f):
    @wraps(f)
    def wrap(*args, **kw):
        t = time.time()
        result = f(*args, **kw)
        t = 1000 * (time.time() - t)
        print(f'method:{f.__name__} took: {t:2.4f} ms')
        return result
    return wrap
