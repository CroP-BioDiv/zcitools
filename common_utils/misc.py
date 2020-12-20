import os
import re
import time
import datetime
import importlib
from functools import wraps
from .exceptions import ZCItoolsValueError


def find_registered(init_filename, package, dir_list=None, exclude_dirs=None):
    # Checks submodules __init__.py for attributes (lists) registered_commands and registered_steps
    # If attributes are not set, checks for classes in files commands.py and steps.py
    commands, steps = [], []
    _dir = os.path.dirname(init_filename)
    exclude_dirs = exclude_dirs or []

    for f in (dir_list or os.listdir(_dir)):
        if f in exclude_dirs:
            continue
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


#
_char_to_ext = dict(T=12, G=9, M=6, K=3)
_chat_to_test = sorted(_char_to_ext.items(), reverse=True, key=lambda x: x[1])


def int_2_human(n):
    for c, e in _chat_to_test:
        if n >= 10**e:
            if n >= 10**(e + 2):  # 3 digits
                return f"{n // 10**e}{c}"
            return f'{round(n / 10**e, 1)}{c}'
    return str(n)


def human_2_int(s):
    if isinstance(s, int):
        return s
    if s.isdigit():
        return int(s)
    e = _char_to_ext.get(s[-1].upper())
    if e:
        return round(float(s[:-1]) * 10**e)
    raise ZCItoolsValueError(f"String {s} is not an integer!")


def coverage_2_human(x):
    if x:
        if x >= 100:
            return f'{round(x)}x'
        return f'{round(x, 1)}x'.replace('.0', '')
    return x


def time_it(f):
    @wraps(f)
    def wrap(*args, **kw):
        t = time.time()
        result = f(*args, **kw)
        t = 1000 * (time.time() - t)
        print(f'method:{f.__name__} took: {t:2.4f} ms')
        return result
    return wrap
