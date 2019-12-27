import os
import importlib
from .exceptions import ZCItoolsValueError


def find_registered(init_filename, package):
    commands, steps = [], []
    _dir = os.path.dirname(init_filename)
    for f in os.listdir(_dir):
        if not f.startswith('_') and os.path.isdir(os.path.join(_dir, f)):
            r = importlib.import_module(f'.{f}', package=package)
            data = r.__dict__
            #
            if 'registered_commands' in data:
                commands.extend(data['registered_commands'])
            else:
                print(f"Warning: module {package}.{f} doesn't have registeded commands!")
            #
            if 'registered_steps' in data:
                steps.extend(data['registered_steps'])
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
