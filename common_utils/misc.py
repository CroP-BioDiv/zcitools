from .exceptions import ZCItoolsValueError


def find_registered():
    return [], []


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
