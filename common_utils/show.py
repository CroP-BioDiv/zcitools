import shutil
import itertools
import math
from .terminal_layout import StringColumns


def print_table(header, rows, sort=False, show_limit=None, ident=None):
    if sort:
        rows = sorted(rows)
    if show_limit and len(rows) > show_limit:
        n = (show_limit - 1) // 2
        rows = rows[:n] + [['...'] * len(rows[0])] + rows[-n:]
    if ident:
        rows = [ident + r for r in rows]
    print(StringColumns(rows, header=header))


def print_hierarchical_table(header, rows, show_limit=None):
    # Note: header consists of more rows
    if show_limit:
        first_order_rows = [i for i, r in enumerate(rows) if r[0]]
        if len(first_order_rows) > show_limit:
            h_depth = len(header)
            n = (show_limit - 1) // 2
            rows = rows[:first_order_rows[n]] + ([['...'] * len(rows[0])] * h_depth) + rows[first_order_rows[-n]:]
    print(StringColumns(rows, header=list(zip(*header))))


def print_ls_like_list(title, data, sort=False, width=None, min_rows_to_split=None):
    # Prints given list of strings in more columns
    # Data is list of strings.
    if sort:
        data = sorted(data)
    if not width:
        width = shutil.get_terminal_size((120, 20)).columns

    # Find number of data columns.
    num_data = len(data)
    if min_rows_to_split and num_data < min_rows_to_split:
        num_c = 1
    else:
        padding = 1
        lengths = [len(d) for d in data]
        max_l = max(lengths)
        min_l = min(lengths)
        if max_l + min_l + 4 * padding > width:
            num_c = 1
        else:
            average_l = math.ceil(sum(lengths) / num_data)
            num_c = math.floor(width / (average_l + 2 * padding))  # num_c >= 1
            if num_c > 1:
                def _can_fit(nc):
                    nr = math.ceil(num_data / nc)
                    needed_width = sum(max(lengths[i:i + nr]) for i in range(0, num_data, nr))
                    needed_width += nc * 2 * padding
                    return needed_width <= width

                if _can_fit(num_c):
                    # Try with more columns
                    for nc in range(num_c + 1, math.floor(width / (min_l + 2 * padding))):
                        if _can_fit(nc):
                            num_c = nc
                        else:
                            break
                else:
                    # Try with less columns
                    for num_c in range(num_c - 1, 1, -1):
                        if _can_fit(num_c):
                            break
                    else:
                        num_c = 1

    # Format rows
    if num_c > 1:
        num_rows = math.ceil(num_data / num_c)
        add_data = num_c * num_rows - num_data
        if add_data:
            data = data + [''] * add_data
        rows = [data[i::num_rows] for i in range(num_rows)]
    else:
        rows = [[d] for d in data]

    print(StringColumns(rows, header=[title] + [''] * (num_c - 1), padding=padding))
