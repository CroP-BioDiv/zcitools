from ..utils.terminal_layout import StringColumns


def print_table(header, rows, sort=False, show_limit=None):
    if sort:
        rows = sorted(rows)
    if show_limit and len(rows) > show_limit:
        n = (show_limit - 1) // 2
        rows = rows[:n] + [['...'] * len(rows[0])] + rows[-n:]
    print(StringColumns(rows, header=header))
