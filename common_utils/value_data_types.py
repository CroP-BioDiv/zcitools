from datetime import date
from decimal import Decimal


KNOWN_DATA_TYPES = frozenset(['seq_ident', 'str', 'int', 'date', 'decimal'])


def _ident_cast(x):
    return x


def fromisoformat(d):
    # datetime.date.fromisoformat(data['Date']),  # In python 3.7 :-/
    return date(*map(int, d.split('-')))


# Note: input value is string or None
data_type_cast = dict(
    seq_ident=_ident_cast,
    str=_ident_cast,
    int=lambda x: int(x) if x else None,
    date=lambda x: fromisoformat(x) if x else None,
    decimal=lambda x: Decimal(x) if x else None,
)

needs_casting = dict(
    seq_ident=False,
    str=False,
    int=True,
    date=True,
    decimal=True,
)

data_type_python = dict(
    seq_ident=str,
    str=str,
    int=int,
    date=date,
    decimal=Decimal,
)


# Columns are in format [(name, data_type), ...]
def columns_needs_casting(data_types):
    return any(needs_casting[dt] for dt in data_types)


def columns_cast(data_types):
    return [data_type_cast[dt] for dt in data_types]


def cast_table(data_types, rows):
    if columns_needs_casting(data_types):
        cast = columns_cast(data_types)
        return [[c(v) for c, v in zip(cast, r)] for r in rows]
    return rows


#
def check_table_data_types(columns, data_types, rows):
    dts = [data_type_python[dt] for dt in data_types]
    for row_idx, r in enumerate(rows):
        wrong_cls = [c for c, dt, data in zip(columns, dts, r) if data is not None and not isinstance(data, dt)]
        if wrong_cls:
            print(f"Warning: row {row_idx + 1} has wrong data type in columns: {','.join(wrong_cls)}")
