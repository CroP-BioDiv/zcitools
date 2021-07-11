from datetime import date
from decimal import Decimal
from common_utils.import_method import import_pandas


KNOWN_DATA_TYPES = frozenset(['seq_ident', 'str', 'int', 'date', 'decimal'])


def _ident_cast(x):
    return x


def fromisoformat(d):
    return datetime.date.fromisoformat(data['Date'])


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

data_type_sqlite = dict(
    seq_ident='text',
    str='text',
    int='integer',
    date='text',
    decimal='numeric',  # 'real',
)


#
def check_table_data_types(columns, data_types, rows):
    dts = [data_type_python[dt] for dt in data_types]
    for row_idx, r in enumerate(rows):
        wrong_cls = [c for c, dt, data in zip(columns, dts, r) if data is not None and not isinstance(data, dt)]
        if wrong_cls:
            print(f"Warning: row {row_idx + 1} has wrong data type in columns: {','.join(wrong_cls)}")


def columns_needs_casting(data_types):
    return any(needs_casting[dt] for dt in data_types)


def columns_cast(data_types):
    return [data_type_cast[dt] for dt in data_types]


def cast_table_rows(data_types, rows):
    if columns_needs_casting(data_types):
        cast = columns_cast(data_types)
        return [[c(v) for c, v in zip(cast, r)] for r in rows]
    return rows


def cast_table_data(data_types, rows):
    if columns_needs_casting(data_types):
        return cast_table_rows(data_types, rows)
    return rows


#
def column_name_2_type(name):
    name = name.lower()
    if name == 'seq_ident':
        return 'seq_ident'
    if name in ('date', 'datum'):
        return 'date'
    # ToDo: ...
    return 'str'


#
def table_data_2_pandas(column_data_types, rows):
    columns, data_types = zip(*column_data_types)
    rows = cast_table_data(data_types, rows)
    return import_pandas().DataFrame(rows, columns=columns)


def table_data_2_excel(filename, column_data_types, rows):
    table_data_2_pandas(column_data_types, rows).to_excel(filename, index=False)


def rows_2_excel(filename, columns, rows, index=False):
    import_pandas().DataFrame(rows, columns=columns).to_excel(filename, index=index)


def sheets_2_excel(filename, sheets, index=False):
    # Sheets is iterable of tuples (sheet_name, columns, rows)
    pandas = import_pandas()
    writer = pandas.ExcelWriter(filename)
    for sheet_name, columns, rows in sheets:
        df = pandas.DataFrame(rows, columns=columns)
        df.to_excel(writer, sheet_name=sheet_name, index=index)
    writer.save()
