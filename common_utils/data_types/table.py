from ..value_data_types import columns_needs_casting, columns_cast, cast_table_rows
from ..show import print_hierarchical_table
from ..import_method import import_pandas


def _extend_rows(nc, rows):
    assert all(len(r) <= nc for r in rows)
    return [(r if len(r) == nc else (r + [''] * (nc - len(r)))) for r in rows]


class Table:
    # Data type table.
    # All data is stored into lists to support update methods.
    def __init__(self, columns, data_types, rows=None):
        # columns is list/tuple of columns names
        # data_types is list/tuple of column data types
        self._num_columns = len(columns)
        assert self._num_columns == len(self._data_types), (self._num_columns, len(self._data_types))
        self._columns = list(columns)
        self._data_types = list(data_types)
        self._cast_rows = columns_cast(data_types) if columns_needs_casting(data_types) else None

        #
        self._rows = []
        if rows:
            self.extend(rows)

    # Insert data methods
    def append_row(self, row):
        assert len(row) == self._num_columns, (self._num_columns, len(row), row)
        if self._cast_rows:
            self._rows.append([c(v) for c, v in zip(self._cast_rows, row)])
        else:
            self._rows.append(list(row))

    def extend_rows(self, rows):
        assert all(self._num_columns == len(r) for r in rows), \
            (self._num_columns, [(len(r), r) for r in rows if len(r) != self._num_columns])
        cast = self._cast_rows
        if cast:
            self._rows.extend([c(v) for c, v in zip(cast, row)] for row in rows)
        else:
            self._rows.extend(list(row) for row in rows)

    # ToDo: add_column/remove_column


class HierarchicalTable:
    # Data type that represent hierarchical table. Something like:
    # c1_1 | c1_2 | c1_3 | ...
    #      | c2_1 | c2_2 | ...
    #      |      | c3_1 | ...
    def __init__(self, columns, data_types):
        # Columns is list of lists of column names. Data types is parallel to it.
        self._depth = len(columns)
        self._num_columns = [len(c) for c in columns]
        assert self._depth == len(data_types), (self._depth, len(data_types))
        assert all(n == len(dt) for n, dt in zip(self._num_columns, data_types)), \
            [(n, len(dt)) for n, dt in zip(self._num_columns, data_types)]
        #
        self._columns = [list(c) for c in columns]
        self._data_types = [list(dt) for dt in data_types]
        self._repr_columns = max(i + n for i, n in enumerate(self._num_columns))
        self._set_cast()

        #
        self._rows = []
        self._row_depths = []  # integers [0, depth)

    def _set_cast(self):
        self._cast_rows = [(columns_cast(dt) if columns_needs_casting(dt) else None) for dt in self._data_types]

    # Insert data methods
    def append_row(self, depth, row):
        assert len(row) == self._num_columns[depth], (depth, self._num_columns[depth], len(row), row)
        cast = self._cast_rows[depth]
        self._row_depths.append(depth)
        if cast:
            self._rows.append([c(v) for c, v in zip(cast, row)])
        else:
            self._rows.append(list(row))

    def extend_rows(self, depth, rows):
        n = self._num_columns[depth]
        assert all(n == len(r) for r in rows), (depth, n, [(len(r), r) for r in rows if len(r) != n])
        cast = self._cast_rows[depth]
        self._row_depths.extend([depth] * len(rows))
        if cast:
            self._rows.extend([c(v) for c, v in zip(cast, row)] for row in rows)
        else:
            self._rows.extend(list(row) for row in rows)

    # Change columns data and type
    def update_column(self, depth, column_name, method, update_column_name=None, update_data_type=None):
        c_idx = self._columns[depth].index(column_name)
        for d, row in zip(self._row_depths, self._rows):
            if d == depth:
                row[c_idx] = method(row[c_idx])
        if update_column_name:
            self._columns[depth][c_idx] = update_column_name
        if update_data_type:
            self._data_types[depth][c_idx] = update_data_type
            self._set_cast()

    # Represent data
    def _extend_data(self):
        columns = [([''] * i + c) for i, c in enumerate(self._columns)]
        rows = [([''] * d + r) for d, r in zip(self._row_depths, self._rows)]
        return _extend_rows(self._repr_columns, columns), _extend_rows(self._repr_columns, rows)

    def print(self, show_limit=None):
        columns, rows = self._extend_data()
        print_hierarchical_table(columns, rows, show_limit=show_limit)

    def to_excel(self, filename):
        columns, rows = self._extend_data()
        df = import_pandas().DataFrame(columns + rows, columns=None)
        df.to_excel(filename, index=False)
