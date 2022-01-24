"""
Steps types to store table data.
Table structure is defined with list of columns, where column has name and data type.
Column names and types are stored in description.yml.
Table data can be stored in different formats, depending on step implementation.

Notes:
 - loading of table data is lazy and cached.
 - ?table step stores column types but doesn't force casting data into them

Step interface regarding table data structure:
 - get_rows()
 - get_column_values(column)
 - get_column_values_by_type(data_type)
"""

import itertools
import os.path
from decimal import Decimal
from datetime import date
from step_project.base_step import Step, StepCollection
from common_utils.exceptions import ZCItoolsValueError
from common_utils.show import print_table
from common_utils.file_utils import ensure_directory, append_line_to_file, read_file_as_list, write_csv, read_csv
from common_utils.value_data_types import KNOWN_DATA_TYPES

_type_format = dict(
    int=lambda x: int(x),
    decimal=lambda x: Decimal(x),
    date=lambda x: date.fromisoformat(x[:10]),
)
_cvs_delimiter = ','


def _check_columns(columns):
    if not isinstance(columns, (tuple, list)):
        raise ZCItoolsValueError(f"Column specification is not a list or tuple ({type(columns)})!")
    wrong_columns = [d for d in columns if not isinstance(d, (tuple, list)) or len(d) != 2]
    if wrong_columns:
        raise ZCItoolsValueError(f"Error in columns specification: {wrong_columns}")
    #
    wrong_types = [(n, ct) for n, ct in columns if ct not in KNOWN_DATA_TYPES]
    if wrong_types:
        raise ZCItoolsValueError(f"Unsupported column type(s) for: {wrong_types}")


def _check_rows(columns, rows):
    n_cols = len(columns)
    diff_lengts = [(i, len(row)) for i, row in enumerate(rows) if len(row) != n_cols]
    if diff_lengts:
        raise ZCItoolsValueError(f"Rows have different length than specified columns {n_cols}: {diff_lengts}")
    # ToDo: data types?


# -------------------------------------------------------------
class TableStep(Step):
    """
Stores table data. Table has list of columns, where column has name and data type.
Column names and types are stored in description.yml.
Table data is stored in table.csv with header, separator ;, quote character ".

"""
    _STEP_TYPE = 'table'
    _TABLE_FILENAME = 'table.csv'

    # Init object
    def _init_data(self, type_description):
        self._rows = None
        self._columns = type_description['columns'] if type_description else None

    def _check_data(self):
        _check_columns(self._columns)

    def rows_as_dicts(self):
        idxs = self.get_column_indices()
        return (_RowReader(idxs, r) for r in self.get_rows())

    # Set data
    def set_columns(self, columns):
        _check_columns(columns)
        self._columns = columns

    def set_table_data(self, rows, columns):
        self._rows = rows
        self._columns = columns
        #
        self._check_data()

    # Save/load data
    def _get_table_filename(self):
        return self.step_file(self._TABLE_FILENAME)

    def save(self, completed=True):
        # Store description.yml
        self.save_description(dict(columns=[list(c) for c in self._columns or []]),
                              completed=completed and bool(self._columns))

        # Write csv
        if self._rows:
            write_csv(self._get_table_filename(), self._columns, self._rows, delimiter=_cvs_delimiter)

    # Retrieve data methods
    def has_column(self, column_name):
        return any(c == column_name for c, _ in self._columns)

    def choose_first_column(self, *columns, error=False):
        for c in columns:
            if self.has_column(c):
                return c
        if error:
            raise ZCItoolsValueError(
                f'No column found from: {columns}.\nExisting columns: {[c for c, _ in self._columns]}')

    def get_column_with_data_types(self):
        return self._columns

    def get_column_names(self):
        return [c for c, _ in self._columns]

    def get_data_types(self):
        return [dt for _, dt in self._columns]

    def get_rows(self):
        if self._rows is None:
            filename = self._get_table_filename()
            if os.path.isfile(filename):
                self._rows = read_csv(filename)
                for i, (_, dt) in enumerate(self._columns):
                    if ff := _type_format.get(dt):
                        for r in self._rows:
                            if r[i]:
                                r[i] = ff(r[i])
                _check_rows(self._columns, self._rows)
            else:
                self._rows = []
        return self._rows

    def num_rows(self):
        return len(self.get_rows())

    def column_index(self, column_name):
        for i, (name, _) in enumerate(self._columns):
            if name == column_name:
                return i
        raise ZCItoolsValueError(f"No column named {column_name}! Columns: {self._columns}")

    def get_column_indices(self):
        return dict((c, i) for i, (c, _) in enumerate(self._columns))

    def get_column_values(self, column_name):
        idx = self.column_index(column_name)
        return set(row[idx] for row in self.get_rows())

    def _column_index_by_type(self, data_type, column_name=None):
        for i, (c, dt) in enumerate(self._columns):
            if dt == data_type and (column_name is None or c == column_name):
                return i
        raise ZCItoolsValueError(f"No column with data_type {data_type}! Columns: {self._columns}")

    def get_column_values_by_type(self, data_type, column_name=None):
        # Iterate through column values
        idx = self._column_index_by_type(data_type, column_name=column_name)
        return set(row[idx] for row in self.get_rows())

    def mapping_between_columns(self, from_column, to_column):
        idx_from = self.column_index(from_column)
        idx_to = self.column_index(to_column)
        return dict((r[idx_from], r[idx_to]) for r in self.get_rows())

    def mapping_column_2_columns(self, idx_column, *to_columns):
        idx_from = self.column_index(idx_column)
        idx_to = [self.column_index(c) for c in to_columns]
        return dict((r[idx_from], (r[x] for x in idx_to)) for r in self.get_rows())

    def index_on_table(self, *columns):
        return IndexOnTable(self.get_column_names(), self.get_rows(), *columns)

    def select(self, columns):  # ToDo: other parts of SELECT statement? where, order
        idxs = [self.column_index(c) for c in columns]
        return ([row[i] for i in idxs] for row in self.get_rows())

    # Show data
    def show_data(self, params=None):
        if not self.is_completed():
            print('Table step is not completed!')
            return

        print('Columns:')
        print_table(None, self._columns)

        if 'columns' not in params:
            print('\nData:')
            print_table([c for c, _ in self._columns], self.get_rows(), show_limit=7)

    def to_excel(self, filename, header=True):
        from common_utils.import_method import import_pandas
        df = import_pandas().DataFrame(self.get_rows(), columns=self.get_column_names())
        df.to_excel(filename, index=False, header=header)

    def to_sqlite(self, db_filename):
        from common_utils.step_database import create_db_from_step
        create_db_from_step(db_filename, self)


class TableGroupedStep(TableStep):
    """
First column is used for grouping.
Data for each group is stored as csv into file data/<group>.
Groups that do not have data are store as a list in file no_data.txt
"""
    _STEP_TYPE = 'table_grouped'
    _DATA_SUBDIRECTORY = 'data'
    _NO_DATA_FILE = 'no_data.txt'

    # Init object
    def _init_data(self, type_description):
        self._rows = None
        if type_description:
            self.set_columns(type_description['columns'])
        else:
            self._columns = None

    def _data_subdirectory(self):
        return self.step_file(self._DATA_SUBDIRECTORY)

    def _no_data_filename(self):
        return self.step_file(self._NO_DATA_FILE)

    # Set data
    def set_group_rows(self, group, rows):
        # ToDo: check rows?
        if rows:
            self._rows = None  # Remove get_rows() cache
            data_dir = self._data_subdirectory()
            ensure_directory(data_dir)
            write_csv(os.path.join(data_dir, group), self._columns[1:], rows)
        else:
            append_line_to_file(self._no_data_filename(), group)

    # Save/load data
    def save(self):
        # Store description.yml
        self.save_description(dict(columns=[list(c) for c in self._columns or []]), completed=bool(self._columns))

    # Retrieve data methods
    def get_rows(self):
        # ToDo: dodati prvu kolonu
        if self._rows is None:
            data_dir = self._data_subdirectory()
            if os.path.isdir(data_dir):
                self._rows = list(
                    itertools.chain.from_iterable(
                        [[group] + r for r in read_csv(os.path.join(data_dir, group))]
                        for group in os.listdir(data_dir)))
            else:
                self._rows = []
        return self._rows

    def get_groups_with_rows(self):
        data_dir = self._data_subdirectory()
        return os.listdir(data_dir) if os.path.isdir(data_dir) else []

    def get_groups_without_rows(self):
        no_data_file = self._no_data_filename()
        return read_file_as_list(no_data_file) if os.path.isfile(no_data_file) else []

    def known_groups(self):
        return self.get_groups_with_rows() + self.get_groups_without_rows()


#
class Rows2Table:
    # Transfer raw data (column names and rows) into formated data for table step
    def __init__(self, column_description, column_names=None):
        # column_names are names of input columns
        # column_description are descriptions of possible columns to transfer with description how to format data
        # column description is dict with attrs:
        # * column   : input column name
        # * output   : if omitted, same as input column name
        # * outputs  : split column into more columns. Transfer should return tuple!
        # * optional : default False
        # * type     : column type, from KNOWN_DATA_TYPES. Default 'str'
        # * tranfer  : format method
        # * value    : value of static column
        # * check    : callable that returns True if record has to be added.

        # Same as in TableStep
        self._columns = []
        self._rows = []
        self._formaters = []  # Tuples (column index, None or format callable, None or check value callable)

        #
        if column_names:
            not_in = []
            for d in column_description:
                col = d.get('column')

                # Column with static value
                if not col:
                    output_col = d.get('output')
                    if not output_col:
                        print('Warning: description without column and output columns!', d)
                        continue
                    _type = d.get('type', 'str')
                    assert _type in KNOWN_DATA_TYPES, d
                    self._columns.append((output_col, _type))

                    # Format method
                    self._formaters.append((None, d.get('value'), None))
                    continue

                # Column with value taken from input data
                try:
                    idx = column_names.index(col)
                except ValueError:
                    if not d.get('optional'):
                        not_in.append(col)
                    continue

                # Column name and type
                output_col = d.get('output', col)
                _type = d.get('type', 'str')
                assert _type in KNOWN_DATA_TYPES, d
                self._columns.append((output_col, _type))

                # Format method
                self._formaters.append(
                    (idx, d.get('transfer', _type_format.get(_type)), d.get('check')))

            #
            if not_in:
                raise ZCItoolsValueError(f'Mandatory columns not presented in csv file: {not_in}')

    def set_rows(self, rows):
        for row in rows:
            out_row = []
            for idx, tranfer, check in self._formaters:
                if idx is None:  # Static column
                    out_row.append(tranfer)
                else:
                    d = row[idx]
                    if check and not check(d):
                        out_row = None
                        break
                    out_row.append(tranfer(d) if tranfer else d)
            if out_row:
                self._rows.append(out_row)

    def in_table_step(self, table_step):
        table_step.set_table_data(self._rows, self._columns)


class IndexOnTable:
    def __init__(self, columns, rows, *index):
        self._column_idxs = dict((c, i) for i, c in enumerate(columns))
        if len(index) == 1:
            idx = self._column_idxs[index[0]]
            self._data = dict((r[idx], r) for r in rows)
        else:
            idxs = [self._column_idxs[c] for c in columns]
            self._data = dict((tuple(r[i] for i in idxs), r) for r in rows)

    def __len__(self):
        return len(self._data)

    def get_row(self, index_value):
        return self._data[index_value]

    def get_cell(self, index_value, column_name):
        return self._data[index_value][self._column_idxs[column_name]]

    def get_cells(self, index_value, *column_names):
        row = self._data[index_value]
        return tuple(row[self._column_idxs[c]] for c in column_names)

    def iterate_column(self, column_name):
        idx = self._column_idxs[column_name]
        for k, row in self._data.items():
            yield k, row[idx]


class _RowReader:
    def __init__(self, idxs, row):
        self._idxs = idxs
        self._row = row

    def __getitem__(self, c):
        return self._row[self._idxs[c]]
