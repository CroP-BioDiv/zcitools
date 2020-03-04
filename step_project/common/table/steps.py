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

import csv
import itertools
import os.path
from step_project.base_step import Step, StepCollection
from common_utils.exceptions import ZCItoolsValueError
from common_utils.show import print_table
from common_utils.file_utils import ensure_directory, append_line_to_file, read_file_as_list
from common_utils.value_data_types import KNOWN_DATA_TYPES


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


def _write_csv(filename, columns, rows):
    with open(filename, 'w', newline='') as outcsv:
        writer = csv.writer(outcsv, delimiter=';', quotechar='"')
        writer.writerow([n for n, _ in columns])  # Header
        writer.writerows(rows)


def _read_csv(filename):
    if os.path.isfile(filename):
        with open(filename, 'r') as incsv:
            reader = csv.reader(incsv, delimiter=';', quotechar='"')
            next(reader)  # Skip header
            return list(reader)
    else:
        return []


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
        if type_description:
            self._columns = type_description['columns']
        else:
            self._columns = None

    def _check_data(self):
        _check_columns(self._columns)

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
            _write_csv(self._get_table_filename(), self._columns, self._rows)

    # Retrieve data methods
    def has_column(self, column_name):
        return any(c == column_name for c, _ in self._columns)

    def get_column_with_data_types(self):
        return self._columns

    def get_column_names(self):
        return [c for c, _ in self._columns]

    def get_data_types(self):
        return [dt for _, dt in self._columns]

    def get_rows(self):
        if self._rows is None:
            self._rows = _read_csv(self._get_table_filename())
            _check_rows(self._columns, self._rows)
        return self._rows

    def _column_index(self, column_name):
        for i, (name, _) in enumerate(self._columns):
            if name == column_name:
                return i
        raise ZCItoolsValueError(f"No column named {column_name}! Columns: {self._columns}")

    def get_column_values(self, column_name):
        idx = self._column_index(column_name)
        return set(row[idx] for row in self.get_rows())

    def _column_index_by_type(self, data_type):
        for i, (_, dt) in enumerate(self._columns):
            if dt == data_type:
                return i
        raise ZCItoolsValueError(f"No column with data_type {data_type}! Columns: {self._columns}")

    def get_column_values_by_type(self, data_type):
        # Iterate through column values
        idx = self._column_index_by_type(data_type)
        return set(row[idx] for row in self.get_rows())

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
            _write_csv(os.path.join(data_dir, group), self._columns[1:], rows)
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
                        [[group] + r for r in _read_csv(os.path.join(data_dir, group))]
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
