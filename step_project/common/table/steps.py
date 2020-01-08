import csv
import itertools
import os.path
from step_project.base_step import Step, StepCollection
from common_utils.exceptions import ZCItoolsValueError
from common_utils.show import print_table
from .table_data import TableData

_COLUMN_TYPES = frozenset(['seq_ident', 'str', 'int', 'date', 'decimal'])


def _check_columns(columns, rows):
    if not isinstance(columns, (tuple, list)):
        raise ZCItoolsValueError(f"Column specification is not a list or tuple ({type(columns)})!")
    wrong_columns = [d for d in columns if not isinstance(d, (tuple, list)) or len(d) != 2]
    if wrong_columns:
        raise ZCItoolsValueError(f"Error in columns specification: {wrong_columns}")
    #
    wrong_types = [(n, ct) for n, ct in columns if ct not in _COLUMN_TYPES]
    if wrong_types:
        raise ZCItoolsValueError(f"Unsupported column type(s) for: {wrong_types}")
    #
    if rows:
        n_cols = len(columns)
        diff_lengts = [(i, len(row)) for i, row in enumerate(rows) if len(row) != n_cols]
        if diff_lengts:
            raise ZCItoolsValueError(f"Rows have different length than specified columns {n_cols}: {diff_lengts}")


class TableStep(Step):
    """
Stores table data. Table has list of columns, where column has name and data type.
Column names and types are stored in description.yml.
Table data is stored in table.csv with header, separator ;, quote character ".
ToDo: store original file?
"""
    _STEP_TYPE = 'table'
    _TABLE_FILENAME = 'table.csv'

    # Init object
    def _init_data(self, type_description):
        self._rows = None
        if type_description:
            self._columns = type_description['columns']
            self._orig_filename = type_description.get('orig_filename')
            self._data_format = type_description.get('data_format')
        else:
            self._columns = None
            self._orig_filename = None
            self._data_format = None

    def _check_data(self):
        _check_columns(self._columns, self._rows)

    # Set data
    def set_table_data(self, rows, columns, orig_filename=None, data_format=None):
        self._rows = rows
        self._columns = columns
        self._orig_filename = orig_filename
        self._data_format = data_format
        #
        self._check_data()

    # Save/load data
    def _get_table_filename(self):
        return self.step_file(self._TABLE_FILENAME)

    def save(self):
        # Store description.yml
        desc = dict(columns=[list(c) for c in self._columns])
        if self._orig_filename:
            desc['orig_filename'] = self._orig_filename
            if self._data_format:
                desc['data_format'] = self._data_format

        self.save_description(desc, completed=bool(self._columns))

        # Write csv
        if self._rows:
            with open(self._get_table_filename(), 'w', newline='') as outcsv:
                writer = csv.writer(outcsv, delimiter=';', quotechar='"')
                writer.writerow([n for n, _ in self._columns])  # Header
                writer.writerows(self._rows)

    # Retrieve data methods
    def get_rows(self):
        if self._rows is None:
            if os.path.isfile(self._get_table_filename()):
                with open(self._get_table_filename(), 'r') as incsv:
                    reader = csv.reader(incsv, delimiter=';', quotechar='"')
                    next(reader)  # Skip header
                    self._rows = [row for row in reader]
                    _check_columns(self._columns, self._rows)
            else:
                self._rows = []
        return self._rows

    def get_table_data(self):
        return TableData(self._columns, self.get_rows())

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

        print('\nData:')
        print_table([c for c, _ in self._columns], self.get_rows(), show_limit=7)


class GroupTablesStep(StepCollection):
    """
Stores table data in more substeps (subdirectories).
First column is used for grouping.
Column names and types are stored in description.yml.
"""
    _STEP_TYPE = 'group_tables'
    _SUBSTEP_CLASS = TableStep

    # Init object
    def _init_data(self, type_description):
        if type_description:
            self.set_columns(type_description['columns'])
        else:
            self._columns = None

    # Set data
    def set_columns(self, columns):
        _check_columns(columns, None)
        self._columns = columns

    def set_group_data(self, substep, data):
        substep.set_table_data(data, self._columns[1:])

    # Save/load data
    def save(self):
        # Store description.yml
        self.save_description(dict(columns=[list(c) for c in self._columns or []]), completed=bool(self._columns))

    # Retrieve data methods
    def get_rows(self):
        return list(itertools.chain.from_iterable(substep.get_rows() for substep in self.step_objects()))

    def get_table_data(self):
        return TableData(self._columns, self.get_rows())
