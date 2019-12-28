import csv
from step_project.base.step import Step
from common_utils.exceptions import ZCItoolsValueError
from common_utils.show import print_table


class TableStep(Step):
    """
Stores table data. Table has list of columns, where column has name and data type.
Column names and types are stored in description.yml.
Table data is stored in table.csv with header, separator ;, quote character ".
ToDo: store original file?
"""
    _STEP_TYPE = 'table'
    _COLUMN_TYPES = frozenset(['seq_ident', 'str', 'int', 'date', 'decimal'])
    _TABLE_FILENAME = 'table.csv'

    # Init object
    def _init_data(self, type_description):
        if type_description:
            self._columns = type_description['columns']
            self._orig_filename = type_description.get('orig_filename')
            self._data_format = type_description.get('data_format')
            if self._columns:
                self._data = self._load_table()
            else:
                self._data = None
        else:
            self._data = None
            self._columns = None
            self._orig_filename = None
            self._data_format = None

    def _check_data(self):
        if not isinstance(self._columns, (tuple, list)):
            raise ZCItoolsValueError(f"Column specification is not a list or tuple ({type(self._columns)})!")
        wrong_columns = [d for d in self._columns if not isinstance(d, (tuple, list)) or len(d) != 2]
        if wrong_columns:
            raise ZCItoolsValueError(f"Error in columns specification: {wrong_columns}")
        #
        n_cols = len(self._columns)
        diff_lengts = [(i, len(row)) for i, row in enumerate(self._data) if len(row) != n_cols]
        if diff_lengts:
            raise ZCItoolsValueError(f"Rows have different length than specified columns {n_cols}: {diff_lengts}")
        #
        wrong_types = [(n, ct) for n, ct in self._columns if ct not in self._COLUMN_TYPES]
        if wrong_types:
            raise ZCItoolsValueError(f"Unsupported column type(s) for: {wrong_types}")

    # Set data
    def set_table_data(self, data, columns, orig_filename=None, data_format=None):
        self._data = data
        self._columns = columns
        self._orig_filename = orig_filename
        self._data_format = data_format
        #
        self._check_data()

    # Save/load data
    def _get_table_filename(self):
        return self.step_file(self._TABLE_FILENAME)

    def _load_table(self):
        with open(self._get_table_filename(), 'r') as incsv:
            reader = csv.reader(incsv, delimiter=';', quotechar='"')
            next(reader)
            return [row for row in reader]

    def save(self):
        # Store description.yml
        desc = dict(columns=self._columns)
        if self._orig_filename:
            desc['orig_filename'] = self._orig_filename
            if self._data_format:
                desc['data_format'] = self._data_format

        self.save_description(desc, completed=bool(self._columns))

        # Write csv
        if self._columns:
            with open(self._get_table_filename(), 'w', newline='') as outcsv:
                writer = csv.writer(outcsv, delimiter=';', quotechar='"')
                writer.writerow([n for n, _ in self._columns])  # Header
                writer.writerows(self._data)

    # Retrieve data methods
    def _column_index_by_type(self, data_type):
        for i, (_, dt) in enumerate(self._columns):
            if dt == data_type:
                return i

    def get_column_by_type(self, data_type):
        # Iterate through column values
        idx = self._column_index_by_type(data_type)
        if idx is None:
            raise ZCItoolsValueError(f"No column with data_type {data_type}! Columns: {self._columns}")
        return (row[idx] for row in self._data)

    # Show data
    def show_data(self, params=None):
        if not self.is_completed():
            print('Table step is not completed!')
            return

        print('Columns:')
        print_table(None, self._columns)

        print('\nData:')
        print_table([c for c, _ in self._columns], self._data, show_limit=7)
