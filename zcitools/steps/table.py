import csv
import os.path
from .step import Step
from ..utils.exceptions import ZCItoolsValueError


class TableStep(Step):
    _STEP_TYPE = 'table'
    _TABLE_FILENAME = 'table.csv'

    def __init__(self, step_data, data=None, columns=None, orig_filename=None, data_format=None):
        # Arguments:
        #  * data    : list/tuple of lists/tuples
        #  * columns : list of tuples (column name, data type)
        assert bool(data) == bool(columns), (step_data, data, columns)
        super().__init__(step_data)

        self._orig_filename = orig_filename
        self._data_format = data_format

        if data:
            self._data = data
            self._columns = columns
        else:
            self._data = self._load_table()
            desc = self.get_type_desciption()
            self._columns = desc['columns']
            self._orig_filename = desc.get('orig_filename')
            self._data_format = desc.get('data_format')
        #
        self._check_data()

    def _check_data(self):
        n_cols = len(self._columns)
        diff_lengts = [(i, len(row)) for i, row in enumerate(self._data) if len(row) != n_cols]
        if diff_lengts:
            raise ZCItoolsValueError(f"Rows have different length than specified columns {n_cols}: {diff_lengts}")
        self.check_columns(self._columns)

    @classmethod
    def check_columns(cls, columns):
        wrong_types = [(n, ct) for n, ct in columns if ct not in cls._COLUMN_TYPES]
        if wrong_types:
            raise ZCItoolsValueError(f"Unsupported column type(s) for: {wrong_types}")

    #
    def _get_table_filename(self):
        return os.path.join(self._step_name, self._TABLE_FILENAME)

    def _load_table(self):
        with open(self._get_table_filename(), 'r') as incsv:
            reader = csv.reader(incsv)
            next(csvreader)
            return [row for row in csvreader]

    def save(self):
        # Store description.yml
        desc = dict(columns=self._columns)
        if self._orig_filename:
            desc['orig_filename'] = self._orig_filename
            if self._data_format:
                desc['data_format'] = self._data_format

        self._store_description(desc)

        # Write csv
        with open(self._get_table_filename(), 'w', newline='') as outcsv:
            writer = csv.writer(outcsv)
            writer.writerow([n for n, _ in self._columns])  # Header
            writer.writerows(self._data)
