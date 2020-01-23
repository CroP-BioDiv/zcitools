from common_utils.import_method import import_pandas
from common_utils.exceptions import ZCItoolsValueError


class CorrelationMatrix:
    def __init__(self, columns, list_values=None):
        self._columns = columns
        self._columns_lower = [c.lower() for c in columns]
        self._values = dict()  # tuple(sorted(c1, c2)) -> value. Note: value can be missing
        if list_values:
            assert len(columns) - 1 == len(list_values), (len(columns), len(list_values))
            for i, (c1, vs) in enumerate(zip(self._columns, list_values)):
                c2s = self._columns[i+1:]
                assert len(c2s) == len(vs), (i, c1, len(c2s), len(vs))
                for c2, v in zip(c2s, vs):
                    if v is not None:
                        self.set(c1, c2, v)

    def num_columns(self):
        return len(self._columns)

    def check_column(self, c, to_assert=False):
        c = c.lower()
        if c in self._columns_lower:
            return c
        if to_assert:
            assert False, (c, self._columns)

    #
    def set(self, c1, c2, v):
        c1 = self.check_column(c1, to_assert=True)
        c2 = self.check_column(c2, to_assert=True)
        k = (c1, c2) if c1 < c2 else (c2, c1)
        if v is None:
            self._values.pop(k)
        else:
            self._values[k] = v

    def get(self, c1, c2):
        c1 = self.check_column(c1, to_assert=True)
        c2 = self.check_column(c2, to_assert=True)
        k = (c1, c2) if c1 < c2 else (c2, c1)
        return self._values.get(k)

    @staticmethod
    def from_excel(filename, triangular='L'):
        df = import_pandas().read_excel(filename, sheetname='Sheet1')
        columns = list(df.columns[1:])
        if triangular.upper() == 'L':
            list_values = [list(df[c1][i+1:]) for i, c1 in enumerate(columns[:-1])]
        else:
            raise NotImplementedError('Upper triangular')
        return CorrelationMatrix(columns, list_values=list_values)

    @staticmethod
    def from_file(filename, triangular='L'):  # Lower/Upper triangular
        if filename.endswith('.xlsx'):
            return CorrelationMatrix.from_excel(filename, triangular=triangular)
        raise ZCItoolsValueError(f"Can't import correlation data from file {filename}!")
