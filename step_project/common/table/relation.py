
from common_utils.cache import cache
from common_utils.value_data_types import cast_table, check_table_data_types
from common_utils.import_method import import_pandas
from common_utils.show import print_table
from common_utils.exceptions import ZCItoolsValueError


class Relation:  # To have different name than 'table' :-)
    def __init__(self, rows, columns=None, data_types=None, columns_data_types=None,
                 cast_data=True, check_data_types=True):
        if columns_data_types:
            self._columns, self._data_types = zip(*columns_data_types)
            self._columns = [c.lower() for c in self._columns]
        else:
            assert columns and data_types, (columns, data_types)
            self._columns = columns
            self._data_types = data_types
        #
        if cast_data:
            self._rows = cast_table(self._data_types, rows)
        else:
            self._rows = rows
        #
        if check_data_types:
            check_table_data_types(self._columns, self._data_types, self._rows)

        #
        self._column_idxs = dict((c, i) for i, c in enumerate(self._columns))

    def is_columns(self, c):
        return c in self._column_idxs

    def column_index(self, c):
        return self._column_idxs[c]

    def to_excel(self, filename):
        self.to_pandas_df().to_excel(filename, index=False)

    @cache
    def to_pandas_df(self):
        return import_pandas().DataFrame(self._rows, columns=self._columns)

    #
    def select(self, columns):
        columns = [c.lower() for c in columns]
        not_in = [c for c in columns if c not in self._column_idxs]
        if not_in:
            raise ZCItoolsValueError(f"Column(s) {', '.join(columns)} not in the relation!")
        idxs = [self._column_idxs[c] for c in columns]
        rows = [[r[i] for i in idxs] for r in self._rows]
        data_types = [self._data_types[i] for i in idxs]
        return Relation(rows, columns=columns, data_types=data_types, cast_data=False, check_data_types=False)

    def join(self, column, relation, rel_column):
        return JoinRelation(self, column, relation, rel_column)

    def where(self, filters):
        # Filters is python expression with columns as identifiers
        c = compile(filters, '<string>', 'eval')
        return Relation([r for r in self._rows if eval(c, None, dict(zip(self._columns, r)))],
                        columns=self._columns, data_types=self._data_types,
                        cast_data=False, check_data_types=False)

    def print_data(self):
        print_table(self._columns, self._rows, show_limit=7)


class JoinRelation(Relation):
    def __init__(self, relation_1, column_1, relation_2, column_2):
        if not relation_1.is_columns(column_1):
            raise ZCItoolsValueError(f"Column {column_1} not in JOIN left relation!")
        if not relation_2.is_columns(column_2):
            raise ZCItoolsValueError(f"Column {column_2} not in JOIN right relation!")

        #
        idx_1 = relation_1.column_index(column_1)
        idx_2 = relation_2.column_index(column_2)
        # ToDo: izbaciti kolonu iz druge tablice
        columns = relation_1._columns + relation_2._columns
        data_types = relation_1._data_types + relation_2._data_types
        rows = []
        super().__init__(rows, columns=columns, data_types=data_types, cast_data=False, check_data_types=False)
