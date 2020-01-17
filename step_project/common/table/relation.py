
from collections import defaultdict
from common_utils.cache import cache, cache_args
from common_utils.value_data_types import cast_table, check_table_data_types
from common_utils.import_method import import_pandas
from common_utils.show import print_table
from common_utils.exceptions import ZCItoolsValueError


class Relation:  # To have different name than 'table' :-)
    def __init__(self, rows, columns=None, data_types=None, columns_data_types=None,
                 cast_data=True, check_data_types=True):
        if columns_data_types:
            cs, dts = zip(*columns_data_types)
            self._columns = tuple(c.lower() for c in cs)
            self._data_types = tuple(dts)
        else:
            assert columns and data_types, (columns, data_types)
            self._columns = tuple(columns)
            self._data_types = tuple(data_types)
        #
        if len(set(self._columns)) != len(self._columns):
            raise ZCItoolsValueError(f"More columns with same name in the relation!")
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

    def num_columns(self):
        return len(self._columns)

    #
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
        print('Rows:', len(self._rows))

    def print_data_grouped(self, column_levels, ident=''):
        # ToDo: this is not possible to make general and to work in all cases :-/
        # column_levels is list (join column name, column names of joined table)
        for row in self._rows:
            pass

    #
    @cache_args
    def index_by_column_idx(self, idx):
        d = defaultdict(set)
        for i, row in enumerate(self._rows):
            d[row[idx]].add(i)
        return d

    def index_by_column(self, c):
        return self.index_by_column_idx(self._column_idxs[c])


class JoinRelation(Relation):
    def __init__(self, relation_1, column_1, relation_2, column_2):
        if not relation_1.is_columns(column_1):
            raise ZCItoolsValueError(f"Column {column_1} not in JOIN left relation!")
        if not relation_2.is_columns(column_2):
            raise ZCItoolsValueError(f"Column {column_2} not in JOIN right relation!")

        #
        idx_1 = relation_1.column_index(column_1)
        idx_2 = relation_2.column_index(column_2)
        # ToDo: remove ambiguous columns from second relation
        rem_cls = [i for i, c in enumerate(relation_2._columns) if relation_1.is_columns(c)]
        rem_cls.append(idx_2)

        def _rm_cl2(row):
            return [d for i, d in enumerate(row) if i not in rem_cls]

        c2_index = relation_2.index_by_column_idx(idx_2)
        rows_2 = relation_2._rows
        rows = []
        for row_1 in relation_1._rows:
            i2 = c2_index.get(row_1[idx_1])
            if i2:
                rows.extend(row_1 + _rm_cl2(rows_2[i]) for i in i2)
        # Store for later
        # self._relation_2 = relation_2
        # self._c2_index = c2_index
        self._added_columns = tuple(_rm_cl2(relation_2._columns))
        #
        super().__init__(
            rows,
            columns=relation_1._columns + self._added_columns,
            data_types=relation_1._data_types + tuple(_rm_cl2(relation_2._data_types)),
            cast_data=False, check_data_types=False)
