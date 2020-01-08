class TableData:
    def __init__(self, columns, rows):
        self._columns = columns
        self._rows = rows
        self._column_idxs = dict((c, i) for i, (c, _) in enumerate(columns))
