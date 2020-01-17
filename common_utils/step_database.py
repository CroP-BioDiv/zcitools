import sqlite3
from collections import defaultdict
from .value_data_types import data_type_sqlite
from .exceptions import ZCItoolsValueError


class StepDatabase:
    def __init__(self, steps):
        self.conn = sqlite3.connect(':memory:')
        self.cursor = self.conn.cursor()
        self._open = True

        #
        self.table_2_step = dict()                # table_name -> step_object
        self.column_2_tables = defaultdict(list)  # column_name -> table_names
        self.table_columns = dict()               # table_name -> list of tuples (column name, data type)

        table_name = 'a'
        for step in steps:
            self.table_2_step[table_name] = step
            self.table_columns[table_name] = c_dts = step.get_column_with_data_types()

            # Table data
            for c, dt in c_dts:
                self.column_2_tables[c].append((table_name, dt))

            # Create table
            cls = ', '.join(f'{c} {data_type_sqlite[dt]}' for c, dt in c_dts)
            self.cursor.execute(f'CREATE TABLE {table_name} ({cls})')
            # Insert data
            self.cursor.executemany(
                f'INSERT INTO {table_name} VALUES ({",".join(["?"] * len(c_dts))})', step.get_rows())

            # Change table name
            table_name = chr(ord(table_name) + 1)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def close(self):
        if self._open:
            # self.conn.commit()  # ?
            self.conn.close()
            self._open = False

    def all_tables(self):
        return sorted(self.table_2_step.keys())

    def exact_column_name(self, c):
        # Returns tuple (table_name.column_name, data type)
        fields = c.split('.')
        if len(fields) == 1:
            tables = self.column_2_tables.get(c)
            if not tables:
                raise ZCItoolsValueError(f'Column {c} is not good specified!')
            if len(tables) > 1:
                raise ZCItoolsValueError(f'Column {c} is in more tables!', ', '.join(tables))
            #
            table_name, dt = tables[0]
            return f'{table_name}.{c}', dt
        #
        elif len(fields) == 2:
            t_name, c_name = fields
            if t_name not in self.table_2_step:
                raise ZCItoolsValueError(f'No table {t_name} for column {c}!')
            dts = [dt for t, dt in self.column_2_tables.get(c_name, []) if t == t_name]
            if len(dts) != 1:
                raise ZCItoolsValueError(f'Column {c_name} is not in table {t_name}!')
            return c, dts[0]
        #
        raise ZCItoolsValueError(f'Column {c} is not good specified!')

    def select_result(self, sql):
        self.cursor.execute(sql)
        return self.cursor.fetchall()
