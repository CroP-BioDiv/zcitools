import sqlite3
from collections import defaultdict
from .value_data_types import data_type_sqlite
from .exceptions import ZCItoolsValueError
from .file_utils import silent_remove_file
# from .misc import time_it


def create_table_from_step(cursor, step, table_name):
    # Create table
    c_dts = step.get_column_with_data_types()
    cls = ', '.join(f'{c} {data_type_sqlite[dt]}' for c, dt in c_dts)
    cursor.execute(f'CREATE TABLE {table_name} ({cls})')
    # Insert data
    cursor.executemany(f'INSERT INTO {table_name} VALUES ({",".join(["?"] * len(c_dts))})', step.get_rows())


def create_db_from_step(db_filename, step, table_name='t1'):
    silent_remove_file(db_filename)
    conn = sqlite3.connect(db_filename)
    create_table_from_step(conn.cursor(), step, table_name=table_name)
    conn.commit()  # ?
    conn.close()


class StepDatabase:
    # @time_it
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

            create_table_from_step(self.cursor, step, table_name)

            # Change table name
            table_name = chr(ord(table_name) + 1)
        # self.conn.commit()  # ?

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

    def exact_column_name(self, c, c_idx):
        # Returns tuple (table_name.column_name, column name, data type)
        # ToDo: better parsing to support functions. Right now functions can work if there are no spaces in them!
        if '(' in c:
            return c, f'c_{c_idx}', 'int'
        #
        fields = c.split('.')
        if len(fields) == 1:
            tables = self.column_2_tables.get(c)
            if not tables:
                raise ZCItoolsValueError(f'Column {c} is not good specified!')
            if len(tables) > 1:
                raise ZCItoolsValueError(f'Column {c} is in more tables!', ', '.join(t for t, _ in tables))
            #
            table_name, dt = tables[0]
            return f'{table_name}.{c}', c, dt
        #
        elif len(fields) == 2:
            t_name, c_name = fields
            if t_name not in self.table_2_step:
                raise ZCItoolsValueError(f'No table {t_name} for column {c}!')
            dts = [dt for t, dt in self.column_2_tables.get(c_name, []) if t == t_name]
            if len(dts) != 1:
                raise ZCItoolsValueError(f'Column {c_name} is not in table {t_name}!')
            return c, c_name, dts[0]
        #
        raise ZCItoolsValueError(f'Column {c} is not good specified!')

    def select_result(self, sql):
        self.cursor.execute(sql)
        return self.cursor.fetchall()

    #
    # @time_it
    def select_all_tables(self, select, where_part='', group_by_part='', having_part='', order_by_part='', info=False):
        # Returns table data generated as result of SELECT statement generated with given data
        # Returns tuple (column_data_types, rows)
        # select is None or string of format {[<table>.]column [AS name],}+
        select_part = []
        column_data_types = []
        if select:
            # ToDo: better parsing to support functions. Right now functions can work if there are no spaces in them!
            for c_idx, s in enumerate(select.split(',')):
                fields = s.strip().split()
                sel_name, c_name, data_type = self.exact_column_name(fields[0], c_idx)
                select_part.append(sel_name)
                if len(fields) == 1:
                    column_data_types.append((c_name, data_type))
                elif len(fields) == 3:
                    assert fields[1].lower() == 'as'
                    column_data_types.append((fields[2].lower(), data_type))
                else:
                    raise ZCItoolsValueError(f"Wrong column name: {s}")

        else:
            for t, column_dts in sorted(self.table_columns.items()):
                select_part.extend(f'{t}.{c}' for c, _ in column_dts)
                column_data_types.extend(column_dts)

        # ToDo: check for same column names. Rename them with table prefix!
        select_part = ', '.join(select_part)
        from_part = ', '.join(self.all_tables())
        where_part = f'WHERE {where_part}' if where_part else ''
        group_by_part = f'GROUP BY {group_by_part}' if group_by_part else ''
        having_part = f'HAVING {having_part}' if having_part else ''
        order_by_part = f'ORDER BY {order_by_part}' if order_by_part else ''
        sql = f'SELECT {select_part} FROM {from_part} {where_part} {group_by_part} {having_part} {order_by_part}'
        if info:
            print(f'SELECT {select_part}\nFROM {from_part}\n{where_part} {group_by_part} {having_part} {order_by_part}')
        return column_data_types, self.select_result(sql)
