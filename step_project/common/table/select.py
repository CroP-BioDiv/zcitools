import itertools
from common_utils.step_database import StepDatabase


def select_data(step_data, result, steps, select, where_part, group_by_part, having_part, order_by_part):
    # Fetch data
    with StepDatabase(steps) as db:
        select_columns = []  # tuples (select column name, column name, data type)
        if select:
            for s in select.split(','):
                fields = s.strip().split()
                sel_name, data_type = db.exact_column_name(fields[0])
                if len(fields) == 1:
                    select_columns.append((sel_name, sel_name.split('.')[1], data_type))
                elif len(fields) == 2:
                    select_columns.append((sel_name, fields[1].lower(), data_type))
                elif len(fields) == 3:
                    assert fields[1].lower() == 'as'
                    select_columns.append((sel_name, fields[2].lower(), data_type))

        else:
            select_columns = list(
                itertools.chain.from_iterable(((f'{t}.{c}', c, dt) for c, dt in column_dts)
                                              for t, column_dts in sorted(db.table_columns.items())))

        # select_part = ', '.join(db.exact_column_name(s.strip()) for s in select.split(','))
        select_part = ', '.join(x[0] for x in select_columns)
        from_part = ', '.join(db.all_tables())
        where_part = f'WHERE {where_part}' if where_part else ''
        group_by_part = f'GROUP BY {group_by_part}' if group_by_part else ''
        having_part = f'HAVING {having_part}' if having_part else ''
        order_by_part = f'ORDER BY {order_by_part}' if order_by_part else ''
        sql = f'SELECT {select_part} FROM {from_part} {where_part} {group_by_part} {having_part} {order_by_part}'
        # print(f'SELECT {select_part}\nFROM {from_part}\n{where_part} {group_by_part} {having_part} {order_by_part}')
        data = db.select_result(sql)

    #
