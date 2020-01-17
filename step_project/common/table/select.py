import itertools
from common_utils.step_database import StepDatabase
from common_utils.show import print_table
from common_utils.value_data_types import table_data_2_excel
from .steps import TableStep


def select_data(step_data, result, steps, select, where_part, group_by_part, having_part, order_by_part,
                output_filename=None):
    # Fetch data
    with StepDatabase(steps) as db:
        column_data_types, rows = db.select_all_tables(select, where_part, group_by_part, having_part, order_by_part)

    result = result.lower()
    if result[0] == 'e':    # Excel
        table_data_2_excel(output_filename or 'out.xslx', column_data_types, rows)

    elif result[0] == 's':  # Step
        step = TableStep(steps[0].project, step_data, remove_data=True)
        step.set_table_data(rows, column_data_types)
        step.save()
        return step

    else:                   # Print
        print('Columns:')
        print_table(None, column_data_types)

        print('\nData:')
        print_table([c for c, _ in column_data_types], rows, show_limit=7)
        print(f'\nNumber of rows: {len(rows)}')
