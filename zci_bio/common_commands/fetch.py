import os.path
from common_utils.exceptions import ZCItoolsValueError


def fetch_common_db_data(step_data, table_step, step_type, common_db):
    step = table_step.project.new_step_by_type(step_type, step_data, remove_data=True)

    for seq_ident in table_step.get_column_values_by_type('seq_ident'):
        f = common_db.get_record(seq_ident, step.directory, info=True)
        if not f:
            raise ZCItoolsValueError(f"There is not CommonDB record for seq ident {seq_ident}!")
        step.add_sequence_file(os.path.basename(f))

    step.save()
    return step
