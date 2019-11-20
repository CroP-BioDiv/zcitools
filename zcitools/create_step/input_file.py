import os.path
from ..steps.table import TableStep
from ..utils.file_utils import filetype_from_ext
from ..utils.exceptions import ZCItoolsValueError


def create_table_step(step_data, filename, data_format=None, columns=None):
    if not os.path.isfile(filename):
        raise ZCItoolsValueError(f"Table file {filename} doesn't exist.")

    # Find how to read data
    if data_format is None:
        data_format = filetype_from_ext(filename)
    if not data_format:
        raise ZCItoolsValueError(f"Data format for input table is not specified or found! Filename {filename}.")

    if columns:
        columns = [x.split(',') for x in columns.split(':')]
        TableStep.check_columns(columns)

    # Read data.
    data = None
    if data_format == 'text':
        # ToDo: separator for more columns. For now only list supported
        with open(filename, 'r') as r:
            data = [[line] for line in filter(None, (_l.strip() for _l in r.readlines()))]
    else:
        raise ZCItoolsValueError(f'Data format {data_format} is not supported!')

    if not columns:
        raise ZCItoolsValueError(f"Columns are not specified for input table! Filename {filename}.")

    # Store step data
    step = TableStep(step_data, data=data, columns=columns, orig_filename=filename, data_format=data_format)
    step.save()
    return step
