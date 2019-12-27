import os.path
from zcitools.steps.table import TableStep
from zcitools.utils.file_utils import filetype_from_ext
from zcitools.utils.exceptions import ZCItoolsValueError


def create_table_step(zcit, step_data, filename, data_format=None, columns=None):
    if not os.path.isfile(filename):
        raise ZCItoolsValueError(f"Table file {filename} doesn't exist.")

    # Find how to read data
    if data_format is None:
        data_format = filetype_from_ext(filename)
    if not data_format:
        raise ZCItoolsValueError(f"Data format for input table is not specified or found! Filename {filename}.")

    if columns:
        columns = [x.split(',') for x in columns.split(':')]

    # Read data.
    data = None
    if data_format == 'text':
        # ToDo: separator for more columns. For now only list supported
        with open(filename, 'r') as r:
            data = [[line] for line in filter(None, (_l.strip() for _l in r.readlines()))]
        data = sorted(data)
    else:
        raise ZCItoolsValueError(f'Data format {data_format} is not supported!')

    if not columns:
        raise ZCItoolsValueError(f"Columns are not specified for input table! Filename {filename}.")

    # Store (or overwrite) step data
    step = TableStep(zcit, step_data, remove_data=True)
    step.set_table_data(data, columns, orig_filename=filename, data_format=data_format)
    step.save()
    return step
