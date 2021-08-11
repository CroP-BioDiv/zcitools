import os.path
import csv
from .steps import TableStep
from common_utils.file_utils import filetype_from_ext
from common_utils.exceptions import ZCItoolsValueError
from common_utils.value_data_types import column_name_2_type


def create_table_step(project, step_data, params):
    if not os.path.isfile(params.filename):
        raise ZCItoolsValueError(f"Table file {params.filename} doesn't exist.")

    # Find how to read data
    data_format = params.data_format
    if data_format is None:
        data_format = filetype_from_ext(params.filename)
    if not data_format:
        raise ZCItoolsValueError(f"Data format for input table is not specified or found! Filename {params.filename}.")

    columns = [x.split(',') for x in params.columns.split(':')] if params.columns else None
    if columns:
        columns = [(c[0], column_name_2_type(c[0])) if len(c) == 1 else c for c in columns]

    # Read data.
    data = None
    if data_format == 'text':
        # ToDo: separator for more columns. For now only list supported
        with open(params.filename, 'r') as r:
            data = [[line] for line in filter(None, (_l.strip() for _l in r.readlines()))]
        data = sorted(data)
    elif data_format == 'csv':
        with open(params.filename, 'r') as incsv:
            reader = csv.reader(incsv, delimiter=params.delimiter, quotechar='"')
            if params.has_header:
                header = next(reader)  # Skip header
                if not columns:
                    columns = [(c, column_name_2_type(c)) for c in header]  # Default
            data = sorted(reader)
    elif data_format == 'raw_data':
        data = [[line] for line in filename.split(';') if line]
    else:
        raise ZCItoolsValueError(f'Data format {data_format} is not supported!')

    if not columns:
        raise ZCItoolsValueError(f"Columns are not specified for input table! Filename {params.filename}.")

    # Store (or overwrite) step data
    step = TableStep(project, step_data, remove_data=True)
    step.set_table_data(data, columns)
    step.save()
    return step
