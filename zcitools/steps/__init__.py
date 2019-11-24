import os.path
from ..utils.file_utils import read_yaml
from .table import TableStep
from .sequences import SequencesStep
from ..utils.exceptions import ZCItoolsValueError

_type_2_step_cls = dict((cls._STEP_TYPE, cls) for cls in locals().values() if hasattr(cls, '_STEP_TYPE'))


def read_step(step_name, check_data_type=None):
    desc_data = read_yaml(os.path.join(step_name, 'description.yml'))
    if not desc_data:
        raise ZCItoolsValueError(f"'{step_name}' is not a step!")

    data_type = desc_data['data_type']

    if check_data_type and check_data_type != data_type:
        raise ZCItoolsValueError(f"Step {ste_name} is not of data type '{check_data_type}'!")

    cls = _type_2_step_cls.get(data_type)
    if not cls:
        raise ZCItoolsValueError(f"No step class for data type {data_type}!")

    return cls(desc_data['project'])
