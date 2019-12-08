import os.path
from ..utils.file_utils import read_yaml
from ..utils.exceptions import ZCItoolsValueError

# Step classes. Used for filling _type_2_step_cls dict
from .table import TableStep
from .images import ImagesStep
from .sequences import SequencesStep
from .annotations import AnnotationsStep
from .alignments import AlignmentStep, AlignmentsStep
from .phylogenetics import RAxMLStep, RAxMLSteps, MrBayesStep, MrBayesSteps

_type_2_step_cls = dict((cls._STEP_TYPE, cls) for cls in locals().values() if hasattr(cls, '_STEP_TYPE'))


def read_step(step_name, check_data_type=None, update_mode=False):
    if isinstance(step_name, str):
        desc_data = read_yaml(os.path.join(step_name, 'description.yml'))
    else:
        assert isinstance(step_name, list), type(step_name)
        desc_data = read_yaml(os.path.join(*step_name, 'description.yml'))
    if not desc_data:
        raise ZCItoolsValueError(f"'{step_name}' is not a step!")

    data_type = desc_data['data_type']

    if check_data_type:
        if isinstance(check_data_type, str):
            if check_data_type != data_type:
                raise ZCItoolsValueError(f"Step {step_name} is not of data type '{check_data_type}'!")
        else:
            if data_type not in check_data_type:
                raise ZCItoolsValueError(f"Step {step_name} is not of data types: {', '.join(check_data_type)}!")

    cls = _type_2_step_cls.get(data_type)
    if not cls:
        raise ZCItoolsValueError(f"No step class for data type {data_type}!")

    return cls(desc_data['project'], update_mode=update_mode)
