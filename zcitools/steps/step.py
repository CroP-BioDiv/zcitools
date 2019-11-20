import os.path
from zcitools.utils.file_utils import write_yaml, read_yaml


class Step:
    _STEP_TYPE = None
    _COLUMN_TYPES = frozenset(['ncbi_ident', 'str', 'int'])

    def __init__(self, step_data):
        self._step_data = step_data
        self._step_name = step_data['step_name']

    def get_step_type(self):
        return self._STEP_TYPE

    # Description methods
    def _store_description(self, type_description):
        description = dict(project=self._step_data)
        description['data'] = type_description
        write_yaml(description, os.path.join(self._step_name, 'description.yml'))

    def get_desription(self):
        return read_yaml(os.path.join(self._step_name, 'description.yml'))

    def get_type_desciption(self):
        return self.get_desription()['data']
