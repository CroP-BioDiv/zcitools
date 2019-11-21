import os.path
from zcitools.utils.file_utils import ensure_directory, remove_directory, write_yaml, read_yaml
from ..utils.exceptions import ZCItoolsValueError

"""
Step object is created with step data that specfies step subdirectory and other environment info.
Subdirectory doesn't have to exist. In that case base class will create it.
"""


class Step:
    _STEP_TYPE = None
    _COLUMN_TYPES = frozenset(['ncbi_ident', 'str', 'int'])

    def __init__(self, step_data, remove_data=False):
        self._step_data = step_data
        self._step_name = step_data['step_name']

        # Call init data method
        if remove_data:
            remove_directory(self._step_name, create=True)
            d = None
        else:
            d = self.get_desription()
            if not d:
                ensure_directory(self._step_name)
        if d:
            if d['data_type'] != self._STEP_TYPE:
                raise ZCItoolsValueError(
                    f"Step class of tyep '{self._STEP_TYPE}' created with data of type '{d['data_type']}'!")
            type_desc = d['data']
        else:
            type_desc = None
        #
        self._init_data(type_desc)

    def _init_data(self, type_description):
        raise NotImplementedError(f'Method {self.__class__.__name__}._init_data() is not implemented!')

    #
    def step_file(self, f):
        return os.path.join(self._step_name, f)

    def step_files(self):
        return os.listdir(self._step_name)

    #
    def get_step_type(self):
        return self._STEP_TYPE

    # Description methods
    def _store_description(self, type_description):
        write_yaml(dict(data_type=self.get_step_type(), data=type_description, project=self._step_data),
                   self.step_file('description.yml'))

    def get_desription(self):
        return read_yaml(self.step_file('description.yml'))

    def get_type_desciption(self):
        d = self.get_desription()
        if d:
            return d['data']
