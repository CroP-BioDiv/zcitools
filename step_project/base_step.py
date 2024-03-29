import os.path
import datetime
import re
from common_utils.file_utils import ensure_directory, remove_directory, silent_remove, write_yaml, read_yaml
from common_utils.exceptions import ZCItoolsValueError
from common_utils.show import print_ls_like_list

"""
Step object is created with step data that specfies step subdirectory and other environment info.
Subdirectory doesn't have to exist. In that case base class will create it.

Step subdirectory contains file description.yml which describes data step contains and environment project data.
description.yml contains data:
 - data_type : Step class type identifier. Value from <step class>._STEP_TYPE's
 - data      : Contains specific data description. Depends on step data type.
 - project   : Environment data, propagated from zcit.py script.
   - step_name  : Subdirectory name
   - prev_steps : List of step names used to produce (calculate) step's data
   - command    : zcit.py script command used to produce the step
   - cmd        : Arguments part of shell command (after zcit.py) used to produce the step.
"""


class Step:
    _STEP_TYPE = None
    _CACHE_PREFIX = '_c_'  # Cache files are prefixed with '_c_'
    _IS_COLLECTION = False

    def __init__(self, project, step_data, remove_data=False, update_mode=False, no_check=False, step_directory=None):
        assert project.__class__.__name__ == 'RunCommand', self.__class__.__name__  # For now
        self.project = project
        self._step_data = step_data
        # step_data['step_name'] is string or list of strings for substeps
        self._step_name_list = step_directory or \
            ([step_data['step_name']] if isinstance(step_data['step_name'], str) else step_data['step_name'])
        self.directory = os.path.join(*self._step_name_list)
        self._update_mode = update_mode

        # Call init data method
        if remove_data:
            remove_directory(self.directory, create=True)
            d = None
        else:
            d = self.get_description()
            if not d:
                ensure_directory(self.directory)
        if d:
            if d['data_type'] != self._STEP_TYPE:
                raise ZCItoolsValueError(
                    f"Step class of tyep '{self._STEP_TYPE}' created with data of type '{d['data_type']}'!")
            type_desc = d['data']
            # Update project data
            self._step_data.update((k, v) for k, v in d['project'].items() if k not in self._step_data)
        else:
            type_desc = None
        #
        self._init_data(type_desc)

        # Check data if exists and step not set in update mode
        if type_desc and not self._update_mode and self.is_completed() and not no_check:
            self._check_data()

    step_data_type = property(lambda self: self._STEP_TYPE)

    def _init_data(self, type_description):
        raise NotImplementedError(f'Method {self.__class__.__name__}._init_data() is not implemented!')

    def _check_data(self):
        print(f'Warning: method {self.__class__.__name__}._check_data() not implemented!')

    #
    def get_command(self):
        return self._step_data['command']

    def get_command_args(self):
        return self._step_data['command_args']

    def get_prev_steps(self):
        return self._step_data['prev_steps']

    def get_data_type(self):
        return self._STEP_TYPE

    #
    def get_step_name_prefix(self):
        return self._step_data.get('step_name_prefix')

    def set_step_name_prefix(self, n):
        self._step_data['step_name_prefix'] = n

    def get_base_step_name(self):
        if n := self._step_data.get('step_name'):
            if p := self.get_step_name_prefix():
                if p in n:
                    i = n.index(p)
                    return n[(i + len(p)):].strip('_')
            fields = n.split('_')
            if fields[0].isdigit():
                i = n.index(fields[0])
                return n[(i + len(fields[0])):].strip('_')
            return n

    def propagate_step_name_prefix(self, to_step):
        name_prefix = self.get_step_name_prefix()
        if name_prefix:
            to_step.set_step_name_prefix(name_prefix)

    #
    def common_db_identifier(self):
        c = self._step_data.get('common_db_identifier')
        if c:
            return tuple(c)

    def seq_ident_of_our_change(self, seq_ident, change):
        # Our change method is specified with lower case letter!
        assert isinstance(change, str) and len(change) == 1 and change.islower(), change
        fields = seq_ident.split('_')
        if fields[0].islower():  # Our changes already prepended, just add one more at the start
            return f'{change}{seq_ident}'
        if fields[0] == 'NC':  # Remove NC to have shorter ident. Longer than 10 chars can cause problems :-/
            return '_'.join([change] + fields[1:])
        return f'{change}_{seq_ident}'

    def is_completed(self):
        return self._step_data['completed']

    def can_be_completed(self):
        return self.is_file('output.zip')

    def get_local_name(self):
        return self._step_name_list[-1]

    # Description methods
    def save_description(self, type_description, create=True, completed=True):
        self._step_data['completed'] = completed
        pd = dict(self._step_data)
        if create:
            pd['created'] = datetime.datetime.now().isoformat()
            pd['updated'] = None
        else:
            pd['updated'] = datetime.datetime.now().isoformat()
        write_yaml(dict(data_type=self._STEP_TYPE, data=type_description, project=pd),
                   self.step_file('description.yml'))

    def get_description(self):
        return read_yaml(self.step_file('description.yml'))

    def get_type_description(self):
        d = self.get_description()
        if d:
            return d['data']

    def get_type_description_elem(self, attr, default=None):
        d = self.get_description()
        if d:
            return d['data'].get(attr, default)
        return default

    # Summary data
    def save_summary_data(self, data):
        assert isinstance(data, dict), data
        write_yaml(data, self.step_file('summary.yml'))

    def make_summary_data(self):
        return

    def get_summary_data(self):
        if os.path.isfile(f := self.step_file('summary.yml')):
            return read_yaml(f)
        # 'Cached' version
        if self.is_completed() and (d := self.make_summary_data()):
            self.save_summary_data(d)
            return d

    #
    def get_finish_data(self):
        if os.path.isfile(f := self.step_file('finish.yml')):
            return read_yaml(f)

    # Substep methods
    def get_substep_step_data(self, step_name):
        return dict(step_name=step_name)  # , prev_steps=None, command=None, command_args=None, cmd=None)

    def create_substep(self, step_cls, local_step_name, remove_data=False, update_mode=False):
        return step_cls(self.project,
                        self.get_substep_step_data(self._step_name_list + [local_step_name]),
                        remove_data=remove_data, update_mode=update_mode)

    def read_substep(self, local_step_name):
        return self.project.read_step(self._step_name_list + [local_step_name])

    # Commonn file methods
    def absolute_path(self):
        return os.path.abspath(self.directory)

    def step_file(self, *f):
        return os.path.join(self.directory, *f)

    def is_file(self, *f):
        return os.path.isfile(self.step_file(*f))

    def strip_step_dir(self, f):
        assert f.startswith(self.directory), (f, self.directory)
        return f[(len(self.directory) + 1):]

    def strip_step_dir_files(self, files):
        return [self.strip_step_dir(f) for f in files]

    @classmethod
    def _is_cache_file(cls, f):
        return f.startswith(cls._CACHE_PREFIX)

    def cache_file(self, f):
        return self.step_file(self._CACHE_PREFIX + f)

    def step_files(self, not_cached=False, matches=None):
        # Returns list of step's filenames relative to step subdirectory
        if not_cached:
            files = [f for f in os.listdir(self.directory) if not self._is_cache_file(f)]
        else:
            files = os.listdir(self.directory)
        if matches:
            pattern = re.compile(matches)
            files = [f for f in files if pattern.search(f)]
        return sorted(files)

    def step_dir_files(self, *dirs):
        d = self.step_file(*dirs)
        if os.path.isdir(d):
            return os.listdir(d)
        return []

    def step_subdirectories(self):
        return sorted(f for f in os.listdir(self.directory) if os.path.isdir(os.path.join(self.directory, f)))

    def run_module_name(self):
        pattern = re.compile(r'run_.*\.py')
        py_files = [f for f in os.listdir(self.directory) if pattern.search(f)]
        if len(py_files) == 1:
            return py_files[0]
        elif py_files:
            print(f'Warning: step {self.directory} contains more python run modules!', py_files)

    def remove_cache_files(self):
        for f in self.step_files():
            if self._is_cache_file(f):
                silent_remove(self.step_file(f))

    def get_common_db_records(self, common_db, record_idents, info=False):
        # Returns list of records that where not processed
        if common_db:
            return common_db.get_records(record_idents, self.directory, info=info)
        return record_idents

    def clean_files(self):
        print(f"Warning: step {self.directory} ({self.__class__.__name__}) doesn't have clean method!")

    # Common step methods/commands
    def show_data(self, params=None):
        print(f'{self.__class__.__name__}.show_data() not implemented!')

    def export_data(self, filename, filters=None, params=None):
        print(f'{self.__class__.__name__}.export_data() not implemented!')


#
class StepCollection(Step):
    """
Stores collection of steps of same step type.
Substep is identified by it's directory name.
Note: list of substeps is not stored in description.yml.
"""
    _IS_COLLECTION = True
    _SUBSTEP_CLASS = None  # Class of substep

    # Init object
    def _init_data(self, type_description):
        pass

    def _check_data(self):
        # Check all subobjects
        for obj in self.step_objects():
            obj._check_data()

    def step_objects(self):
        # Used to retrieve all step object encapsulated in this object.
        return (self.read_substep(n) for n in self.substep_names())

    # Save/load data
    def save(self, create=True, completed=True):
        # Store description.yml
        self.save_description(dict(), create=create, completed=completed)

    # Retrieve data methods
    def substep_names(self):
        return self.step_subdirectories()

    # Substep methods
    def create_substep(self, local_step_name, remove_data=False, update_mode=False):
        return self._SUBSTEP_CLASS(self.project, self.get_substep_step_data(self._step_name_list + [local_step_name]),
                                   remove_data=remove_data, update_mode=update_mode)

    def clean_files(self):
        for obj in self.step_objects():
            obj.clean_files()

    # Show data
    def show_data(self, params=None):
        print_ls_like_list(self._STEP_TYPE, self.substep_names(), sort=True, min_rows_to_split=20)
        if not self.is_completed():
            print('Step is not completed!')


class SimpleStep(Step):
    # Stores data without additional methods

    # Init object
    def _init_data(self, type_description):
        pass

    def _check_data(self):
        pass

    # Save/load data
    def save(self, data, create=True, completed=True):
        if data is None:
            data = self.get_type_description()
        self.save_description(data, create=create, completed=completed)

    # Show data
    def show_data(self, params=None):
        print(self.__class__.__name__, self.directory)
