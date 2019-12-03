import os.path
from zipfile import ZipFile, ZIP_BZIP2
from .file_utils import ensure_directory
from .exceptions import ZCItoolsValueError

"""
Cache object represent cache storage for one step data type, possible with specific parameters or environment.
Cache object point to specific directory, and stores records.
Record is named by step data identifier, which depends on step type.
Record is implemented as zip file. It can contain any kind of data: any number of files and/or directories.
"""


class Cache:
    def __init__(self, cache_dir):
        self.cache_dir = cache_dir
        self._dir_exists = os.path.exists(cache_dir)
        if self._dir_exists and os.path.isfile(cache_dir):
            raise ZCItoolsValueError(f'Cache location ({cache_dir}) is a file!')

    def _save_file(self, zip_f, f):
        # Note: step filenames are relative to project main directory. Zip filenames are stored without step_name.
        zip_name = os.path.join(*(f.split(os.path.sep)[1:]))
        zip_f.write(f, arcname=zip_name)

    def _save_directory(self, zip_f, d):
        # ToDo: is it needed to save directory into the zip?
        for f in os.listdir(d):
            f = os.path.join(d, f)
            if os.path.isfile(f):
                self._save_file(zip_f, f)
            elif os.path.isdir(f):
                self._save_directory(zip_f, d)

    def set_record(self, record_ident, *step_files):
        if not self._dir_exists:
            ensure_directory(self.cache_dir)
            self._dir_exists = True

        with ZipFile(os.path.join(self.cache_dir, record_ident), mode='w', compression=ZIP_BZIP2) as zip_f:
            for f in step_files:
                if os.path.isfile(f):
                    self._save_file(zip_f, f)
                elif os.path.isdir(f):
                    self._save_directory(zip_f, d)
                else:
                    raise ZCItoolsValueError(f"Cache: file/directory to store ({f}) doesn't exist!")

    def has_record(self, record_ident):
        return self._dir_exists and os.path.isfile(os.path.join(self.cache_dir, record_ident))

    def get_record(self, record_ident, save_location, info=False):
        # Extract record data into given location.
        zip_filename = os.path.join(self.cache_dir, record_ident)
        if os.path.isfile(zip_filename):
            with ZipFile(zip_filename, 'r') as zip_f:
                for z_i in zip_f.infolist():
                    if not z_i.is_dir():
                        save_filename = os.path.join(save_location, z_i.filename)
                        ensure_directory(os.path.dirname(save_filename))
                        if info:
                            print(f"  cache fetch: {zip_filename}[{z_i.filename}] -> {save_filename}")
                        with open(save_filename, 'wb') as save_f:
                            save_f.write(zip_f.read(z_i.filename))
            return True
        return False

    def get_records(self, record_idents, save_location, info=False):
        not_found = []
        for record_ident in record_idents:
            if not self.get_record(record_ident, save_location, info=info):
                not_found.append(record_ident)
        return not_found
