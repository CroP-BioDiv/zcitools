import os
from io import StringIO
from zipfile import ZipFile, ZIP_DEFLATED
from .file_utils import ensure_directory, silent_remove
from .exceptions import ZCItoolsValueError

"""
Common DB stores data used in more projects.
DB object point to specific directory, and stores records.
Record is named by step data identifier, which depends on step type.
Record is implemented as zip file. It can contain any kind of data: any number of files and/or directories.
"""

_ZCI_COMMON_DB_DIR = os.environ.get('ZCI_COMMON_DB')


class CommonDB:
    def __init__(self, db_dir, base_dir=None):
        assert not base_dir or db_dir.startswith(base_dir), (db_dir, base_dir)
        self.db_dir = db_dir
        self._idents = tuple(db_dir[(len(base_dir) + 1):].split(os.path.sep)) if base_dir else tuple()
        self._base_dir = base_dir
        #
        self._strip_length = (len(base_dir) + 1) if base_dir else 0
        self._dir_exists = os.path.exists(db_dir)
        if self._dir_exists and os.path.isfile(db_dir):
            raise ZCItoolsValueError(f'Common DB location ({db_dir}) is a file!')

    @staticmethod
    def get_zci_db(idents):
        return CommonDB(os.path.join(_ZCI_COMMON_DB_DIR, *idents), base_dir=_ZCI_COMMON_DB_DIR)

    #
    def get_relative_db(self, *path):
        db_dir = os.path.normpath(os.path.join(self.db_dir, *path))
        return CommonDB(db_dir, base_dir=self._base_dir)

    def _save_file(self, zip_f, f, remove_directories=False):
        # Note: step filenames are relative to project main directory. Zip filenames are stored without step_name.
        parts = f.split(os.path.sep)
        if remove_directories:
            zip_name = parts[-1]
        else:
            zip_name = f if len(parts) == 1 else os.path.join(*(f.split(os.path.sep)[1:]))
        zip_f.write(f, arcname=zip_name)

    def _save_directory(self, zip_f, d, remove_directories=False):
        # ToDo: is it needed to save directory into the zip?
        for f in os.listdir(d):
            f = os.path.join(d, f)
            if os.path.isfile(f):
                self._save_file(zip_f, f, remove_directories=remove_directories)
            elif os.path.isdir(f):
                self._save_directory(zip_f, d, remove_directories=remove_directories)

    def get_record_filename(self, record_ident):
        if isinstance(record_ident, (list, tuple)):
            return os.path.join(self.db_dir, *record_ident)
        return os.path.join(self.db_dir, record_ident)

    def ensure_location(self):
        ensure_directory(self.db_dir)

    def _check_set_record(self, record_ident, force=False, info=False):
        # Returns reccord filename to set recor into, or None if data already exists
        rec_filename = self.get_record_filename(record_ident)
        if not force and os.path.isfile(rec_filename):
            if info:
                print(f"  CommonDB: record {rec_filename[self._strip_length:]} already exists!")
            return

        if not self._dir_exists:
            ensure_directory(self.db_dir)
            self._dir_exists = True
        return rec_filename

    def set_record(self, record_ident, *step_files, force=False, info=False, remove_directories=False):
        rec_filename = self._check_set_record(record_ident, force=force, info=info)
        if rec_filename:
            with ZipFile(rec_filename, mode='w', compression=ZIP_DEFLATED) as zip_f:
                for f in step_files:
                    if info:
                        print(f"  CommonDB zipping: {f} -> {record_ident}")
                    if os.path.isfile(f):
                        self._save_file(zip_f, f, remove_directories=remove_directories)
                    elif os.path.isdir(f):
                        self._save_directory(zip_f, d, remove_directories=remove_directories)
                    else:
                        raise ZCItoolsValueError(f"Common DB: file/directory to store ({f}) doesn't exist!")

    def set_record_from_stream(self, record_ident, data, arcname, force=False, info=False):
        rec_filename = self._check_set_record(record_ident, force=force, info=info)
        if rec_filename:
            with ZipFile(rec_filename, mode='w', compression=ZIP_DEFLATED) as zip_f:
                zip_f.writestr(arcname, data)

    # Get record data
    def has_record(self, record_ident):
        return self._dir_exists and os.path.isfile(self.get_record_filename(record_ident))

    def has_records(self, record_idents):
        if self._dir_exists:
            return [ri for ri in record_idents if os.path.isfile(self.get_record_filename(ri))]
        return []

    def get_all_record_ident(self, startswith=None):
        n = len(self.db_dir) + 1
        for root, subdirs, files in os.walk(self.db_dir):
            if (r := root[n:]):
                rel_dir = os.path.split(r)
            else:
                rel_dir = tuple()
            if startswith:
                yield from (rel_dir + (f,) for f in files if f.startswith(startswith))
            else:
                yield from (rel_dir + (f,) for f in files)

    def _get_zip_data(self, record_ident):
        # Returns tuple (filename, data) of a record.
        zip_filename = self.get_record_filename(record_ident)
        if os.path.isfile(zip_filename):
            with ZipFile(zip_filename, 'r') as zip_f:
                for z_i in zip_f.infolist():
                    if not z_i.is_dir():
                        yield zip_filename, z_i.filename, zip_f.read(z_i.filename)

    def get_record(self, record_ident, save_location, info=False):
        # Extract record data into given location.
        ret = None
        for zip_filename, filename, data in self._get_zip_data(record_ident):
            save_filename = os.path.join(save_location, filename)
            ensure_directory(os.path.dirname(save_filename))  # Note: filename can have path!
            if info:
                print(f"  CommonDB fetch: {zip_filename[self._strip_length:]}[{filename}] -> {save_filename}")
            with open(save_filename, 'wb') as save_f:
                save_f.write(data)
            ret = save_filename
        return ret

    def get_records(self, record_idents, save_location, info=False):
        not_found = []
        for record_ident in record_idents:
            if not self.get_record(record_ident, save_location, info=info):
                not_found.append(record_ident)
        return not_found

    def get_record_data(self, record_ident):
        for _, _, data in self._get_zip_data(record_ident):
            return data

    def get_record_stringIO(self, record_ident):
        if (d := self.get_record_data(record_ident)):
            return StringIO(d.decode('utf-8'))

    def get_record_str(self, record_ident):
        if (d := self.get_record_data(record_ident)):
            return d.decode()

    def get_one_file(self, record_ident):
        files = list(self._get_zip_data(record_ident))
        if len(files) != 1:
            print(f"Warning: record {record_ident} has zero or more files!")
            return
        return files[0]

    #
    def remove_record(self, record_ident):
        silent_remove(self.get_record_filename(record_ident))
