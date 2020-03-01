import os
from zipfile import ZipFile, ZIP_BZIP2
from .file_utils import ensure_directory
from .exceptions import ZCItoolsValueError

"""
Common DB stores data used in more projects.
DB object point to specific directory, and stores records.
Record is named by step data identifier, which depends on step type.
Record is implemented as zip file. It can contain any kind of data: any number of files and/or directories.
"""

_ZCI_COMMON_DB_DIR = os.environ.get('ZCI_COMMON_DB')
SEQUENCE_DBS_RELATIVE_DIR = 'sequence_dbs'
_SEQUENCE_DBS_DIR = os.path.join(_ZCI_COMMON_DB_DIR, SEQUENCE_DBS_RELATIVE_DIR) if _ZCI_COMMON_DB_DIR else None


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

    @staticmethod
    def get_zci_sequence_dbs():
        if _SEQUENCE_DBS_DIR:
            return [f for f in os.listdir(_SEQUENCE_DBS_DIR) if os.path.isdir(os.path.join(_SEQUENCE_DBS_DIR, f))]
        return []

    @staticmethod
    def get_check_zci_sequence_db(db):
        dbs = CommonDB.get_zci_sequence_dbs()
        if db not in dbs:
            raise ZCItoolsValueError(f"Database {db} doesn't exist. Possible values: {', '.join(dbs)}!")

    def get_sequence_db(self):
        if len(self._idents) >= 2 and self._idents[0] == SEQUENCE_DBS_RELATIVE_DIR:
            return self._idents[1]

    #
    def get_relative_db(self, *path):
        db_dir = os.path.normpath(os.path.join(self.db_dir, *path))
        return CommonDB(db_dir, base_dir=self._base_dir)

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

    def get_record_filename(self, record_ident):
        return os.path.join(self.db_dir, record_ident)

    def ensure_location(self):
        ensure_directory(self.db_dir)

    def _check_set_record(self, record_ident, force=False):
        # Returns reccord filename to set recor into, or None if data already exists
        rec_filename = self.get_record_filename(record_ident)
        if not force and os.path.isfile(rec_filename):
            # print(f"  CommonDB: record {rec_filename[self._strip_length:]} already exists!")
            return

        if not self._dir_exists:
            ensure_directory(self.db_dir)
            self._dir_exists = True
        return rec_filename

    def set_record(self, record_ident, *step_files, force=False):
        rec_filename = self._check_set_record(record_ident, force=force)
        if rec_filename:
            with ZipFile(rec_filename, mode='w', compression=ZIP_BZIP2) as zip_f:
                for f in step_files:
                    if os.path.isfile(f):
                        self._save_file(zip_f, f)
                    elif os.path.isdir(f):
                        self._save_directory(zip_f, d)
                    else:
                        raise ZCItoolsValueError(f"Common DB: file/directory to store ({f}) doesn't exist!")

    def set_record_from_stream(self, record_ident, data, arcname, force=False):
        rec_filename = self._check_set_record(record_ident, force=force)
        if rec_filename:
            with ZipFile(rec_filename, mode='w', compression=ZIP_BZIP2) as zip_f:
                zip_f.writestr(arcname, data)

    def has_record(self, record_ident):
        return self._dir_exists and os.path.isfile(self.get_record_filename(record_ident))

    def get_record(self, record_ident, save_location, info=False):
        # Extract record data into given location.
        zip_filename = self.get_record_filename(record_ident)
        if os.path.isfile(zip_filename):
            with ZipFile(zip_filename, 'r') as zip_f:
                for z_i in zip_f.infolist():
                    if not z_i.is_dir():
                        save_filename = os.path.join(save_location, z_i.filename)
                        ensure_directory(os.path.dirname(save_filename))
                        if info:
                            print(f"  CommonDB fetch: {zip_filename[self._strip_length:]}[{z_i.filename}] " +
                                  f"-> {save_filename}")
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
