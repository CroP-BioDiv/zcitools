import os
import sys
import shutil
import yaml
from zipfile import ZipFile, ZIP_BZIP2


# Files and directories
def ensure_directory(d, check_empty=False):
    if d == '':  # Current directory
        return
    if not os.path.exists(d):
        os.makedirs(d)
        return True
    if os.path.isfile(d):
        raise CreateDirectoryError(f"Can't create directory, file exists with same name ({d})!")
    if check_empty and os.listdir(d):
        print(f'Warning: directory {d} exists but is not empty!')
        return False
    return True


def remove_directory(_dir, create=False):
    # assert not _dir.endswith("/") or _dir.endswith("\\"), _dir
    _dir = os.path.normpath(_dir)

    if os.path.isdir(_dir):
        if sys.platform == "win32":
            temp_path = _dir + "_"

            if os.path.exists(temp_path):
                remove_and_create(temp_path, create=False)

            try:
                os.renames(_dir, temp_path)
            except OSError as exception:
                if exception.errno != errno.ENOENT:
                    raise
            else:
                shutil.rmtree(temp_path)
        else:
            shutil.rmtree(_dir)

    if create:
        os.makedirs(_dir)


def silent_remove_file(filename):
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno != errno.ENOENT:
            raise


def silent_remove(filename):
    if os.path.isfile(filename):
        silent_remove_file(filename)
    elif os.path.isdir(filename):
        remove_directory(filename)


def link_file(source, dest):
    if sys.platform == "win32":
        shutil.copyfile(source, dest)
    else:
        os.symlink(source, dest)


def copy_file(source, dest):
    shutil.copyfile(source, dest)


# Simple read/write from a file
def write_str_in_file(filename, s):
    with open(filename, 'w') as r:
        r.write(s)


def read_file_as_str(filename):
    with open(filename, 'r') as r:
        return r.read()


def read_file_as_list(filename):
    with open(filename, 'r') as r:
        return [line.strip() for line in r.readlines()]


# YAML
def write_yaml(data, filename, mode='w'):
    with open(filename, mode, encoding='utf-8') as r:
        r.write(yaml.dump(data, default_flow_style=False))


def read_yaml(filename):
    if os.path.isfile(filename):
        with open(filename, 'r', encoding='utf-8') as r:
            return yaml.load(r, Loader=yaml.CLoader)


# ZipFile
def extract_from_zip(zip_f, zip_filename, output_filename):
    with open(output_filename, 'wb') as save_f:
        save_f.write(zip_f.read(zip_filename))


def zip_files(output_filename, files):
    with ZipFile(output_filename, 'w', compression=ZIP_BZIP2) as output:
        for f in files:
            output.write(f)


#
_ext_to_filetype = dict(
    txt='text', text='text',
    csv='csv',
    xlsx='excel',
    # ...
)


def filetype_from_ext(filename):
    _, ext = os.path.splitext(filename)
    if ext:
        ext = ext[1:]
    return _ext_to_filetype.get(ext)
