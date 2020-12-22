import os
import errno
import sys
import shutil
import csv
import yaml
from zipfile import ZipFile, ZIP_DEFLATED

# Settings
settings_defaults = dict(
    ps_viewer='evince',
    image_viewer='gpicview',
    email=None,
    workflow=None,
    workflow_parameters=None,
)


def get_settings():
    settings = dict((k, None) for k in settings_defaults.keys())
    s = read_yaml('settings.yml')
    if s:
        settings.update(s)
    return settings


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
        source = os.path.relpath(os.path.abspath(source), start=os.path.dirname(dest))
        os.symlink(source, dest)


def copy_file(source, dest):
    shutil.copyfile(source, dest)


# Simple read/write from a file
def write_str_in_file(filename, s):
    with open(filename, 'w') as r:
        r.write(s)


def write_lines_in_file(filename, lines):
    with open(filename, 'w') as r:
        r.write('\n'.join(lines))


def read_file_as_str(filename):
    with open(filename, 'r') as r:
        return r.read()


def append_line_to_file(filename, s):
    with open(filename, 'a') as r:
        r.write(s)
        r.write('\n')


def read_file_as_list(filename):
    with open(filename, 'r') as r:
        return [line.strip() for line in r.readlines()]


# CSV
def write_csv(filename, columns, rows):
    with open(filename, 'w', newline='') as outcsv:
        writer = csv.writer(outcsv, delimiter=';', quotechar='"')
        writer.writerow([n for n, _ in columns])  # Header
        writer.writerows(rows)


def read_csv(filename):
    if os.path.isfile(filename):
        with open(filename, 'r') as incsv:
            reader = csv.reader(incsv, delimiter=';', quotechar='"')
            next(reader)  # Skip header
            return list(reader)
    else:
        return []


# YAML
def write_yaml(data, filename, mode='w'):
    with open(filename, mode, encoding='utf-8') as r:
        r.write(yaml.dump(data, default_flow_style=False))


def read_yaml(filename):
    if os.path.isfile(filename):
        with open(filename, 'r', encoding='utf-8') as r:
            return yaml.load(r, Loader=yaml.CLoader)


def print_yaml(data):
    print(yaml.dump(data, default_flow_style=False))


# ZipFile
def extract_from_zip(zip_f, zip_filename, output_filename):
    with open(output_filename, 'wb') as save_f:
        save_f.write(zip_f.read(zip_filename))


def zip_files(output_filename, files):
    with ZipFile(output_filename, 'w', compression=ZIP_DEFLATED) as _zip:
        for f in files:
            assert os.path.isfile(f), f
            _zip.write(f)


def unzip_file(zip_filename, into_directory):
    with ZipFile(zip_filename, 'r') as _zip:
        _zip.extractall(into_directory)


def list_zip_files(zip_filename):
    with ZipFile(zip_filename, 'r') as _zip:
        return [z_i.filename for z_i in _zip.infolist() if not z_i.is_dir()]


def merge_zip_files(z_name, zip_filenames):
    if len(zip_filenames) == 1:
        copy_file(zip_filenames[0], z_name)
    elif zip_filenames:
        with ZipFile(z_name, 'w', compression=ZIP_DEFLATED) as _zip:
            for fname in zip_filenames:
                zf = ZipFile(fname, 'r')
                for n in zf.namelist():
                    _zip.writestr(n, zf.open(n).read())


# Fasta file
def write_fasta(filename, data):
    with open(filename, 'w') as fa:
        for ident, seq in data:
            fa.write(f">{ident}\n{seq}\n")


def read_fasta_identifiers(filename):
    with open(filename, 'r') as fa:
        return [line[1:].strip() for line in fa.readlines() if line[0] == '>']


# Run script in step directory
def run_module_script(module, step, threads=None, **kwargs):
    project_dir = os.getcwd()
    os.chdir(step.directory)
    module.run(locale=True, threads=threads, **kwargs)
    os.chdir(project_dir)


def set_run_instructions(module, step, files_to_zip, instructions):
    # Copy run script
    script_name = os.path.basename(module.__file__)
    run_f = step.step_file(script_name)
    copy_file(module.__file__, run_f)

    # Instructions
    inst_f = step.step_file('INSTRUCTIONS.txt')
    write_str_in_file(inst_f, instructions.format(step_name=step.directory, script_name=script_name))

    more_to_zip = [inst_f, run_f]

    # Copy exec_utils script
    if run_module := getattr(module, 'exec_utils', None):
        utils_f = step.step_file(os.path.basename(run_module.__file__))
        copy_file(run_module.__file__, utils_f)
        more_to_zip.append(utils_f)

    # Zip needed files
    zip_files(step.step_file('calculate.zip'), files_to_zip + more_to_zip)


#
_ext_to_filetype = dict(
    txt='text', text='text',
    csv='csv',
    xlsx='excel',
    # ...
)


def extension_no_dot(filename):
    _, ext = os.path.splitext(filename)
    return ext[1:] if ext else ext


def filetype_from_ext(filename):
    return _ext_to_filetype.get(extension_no_dot(filename))


def basename_no_ext(filename):
    return os.path.splitext(os.path.basename(filename))[0]


def files_from_args(files_or_dirs, extension):
    for fd in files_or_dirs:
        if os.path.isfile(fd) and fd.endswith(extension):
            yield fd
        elif os.path.isdir(fd):
            for f in sorted(os.listdir(fd)):
                if f.endswith(extension) and os.path.isfile(f := os.path.join(fd, f)):
                    yield f


def find_executable(exe, dir_or_filename=None):
    if dir_or_filename:
        if os.path.isdir(dir_or_filename):
            exe = os.path.join(dir_or_filename, exe)
            if not os.path.isfile(exe):
                raise ValueError(f"File {exe} doesn't exist!")
            return exe
        if os.path.isfile(dir_or_filename):
            return dir_or_filename
        raise ValueError(f"Directory or filename {dir_or_filename} doesn't exist!")
    if not shutil.which(exe):
        raise ValueError(f'No executable {exe}!')
    return exe
