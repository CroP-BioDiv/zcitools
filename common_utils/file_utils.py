import os
import errno
import sys
import shutil
import yaml
from zipfile import ZipFile, ZIP_BZIP2

# Settings
settings_defaults = dict(
    ps_viewer='evince',
    image_viewer='gpicview',
    email=None,
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
    with ZipFile(output_filename, 'w', compression=ZIP_BZIP2) as _zip:
        for f in files:
            _zip.write(f)


def unzip_file(zip_filename, into_directory):
    with ZipFile(zip_filename, 'r') as _zip:
        _zip.extractall(into_directory)


def list_zip_files(zip_filename):
    with ZipFile(zip_filename, 'r') as _zip:
        return [z_i.filename for z_i in _zip.infolist() if not z_i.is_dir()]


# Fasta file
def write_fasta(filename, data):
    with open(filename, 'w') as fa:
        for ident, seq in data:
            fa.write(f">{ident}\n{seq}\n")


def read_fasta_identifiers(filename):
    with open(filename, 'r') as fa:
        return [line[1:].strip() for line in fa.readlines() if line[0] == '>']


# Run script in step directory
def run_module_script(module, step):
    project_dir = os.getcwd()
    os.chdir(step.directory)
    module.run(locale=True)
    os.chdir(project_dir)


def set_run_instructions(module, step, files_to_zip, instructions):
    # Copy run script
    script_name = os.path.basename(module.__file__)
    run_f = step.step_file(script_name)
    copy_file(module.__file__, run_f)

    # Instructions
    inst_f = step.step_file('INSTRUCTIONS.txt')
    write_str_in_file(inst_f, instructions.format(step_name=step.directory, script_name=script_name))

    # Zip needed files
    zip_files(step.step_file('calculate.zip'), files_to_zip + [inst_f, run_f])


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
