import os
import yaml


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


def write_yaml(data, filename, mode='w'):
    with open(filename, mode, encoding='utf-8') as r:
        r.write(yaml.dump(data, default_flow_style=False))


def read_yaml(filename):
    with open(filename, 'r', encoding='utf-8') as r:
        return yaml.load(r, Loader=yaml.CLoader)


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
