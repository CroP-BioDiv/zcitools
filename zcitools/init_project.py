import os
from zcitools.utils.file_utils import ensure_directory

_readme = """
Tools description:
This is a ZCI tools project! Check git repository https://github.com/CroP-BioDiv/zcitools

Use zcit (zcitools/bin) command line script to work with it.
"""


def init_project(dirname, *args):
    if ensure_directory(dirname):
        with open(os.path.join(dirname, 'README.txt'), 'w') as r:
            if args:
                r.write(f"Project description:\n{' '.join(args)}\n")
            r.write(_readme)
    else:
        print(f'Warning: project {dirname} was not created!')
