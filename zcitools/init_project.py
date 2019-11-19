import os.path
from zcitools.utils.file_utils import ensure_directory

_readme = """
Tools description:
This is a ZCI tools project! Check git repository https://github.com/CroP-BioDiv/zcitools

Use zcit (zcitools/bin) command line script to work with it.
"""


def init_project(dirname, project_desc):
    if os.path.isfile('project_log.yml'):
        print(f'Warning: init project called on existing project!')
        print(f'Warning: project {dirname} was not created!')

    elif ensure_directory(dirname, check_empty=True):
        # Create empty project.log file
        with open(os.path.join(dirname, 'project_log.yml'), 'w') as r:
            pass
        #
        with open(os.path.join(dirname, 'README.txt'), 'w') as r:
            if project_desc:
                r.write(f'Project description:\n{project_desc}\n')
            r.write(_readme)

    else:
        print(f'Warning: project {dirname} was not created!')
