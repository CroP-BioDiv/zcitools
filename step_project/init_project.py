import os.path
from common_utils.file_utils import settings_defaults, ensure_directory, write_yaml

_readme = """
Tools description:
This is a ZCI tools project! Check git repository https://github.com/CroP-BioDiv/zcitools

Use zcit (zcitools/bin) command line script to work with it.
"""

_wf_readme = """
Note: project is set to be of workflow {workflow}. Use command:
zcit workflow <command> <args>
"""


def init_project(dirname, project_desc, workflow):
    if os.path.isfile('project_log.yml'):
        print(f'Warning: init project called on existing project!')
        print(f'Warning: project {dirname} was not created!')

    elif ensure_directory(dirname, check_empty=True):
        # Add setting file
        settings = dict(settings_defaults)
        if workflow:
            settings['workflow'] = workflow
        write_yaml(settings, os.path.join(dirname, 'settings.yml'))

        # Create empty project.log file
        with open(os.path.join(dirname, 'project_log.yml'), 'w') as r:
            pass

        # Set README.txt file
        with open(os.path.join(dirname, 'README.txt'), 'w') as r:
            if project_desc:
                r.write(f'Project description:\n{project_desc}\n')
            r.write(_readme)
            if workflow:
                r.write(_wf_readme.format(workflow=workflow))

    else:
        print(f'Warning: project {dirname} was not created!')
