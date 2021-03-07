import os.path
from common_utils.file_utils import settings_defaults, ensure_directory, write_yaml
from common_utils.exceptions import ZCItoolsValueError

_readme = """
Tools description:
This is a ZCI tools project! Check git repository https://github.com/CroP-BioDiv/zcitools

Use zcit (zcitools/bin) command line script to work with it.
"""

_wf_readme = """
Note: project is set to be of workflow {workflow}. Use command:
zcit workflow <command> <args>
"""


def init_project(project, dirname, project_desc, workflow, workflow_parameters):
    if os.path.isfile('project_log.yml'):
        print(f'Warning: init project called on existing project!')
        print(f'Warning: project {dirname} was not created!')

    elif ensure_directory(dirname, check_empty=True):
        # Add setting file
        settings = dict(settings_defaults)
        if workflow:
            if workflow_parameters:
                w_pars = dict(x.split('=') for x in workflow_parameters.split(';'))
            else:
                w_pars = dict()
            wf_cls = project.get_workflow_cls(workflow)
            if (not_in := [p for p in wf_cls.required_parameters() if p not in w_pars]):
                raise ZCItoolsValueError(f"Workflow's parameters not specified: {', '.join(not_in)}!")

            settings['workflow'] = workflow
            settings['workflow_parameters'] = wf_cls.format_parameters(w_pars)
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
