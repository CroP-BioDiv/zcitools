# Note: importing is done in run() methods to prevent crashes because of not used missing libraries!
import os
from step_project.base_commands import ProjectCommand, NonProjectCommand, CreateStepCommand
from common_utils.exceptions import ZCItoolsValueError


class InitProject(NonProjectCommand):
    _COMMAND = 'init'
    _HELP = "Initialize project in given directory name."

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('dirname', help='Directory name')
        parser.add_argument('-d', '--description', help='Project description text')
        parser.add_argument('-w', '--workflow', help="Set project's workflow")
        parser.add_argument('-p', '--workflow-parameters', default='',
                            help="Workflow's parameters. Format: name1=value1;name2=value2;...")

    def run(self):
        from ..init_project import init_project
        a = self.args
        init_project(self.project, a.dirname, a.description, a.workflow, a.workflow_parameters)


class Unfinish(ProjectCommand):
    _COMMAND = 'unfinish'
    _HELP = "Undo finish step"

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Step name')

    def run(self):
        step = self.project.read_step(self.args.step, update_mode=True)  # Set to be in update mode
        if not step.is_completed():
            print(f"Info: step {self.args.step} is not completed!")
        else:
            step.save_description(step.get_type_description() or dict(), completed=False)


class Finish(ProjectCommand):
    _COMMAND = 'finish'
    _HELP = "Finish step that needed editing."

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Step name')
        parser.add_argument(
            '-f', '--force', action='store_true', help='Force finish even if step is already completed')

    def run(self):
        step = self.project.read_step(self.args.step, update_mode=True)  # Set to be in update mode
        if step.is_completed() and not self.args.force:
            print(f"Info: step {self.args.step} is completed!")
        else:
            from types import SimpleNamespace
            orig_args = SimpleNamespace(**step._step_data['command_args'])
            command_obj = self.project.commands_map[step.get_command()](self.project, orig_args)
            command_obj.finish(step)


class RunAndFinish(ProjectCommand):
    _COMMAND = 'run_and_finish'
    _HELP = "Run locally and finish step"

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Step name')

    def run(self):
        import subprocess
        step_name = self.args.step
        step = self.project.read_step(step_name, no_check=True)  # Set to be in update mode
        if not step.is_completed():
            subprocess.run(['python3', step.run_module_name()], cwd=step_name)
            self.project.run_command_with_args('finish', step_name)


class Clean(ProjectCommand):
    _COMMAND = 'clean'
    _HELP = "Remove not needed step data (cache and processed files)"

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', nargs='*', help='Step name')

    def run(self):
        steps = self.args.step if self.args.step else os.listdir('.')
        for step_name in steps:
            if os.path.isfile(os.path.join(step_name, 'description.yml')):
                step = self.project.read_step(step_name)
                step.remove_cache_files()
                if step.is_completed():
                    step.clean_files()


class CleanCache(ProjectCommand):
    _COMMAND = 'cache'
    _HELP = "Remove cache of given steps"

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', nargs='+', help='Step name')

    def run(self):
        for s in self.args.step:
            step = self.project.read_step(s)
            step.remove_cache_files()


class Show(ProjectCommand):
    _COMMAND = 'show'
    _HELP = "Print step data"

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Step name')
        parser.add_argument('params', nargs='*', help='Additional format option (free format, depends on step type)')

    def run(self):
        step = self.project.read_step(self.args.step)
        step.show_data(params=self.args.params)


class ZipCalculate(ProjectCommand):
    _COMMAND = 'zip_calculate'
    _HELP = "Zip unfinished step's calculate.zip files into calculates_pending.zip"

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step_names', nargs='*', help='Step names')

    def run(self):
        to_zip = []
        for sn in sorted(self.args.step_names or os.listdir()):
            if (step := self.project.read_step_if_in(sn, no_check=True)) and \
               not step.is_completed() and \
               step.is_file('calculate.zip'):
                to_zip.append(step.step_file('calculate.zip'))
        if to_zip:
            from common_utils.file_utils import merge_zip_files
            merge_zip_files('calculates_pending.zip', to_zip, info=True)
        else:
            print('No pending calculations!')


class CopyStepDirectory(CreateStepCommand):
    _COMMAND = 'copy_step'
    _HELP = "Copies step directory into a step. Input step can be from other project"

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step_directory', help='Step directory to copy from')
        parser.add_argument('-t', '--data-type', help='Check step data type')

    def run(self, step_data):
        a = self.args
        if not os.path.isdir(a.step_directory):
            raise ZCItoolsValueError(f"Step directory {a.step_directory} is not an directory!")
        if not (desc_data := read_yaml(os.path.join(a.step_directory, 'description.yml'))):
            raise ZCItoolsValueError(f"Directory {a.step_directory} is not an step!")
        if a.data_type and a.data_type != desc_data['data_type']:
            raise ZCItoolsValueError(
                f"Step {a.step_directory} is not of specified data type ({a.data_type} / {desc_data['data_type']})!")
        assert False, 'ToDo'


#
class Workflow(ProjectCommand):
    _COMMAND = 'workflow'
    _HELP = "Execute workflow command"
    _WF_COMMANDS = ('run', 'graph')

    @staticmethod
    def set_arguments(parser):
        from ..base_workflow import BaseWorkflow
        parser.add_argument('command', choices=BaseWorkflow.all_commands(), help='Command to execute')
        # parser.add_argument('params', nargs='*', help='Additional format option (free format, depends on step type)')

    def run(self):
        self.project.get_workflow().run_command(self.args.command)
