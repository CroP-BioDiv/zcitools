# Note: importing is done in run() methods to prevent crashes because of not used missing libraries!
from types import SimpleNamespace
from step_project.base_commands import ProjectCommand, NonProjectCommand


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
            orig_args = SimpleNamespace(**step._step_data['command_args'])
            command_obj = self.project.commands_map[step.get_command()](self.project, orig_args)
            command_obj.finish(step)


class Clean(ProjectCommand):
    _COMMAND = 'clean'
    _HELP = "Remove not needed step data (cache and processed files)"

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', nargs='*', help='Step name')

    def run(self):
        import os
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
