# Note: importing is done in run() methods to prevent crashes because of not used missing libraries!
from .base import _Command


class InitProject(_Command):
    _COMMAND = 'init'
    _HELP = "Initialize project in given directory name."

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('dirname', help='Directory name')
        parser.add_argument('-d', '--description', help='Project description text')

    def run(self):
        from ..processing.init_project import init_project
        init_project(self.args.dirname, self.args.description)


class CleanCache(_Command):
    _COMMAND = 'cache'
    _HELP = "Remove cache of given steps"

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', nargs='+', help='Step name')

    def run(self):
        from ..steps import read_step
        for s in self.args.step:
            step = read_step(s)
            step.remove_cache_files()


class Show(_Command):
    _COMMAND = 'show'
    _HELP = "Print step(s) data"

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Step name')
        parser.add_argument('params', nargs='*', help='Additional format option (free format, depends on step type)')

    def run(self):
        from ..steps import read_step
        step = read_step(self.args.step)
        step.show_data(params=self.args.params)


class Graph(_Command):
    _COMMAND = 'graph'
    _HELP = "Show project structure as a graph"

    @staticmethod
    def set_arguments(parser):
        pass

    def run(self):
        from ..processing.project_graph import create_graph
        create_graph()
