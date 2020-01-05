from step_project.base.commands import Command


class Graph(Command):
    _COMMAND = 'graph'
    _HELP = "Show project structure as a graph"

    @staticmethod
    def set_arguments(parser):
        pass

    def run(self):
        from .project_graph import create_graph
        create_graph(self.project)
