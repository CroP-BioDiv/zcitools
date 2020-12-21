from step_project.base_commands import ProjectCommand


class Graph(ProjectCommand):
    _COMMAND = 'graph'
    _HELP = "Show project structure as a graph"

    @staticmethod
    def set_arguments(parser):
        pass

    def run(self):
        from .project_graph import create_project_graph
        create_project_graph(self.project)
