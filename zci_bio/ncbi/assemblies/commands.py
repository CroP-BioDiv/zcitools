from step_project.base.commands import CreateStepCommand


class NCBIAssembliesList(CreateStepCommand):
    _COMMAND = 'ncbi_assemblies_fetch'
    _HELP = "Creates table step from NCBI genome assemblies data"
    _STEP_BASE_NAME = 'ncbi_assemblies'

    @staticmethod
    def set_arguments(parser):
        pass

    def run(self, step_data):
        from .fetch_genome_assemblies import fetch_genome_assemblies
        return fetch_genome_assemblies(self.project, step_data)

    def finish(self, step_obj):
        from .fetch_genome_assemblies import finish_fetch_genome_assemblies
        finish_fetch_genome_assemblies(step_obj)
