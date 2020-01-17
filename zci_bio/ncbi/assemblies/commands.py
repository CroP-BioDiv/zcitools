from step_project.base_commands import CreateStepCommand, CreateStepFromStepCommand


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


class NCBISraSummary(CreateStepFromStepCommand):
    _COMMAND = 'ncbi_sra_summary'
    _HELP = "Creates tables step from table of NCBI bio projects"
    _INPUT_STEP_DATA_TYPE = 'table'

    def run(self, step_data):
        from .fetch_sra_summaries import fetch_sra_summaries
        return fetch_sra_summaries(step_data, self._input_step())


class NCBISraGroup(CreateStepFromStepCommand):
    _COMMAND = 'ncbi_sra_group'
    _HELP = "Creates tables step from table of NCBI SRA data"
    _INPUT_STEP_DATA_TYPE = ('table_grouped', 'table')

    def run(self, step_data):
        from .fetch_sra_summaries import group_sra_data
        return group_sra_data(step_data, self._input_step())
