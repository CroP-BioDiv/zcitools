from step_project.base_commands import Command, CreateStepCommand, CreateStepFromStepCommand


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
    _HELP = "Fetch SRA summaries for bio projects in given table"
    _INPUT_STEP_DATA_TYPE = 'table'

    def run(self, step_data):
        from .fetch_sra_summaries import fetch_sra_summaries
        return fetch_sra_summaries(step_data, self._input_step())


class NCBISraGroup(CreateStepFromStepCommand):
    _COMMAND = 'ncbi_sra_group'
    _HELP = "Creates table step with aggregated data from NCBI SRA summary data"
    _INPUT_STEP_DATA_TYPE = ('table_grouped', 'table')

    def run(self, step_data):
        from .fetch_sra_summaries import group_sra_data
        return group_sra_data(step_data, self._input_step())


class NCBIAssembliesReport(Command):
    _COMMAND = 'ncbi_assemblies_report'
    _HELP = "Creates excel filestep with aggregated data from NCBI SRA summary data"

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('assemblies', help='Input table step with assemblies')
        parser.add_argument('sra', help='Input table step with sra summaries')
        # Filters
        parser.add_argument('-f', '--from-genome-size', help='From genome size. Format int or float{T|G|M|K}')
        parser.add_argument('-t', '--to-genome-size', help='To genome size. Format int or float{T|G|M|K}')
        parser.add_argument('-F', '--from-date', help='Project from data')
        parser.add_argument('-T', '--to-date', help='Project to data')
        # Output
        parser.add_argument('-o', '--output', default='assemblies.xlsx', help='Output filename')
        parser.add_argument('-p', '--print', action='store_true', help='Print results instead of exporting them')

    def _prev_steps(self):
        return [self.args.assemblies, self.args.sra]

    def run(self):
        from .fetch_sra_summaries import make_report
        ps = self.args
        assem = self.project.read_step(ps.assemblies, check_data_type=('table_grouped', 'table'))
        sra = self.project.read_step(ps.sra, check_data_type=('table_grouped', 'table'))
        return make_report(assem, sra, ps)
