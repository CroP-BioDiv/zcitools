from step_project.base_commands import ProjectCommand, CreateStepCommand, CreateStepFromStepCommand


class NCBIAssembliesList(CreateStepCommand):
    _COMMAND = 'ncbi_assemblies_fetch'
    _COMMAND_GROUP = 'NCBI'
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
    _COMMAND_GROUP = 'NCBI'
    _HELP = "Fetch SRA summaries for bio projects in given table"
    _INPUT_STEP_DATA_TYPE = 'table'

    def run(self, step_data):
        from .fetch_sra_summaries import fetch_sra_summaries
        return fetch_sra_summaries(step_data, self._input_step())


class NCBISraGroup(CreateStepFromStepCommand):
    _COMMAND = 'ncbi_sra_group'
    _COMMAND_GROUP = 'NCBI'
    _HELP = "Creates table step with aggregated data from NCBI SRA summary data"
    _INPUT_STEP_DATA_TYPE = ('table_grouped', 'table')

    def run(self, step_data):
        from .fetch_sra_summaries import group_sra_data
        return group_sra_data(step_data, self._input_step())


class NCBIAssembliesReport(ProjectCommand):
    _COMMAND = 'ncbi_assemblies_report'
    _COMMAND_GROUP = 'NCBI'
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
        parser.add_argument('-m', '--method', help='Filter by assembly method used in the project')
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


# ---------------------------------------------------------
# Chloroplast
# ---------------------------------------------------------
class NCBIChloroplastList(CreateStepCommand):
    _COMMAND = 'ncbi_chloroplast_list'
    _COMMAND_GROUP = 'NCBI'
    _HELP = "Creates table step from NCBI chloroplast genome data"
    _STEP_BASE_NAME = 'chloroplast_list'

    def step_base_name(self):
        return f'{self.args.family}_chloroplast_list' if self.args.family else 'chloroplast_list'

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('-c', '--csv-filename', help='Filename of downloaded csv data.')
        parser.add_argument('-f', '--family', help='Family (ToDo)')
        parser.add_argument('-o', '--outgroup', action='append', help='Outgroup(s)')
        parser.add_argument('-t', '--taxon', action='append', help='General taxon(s)')
        parser.add_argument('--max-update-date', help='Set max update date. ISO format (yyyy-mm-dd)')
        parser.add_argument('--remove-irl', action='store_true', help='Remove IRL clade species')
        parser.add_argument('-P', '--fetch-plastids', action='store_true', help='Fetch sequences declared as plastids.')

    def run(self, step_data):
        from .fetch_genome_assemblies import fetch_chloroplast_list
        return fetch_chloroplast_list(self.project, step_data, self.args)


class NCBIChloroplastListToExcel(ProjectCommand):
    _COMMAND = 'ncbi_chloroplast_list_to_excel'
    _COMMAND_GROUP = 'NCBI'
    _HELP = "Exports list of sequences into excel"

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('-s', '--step-name', default='01_chloroplast_list', help='Table step name with list of sequences.')
        parser.add_argument('-o', '--output-filename', default='chloroplast_list.xls', help='Output excel filename')
        parser.add_argument('-c', '--num-columns', default=2, help='Number of columns')

    def run(self):
        from .fetch_genome_assemblies import export_chloroplast_list
        a = self.args
        step = self.project.read_step(a.step_name, check_data_type='table', no_check=True)
        return export_chloroplast_list(step, a.output_filename, a.num_columns)
