from step_project.base_commands import ProjectCommand, NonProjectCommand, CreateStepCommand, CreateStepFromStepCommand
from common_utils.exceptions import ZCItoolsValueError
from .constants import DEFAULT_KEEP_OFFSET


class ChloroplastAnalyse(CreateStepFromStepCommand):
    _COMMAND = 'analyse_chloroplast'
    _HELP = "Analyse chloroplast genomes. Output is a table with result data."
    _COMMAND_GROUP = 'Chloroplast'
    _STEP_BASE_NAME = 'AnalyseChloroplast'
    _INPUT_STEP_DATA_TYPE = 'annotations'

    def run(self, step_data):
        from .analyse import analyse_genomes
        return analyse_genomes(step_data, self._input_step(no_data_check=True))


class ChloroplastFixByAnalyse(CreateStepFromStepCommand):
    _COMMAND = 'fix_by_analysis'
    _HELP = "Fix chloroplast genome by done analyse. Output is a sequences step."
    _COMMAND_GROUP = 'Chloroplast'
    _INPUT_STEP_DATA_TYPE = 'table'
    _COMMON_DB_IDENT = ('sequences',)

    @staticmethod
    def set_arguments(parser):
        # Note: method than step
        parser.add_argument('method', help='Fix method. Options: parts, trnF-GAA. Only first character is needed.')
        parser.add_argument('subset', choices=('all', 'sum', 'ge_seq', 'ncbi'),
                            help='Subset to work with')
        CreateStepFromStepCommand.set_arguments(parser)
        parser.add_argument('-o', '--keep-offset', default=DEFAULT_KEEP_OFFSET, type=int,
                            help='Do not rotate genome if offset is less than given value.')

    def step_base_name(self):
        m = self.args.method[0].lower()
        if m == 'p':
            n = self._format_step_name('FixByParts')
        # elif m == 'h':
        #     n = self._format_step_name('FixByTrnH-GUG')
        elif m == 'f':
            n = self._format_step_name('FixByTrnF-GAA')
        else:
            raise ZCItoolsValueError(f'Not known method {self.args.method}!')
        # Add offset or not?
        o = self.args.keep_offset
        return f'{n}_{o}' if o != DEFAULT_KEEP_OFFSET else n

    def run(self, step_data):
        m = self.args.method[0].lower()
        if m == 'p':
            from .fix_by_analysis import fix_by_parts as fix_method
        # elif m == 'h':
        #     from .fix_by_analysis import fix_by_trnH_GUG as fix_method
        elif m == 'f':
            from .fix_by_analysis import fix_by_trnF_GAA as fix_method
        else:
            raise ZCItoolsValueError(f'Not known method {self.args.method}!')
        #
        return fix_method(step_data,
                          self._input_step(no_data_check=True),
                          self.args.subset, self.args.keep_offset,
                          self.get_common_db_object())


# Normalization results
class ChloroplastNormalizationResult(CreateStepCommand):
    _COMMAND = 'normalization_result'
    _HELP = "Analyses phylogeny results of chloroplast normalization."
    _COMMAND_GROUP = 'Chloroplast'
    _STEP_BASE_NAME = 'normalization_result'

    def run(self, step_data):
        from .normalization_result import NormalizationResult
        return NormalizationResult(self.project, '.').run(step_data)


class ChloroplastNormalizationResultGraph(ProjectCommand):
    _COMMAND = 'normalization_result_graph'
    _HELP = "Create graph from results of chloroplast normalization."
    _COMMAND_GROUP = 'Chloroplast'

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Input step')

    def run(self):
        from .normalization_result import NormalizationResult
        step = self.project.read_step(self.args.step, check_data_type='table', no_check=True)
        return NormalizationResult(self.project, '.').create_graph(step, show=True)


class ChloroplastNormalizationResultGraphJoin(NonProjectCommand):
    _COMMAND = 'normalization_result_graph_join'
    _HELP = "Create chloroplast normalization result graphs and join them from"
    _COMMAND_GROUP = 'Chloroplast'

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('steps', nargs='+', help='Input steps')
        parser.add_argument('-t', '--two-columns', action='store_true', help='Graphs in two columns')

    def run(self):
        from os.path import sep
        from .normalization_result import NormalizationResult
        steps = [self.project.read_step(s.split(sep), check_data_type='table', no_check=True, outside_of_project=True)
                 for s in self.args.steps]
        return NormalizationResult.create_graphs(self.project, steps, self.args.two_columns, show=True)


# Stats
class ChloroplastNormalizationStatByTaxonomy(ProjectCommand):
    _COMMAND = 'normalization_stat_by_taxonomy'
    _HELP = "Create report of chloroplast normalization by grouping data in taxonomy lavels"
    _COMMAND_GROUP = 'Chloroplast'

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Input chloroplast analyses step')
        parser.add_argument('-r', '--ranks', help='Taxa ranks to group by. Format "rank1:rank2".')
        parser.add_argument('-n', '--names', help='Taxa names to group by. Format "name1:name2".')
        parser.add_argument('-m', '--minimum-sequences', type=int, help='Minimum sequences to report')
        parser.add_argument('-o', '--output-excel', help='Excel output filename')
        parser.add_argument('-c', '--print-output', action='store_true',
                            help='Print result into console')
        parser.add_argument('-p', '--pie-chart-name', action='append',
                            help='Taxa name to make pie chart. More charts are possible')
        parser.add_argument('-P', '--pie-chart-all-names', action='store_true',
                            help='Make pie charts of all taxa names.')
        parser.add_argument('-M', '--merge-pie-charts', action='store_true',
                            help='Make one figure with all pie charts')

    def run(self):
        from .stats import statistics_by_taxa
        args = self.args
        taxa_names = args.names.split(':') if args.names else []
        pie_chart_names = list(args.pie_chart_name) if args.pie_chart_name else []
        if args.pie_chart_all_names:
            pie_chart_names.extend(taxa_names)
        statistics_by_taxa(self.project,
                           self.project.read_step(args.step, check_data_type='table'),
                           args.ranks.split(':') if args.ranks else None,
                           taxa_names,
                           args.minimum_sequences,
                           args.output_excel,
                           pie_chart_names, args.merge_pie_charts,
                           args.print_output)
