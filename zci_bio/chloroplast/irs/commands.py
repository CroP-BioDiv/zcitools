from step_project.base_commands import CreateStepFromStepCommand
from common_utils.exceptions import ZCItoolsValueError


class IRsCollectNeededData(CreateStepFromStepCommand):
    _COMMAND = 'analyse_irs_collect_needed_data'
    _COMMAND_GROUP = 'Chloroplast'
    _HELP = "Collect needed data for IRs analysis of a given method"
    _STEP_BASE_NAME = 'AnalyseIRsNeededData'
    _INPUT_STEP_DATA_TYPE = 'table'
    _COMMON_DB_IDENT = ('sequences',)

    @staticmethod
    def set_arguments(parser):
        from .analyse_irs import METHOD_NAMES
        CreateStepFromStepCommand.set_arguments(parser)
        parser.add_argument(
            '-m', '--methods', action='append', choices=METHOD_NAMES,
            help=f'Method(s) to check. Available methods: {", ".join(METHOD_NAMES)}')
        parser.add_argument('method', help='Method to collect for (seqs or ge_seq).')
        parser.add_argument('-s', '--seqs-methods', action='append', help='Methods to collect for seqs type')

    def run(self, step_data):
        from .analyse_irs import analyse_irs_collect_needed_data
        a = self.args
        return analyse_irs_collect_needed_data(
            step_data, self._input_step(no_data_check=True), a.method, a.seqs_methods, self.get_common_db_object())


class IRsAnalyse(CreateStepFromStepCommand):
    _COMMAND = 'analyse_irs'
    _COMMAND_GROUP = 'Chloroplast'
    _HELP = "Analyse IRs of chloroplast genomes. Output is a table with result data."
    _STEP_BASE_NAME = 'AnalyseIRs'
    _INPUT_STEP_DATA_TYPE = 'table'

    @staticmethod
    def set_arguments(parser):
        from .analyse_irs import METHOD_NAMES
        CreateStepFromStepCommand.set_arguments(parser)
        parser.add_argument('seqs_step', help='Step with NCBI sequences.')
        parser.add_argument('-g', '--ge-seq-step', help='Step with GeSeq annotated sequences.')
        parser.add_argument('-c', '--chloe-step', help='Step with Chloe annotated sequences.')
        parser.add_argument(
            '-m', '--methods', action='append', choices=METHOD_NAMES,
            help=f'Method(s) to check. Available methods: {", ".join(METHOD_NAMES)}')
        parser.add_argument('-r', '--ranks', action='append', help='Taxa ranks to group by.')
        parser.add_argument('-n', '--names', action='append', help='Taxa names to group by.')

    def run(self, step_data):
        from .analyse_irs import analyse_irs
        args = self.args
        if not args.methods:
            raise ZCItoolsValueError('No methods defined!')

        seqs = self.project.read_complete_step(
            args.seqs_step, check_data_type=('sequences', 'annotations'), no_check=True)
        ge_seq = self.project.read_complete_step(args.ge_seq_step, check_data_type='annotations', no_check=True) \
            if args.ge_seq_step else None
        chloe = self.project.read_complete_step(args.chloe_step, check_data_type='annotations', no_check=True) \
            if args.chloe_step else None

        return analyse_irs(step_data, self._input_step(no_data_check=True), seqs, ge_seq, chloe,
                           args.methods, args.ranks, args.names)
