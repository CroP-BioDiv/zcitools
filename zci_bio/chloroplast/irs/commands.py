from step_project.base_commands import CreateStepFromStepCommand
# ProjectCommand, NonProjectCommand, CreateStepCommand,


class IRsAnalyse(CreateStepFromStepCommand):
    _COMMAND = 'analyse_irs'
    _HELP = "Analyse IRs of chloroplast genomes. Output is a table with result data."
    _COMMAND_GROUP = 'Chloroplast'
    _STEP_BASE_NAME = 'AnalyseIRs'
    _INPUT_STEP_DATA_TYPE = 'annotations'

    @staticmethod
    def set_arguments(parser):
        from .analyse_irs import AnalyseIRs
        CreateStepFromStepCommand.set_arguments(parser)
        choices = AnalyseIRs.get_method_names()
        parser.add_argument(
            '-m', '--methods',
            help=f'Method(s) to check. Format "method1;method2;...". Available methods: {", ".join(choices)}')
        parser.add_argument(
            '-a', '--export-all', action='store_true',
            help='Excel output to contain all data not only problematic.')

    def run(self, step_data):
        from .analyse_irs import analyse_irs
        args = self.args
        return analyse_irs(step_data,
                           self._input_step(no_data_check=True),
                           args.methods.split(';'),
                           args.export_all)
