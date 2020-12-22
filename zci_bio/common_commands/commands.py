from step_project.base_commands import CreateStepFromStepCommand
from common_utils.exceptions import ZCItoolsValueError


class FetchCommonDB(CreateStepFromStepCommand):
    _COMMAND = 'fetch_common_db'
    _HELP = "Creates sequences or annotations step. Mandatory argument is a table step."
    # _STEP_BASE_NAME = 'seqs'
    _INPUT_STEP_DATA_TYPE = 'table'

    @staticmethod
    def set_arguments(parser):
        CreateStepFromStepCommand.set_arguments(parser)               # step
        parser.add_argument('step_type', choices=('sequences', 'annotations', 'images'), help='Step type')
        parser.add_argument('idents', nargs='+', help='Database inside sequence database to use.')

    def step_base_name(self):
        return '_'.join(self.args.idents)

    def common_db_identifier(self):
        return tuple(self.args.idents)

    def run(self, step_data):
        from .fetch import fetch_common_db_data
        return fetch_common_db_data(step_data, self._input_step(), self.args.step_type, self.get_common_db_object())


class SequencesSubset(CreateStepFromStepCommand):
    _COMMAND = 'seq_subset'
    _HELP = "Creates sequences or annotations step, that contains subset of sequences of given step."
    _STEP_BASE_NAME = 'subset'
    _INPUT_STEP_DATA_TYPE = ('sequences', 'annotations')
    _COMMAND_GROUP = 'Bio'

    @staticmethod
    def set_arguments(parser):
        CreateStepFromStepCommand.set_arguments(parser)
        # parser.add_argument('-f', '--from-file', help='File with sequence identifiers to extract')
        parser.add_argument('-s', '--seq-idents', action='append', help='Sequences to take.')
        parser.add_argument('-S', '--seq-idents-re', action='append', help='Sequences to take (regexp)')
        parser.add_argument('-w', '--without-seq-idents', action='append', help='Sequences to omit')
        parser.add_argument('-W', '--without-seq-idents-re', action='append', help='Sequences to omit (regexp)')
        parser.add_argument(
            '--analyses-with-irs',
            help='Chloroplast analyses step where to check for chloroplast sequences that have IRs')

    def run(self, step_data):
        from .fetch import create_subset
        return create_subset(step_data, self._input_step(), self.args)
