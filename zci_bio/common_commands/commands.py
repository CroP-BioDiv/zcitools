from step_project.base_commands import CreateStepFromStepCommand
from common_utils.exceptions import ZCItoolsValueError


class FetchCommonDB(CreateStepFromStepCommand):
    _COMMAND = 'fetch_common_db'
    _HELP = "Creates sequences or annotations step. Mandatory argument is a table step."
    # _STEP_BASE_NAME = 'seqs'
    _INPUT_STEP_DATA_TYPE = 'table'

    @staticmethod
    def set_arguments(parser):
        CreateStepFromStepCommand._sequence_db_set_arguments(parser)  # sequence_db
        CreateStepFromStepCommand.set_arguments(parser)               # step
        parser.add_argument('db_ident', help='Database inside sequence database to use.')
        parser.add_argument('step_type', choices=('sequences', 'annotations', 'images'), help='Step type')

    def step_base_name(self):
        return f"{self.args.sequence_db}_{self.args.db_ident}"

    def common_db_identifier(self):
        return self.sequence_db_identifier(self.args.sequence_db, self.args.db_ident)

    def run(self, step_data):
        from .fetch import fetch_common_db_data
        return fetch_common_db_data(step_data, self._input_step(), self.args.step_type, self.get_common_db_object())
