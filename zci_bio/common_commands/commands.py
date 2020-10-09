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
