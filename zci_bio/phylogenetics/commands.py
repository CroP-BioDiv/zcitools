from step_project.base_commands import CreateStepFromStepCommand
from common_utils.exceptions import ZCItoolsValueError


class RAxML(CreateStepFromStepCommand):
    _COMMAND = 'raxml'
    _HELP = "Run RAxML on alignment(s)"
    _STEP_BASE_NAME = 'RAxML'
    _INPUT_STEP_DATA_TYPE = ('alignment', 'alignments')
    # Note: no global caching

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Input sequences step')
        parser.add_argument('-p', '--no-partitions', action='store_true', help='Do not set partitions')
        parser.add_argument('-r', '--run', action='store_true', help='Run RAxML locale')

    def run(self, step_data):
        from .raxml import create_raxml_data
        return create_raxml_data(step_data, self._input_step(), not self.args.no_partitions, self.args.run)

    def finish(self, step_obj):
        from .raxml import finish_raxml_data
        finish_raxml_data(step_obj)


class MrBayes(CreateStepFromStepCommand):
    _COMMAND = 'mr_bayes'
    _HELP = "Run MrBayes on alignment(s)"
    _STEP_BASE_NAME = 'MrBayes'
    _INPUT_STEP_DATA_TYPE = ('alignment', 'alignments')
    # Note: no global caching

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Input sequences step')
        parser.add_argument('-r', '--run', action='store_true', help='Run MrBayes locale')

    def run(self, step_data):
        from .mr_bayes import create_mr_bayes_data
        return create_mr_bayes_data(step_data, self._input_step(), self.args.run)

    def finish(self, step_obj):
        from .mr_bayes import finish_mr_bayes_data
        finish_mr_bayes_data(step_obj)
