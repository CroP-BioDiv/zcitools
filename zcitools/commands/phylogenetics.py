from zcitools.base.commands import CreateStepFromStepCommand
from zcitools.utils.exceptions import ZCItoolsValueError

# Check base.py for description


class RAxML(CreateStepFromStepCommand):
    _COMMAND = 'raxml'
    _HELP = "Run RAxML on alignment(s)"
    _STEP_BASE_NAME = 'RAxML'
    _INPUT_STEP_DATA_TYPE = ('alignment', 'alignments')
    # Note: no global caching

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Input sequences step')
        parser.add_argument('-r', '--run', action='store_true', help='Run RAxML locale')

    def run(self, step_data):
        from ..processing.phylogenetics.raxml import create_raxml_data
        return create_raxml_data(step_data, self._input_step(), self.get_cache_object(), self.args.run)

    def finish(self, step_obj):
        from ..processing.phylogenetics.raxml import finish_raxml_data
        finish_raxml_data(step_obj, self.get_cache_object())


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
        from ..processing.phylogenetics.mr_bayes import create_mr_bayes_data
        return create_mr_bayes_data(step_data, self._input_step(), self.get_cache_object(), self.args.run)

    def finish(self, step_obj):
        from ..processing.phylogenetics.mr_bayes import finish_mr_bayes_data
        finish_mr_bayes_data(step_obj, self.get_cache_object())
