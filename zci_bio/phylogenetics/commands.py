from step_project.base_commands import CreateStepFromStepCommand
from common_utils.exceptions import ZCItoolsValueError


class RAxML(CreateStepFromStepCommand):
    _COMMAND = 'raxml'
    _HELP = "Run RAxML on alignment(s)"
    _STEP_BASE_NAME = 'RAxML'
    _INPUT_STEP_DATA_TYPE = ('alignment', 'alignments')

    @staticmethod
    def set_arguments(parser):
        CreateStepFromStepCommand.set_arguments(parser)
        parser.add_argument('-a', '--annotations-step', help='Annotations step to use for partitions')
        parser.add_argument('-p', '--no-partitions', action='store_true', help='Do not set partitions')
        parser.add_argument('-r', '--run', action='store_true', help='Run RAxML locale')

    def step_base_name(self):
        return self._format_step_name(f"{self._STEP_BASE_NAME}_{'N' if self.args.no_partitions else 'P'}")

    def _partitions(self):
        from .partitions import Partitions
        a = self.args
        if a.no_partitions:
            annotations_step = None
        elif a.annotations_step:
            annotations_step = self.project.read_step(a.annotations_step, check_data_type='annotations', no_check=True)
        else:
            annotations_step = self.project.find_previous_step_of_type(self._input_step(), 'annotations')
        return Partitions(make_partitions=not a.no_partitions, annotations_step=annotations_step)

    def run(self, step_data):
        from .raxml import create_raxml_data
        return create_raxml_data(step_data, self._input_step(), self._partitions(), self.args.run)

    def finish(self, step_obj):
        from .raxml import finish_raxml_data
        finish_raxml_data(step_obj)


class MrBayes(RAxML):
    _COMMAND = 'mr_bayes'
    _HELP = "Run MrBayes on alignment(s)"
    _STEP_BASE_NAME = 'MrBayes'

    @staticmethod
    def set_arguments(parser):
        RAxML.set_arguments(parser)
        parser.add_argument('--ngen', default=1000000, type=int, help='ngen parameter')
        parser.add_argument('--samplefreq', default=1000, type=int, help='samplefreq parameter')
        parser.add_argument('--burnin', type=int, help='burnin parameter')
        parser.add_argument('--burninfrac', type=float, help='burninfrac parameter')
        parser.add_argument('--nchains', default=4, type=int, help='nchains parameter')

    def run(self, step_data):
        from .mr_bayes import create_mr_bayes_data
        a = self.args
        return create_mr_bayes_data(step_data, self._input_step(), a, self._partitions(), a.run)

    def finish(self, step_obj):
        from .mr_bayes import finish_mr_bayes_data
        finish_mr_bayes_data(step_obj)
