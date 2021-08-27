from step_project.base_commands import ProjectCommand, CreateStepFromStepCommand
from common_utils.exceptions import ZCItoolsValueError


class RAxML(CreateStepFromStepCommand):
    _COMMAND = 'raxml'
    _HELP = "Run RAxML on alignment(s)"
    _COMMAND_GROUP = 'Phylogenetics'
    _STEP_BASE_NAME = 'RAxML'
    _INPUT_STEP_DATA_TYPE = ('alignment', 'alignments')

    @staticmethod
    def set_arguments(parser):
        CreateStepFromStepCommand.set_arguments(parser)
        parser.add_argument('-a', '--annotations-step', help='Annotations step to use for partitions')
        parser.add_argument('-p', '--no-partitions', action='store_true', help='Do not set partitions')
        parser.add_argument(
            '-r', '--run', const=-1, nargs='?', type=int,
            help=f'Run {RAxML._STEP_BASE_NAME} locale. Number of threads to use can be specified.')

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
        return create_raxml_data(step_data, self._input_step(), self._partitions(), self._run_threads())

    def finish(self, step_obj):
        from .raxml import finish_raxml_data
        finish_raxml_data(step_obj)


class MrBayes(RAxML):
    _COMMAND = 'mr_bayes'
    _HELP = "Run MrBayes on alignment(s)"
    _COMMAND_GROUP = 'Phylogenetics'
    _STEP_BASE_NAME = 'MrBayes'

    @staticmethod
    def set_arguments(parser):
        RAxML.set_arguments(parser)
        parser.add_argument('--no-mpi', action='store_true', help='Do not use MPI version')
        parser.add_argument('-R', '--num-runs', type=int, help='Number of runs to perform.')
        parser.add_argument('--ngen', default=1000000, type=int, help='Parameter: ngen')
        parser.add_argument('--samplefreq', default=1000, type=int, help='Parameter: samplefreq')
        parser.add_argument('--burnin', type=int, help='Parameter: burnin')
        parser.add_argument('--burninfrac', type=float, help='Parameter: burninfrac')
        parser.add_argument('--nchains', default=4, type=int, help='Parameter: nchains')

    def step_base_name(self):
        n = f'_R{nr}' if (nr := self.args.num_runs) and nr > 1 else ''
        return self._format_step_name(f"{self._STEP_BASE_NAME}_{'N' if self.args.no_partitions else 'P'}{n}")

    def run(self, step_data):
        from .mr_bayes import create_mr_bayes_data
        return create_mr_bayes_data(step_data, self._input_step(), self.args, self._partitions(), self._run_threads())

    def finish(self, step_obj):
        from .mr_bayes import finish_mr_bayes_data
        finish_mr_bayes_data(step_obj)


# Helper commands
class _RunOnMrBayes(ProjectCommand):
    _COMMAND_GROUP = 'Phylogenetics'

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Input MrBayes step')
        parser.add_argument('-e', '--exe-location', help='Location of executable. Dirname or filename')


class RunTracer(_RunOnMrBayes):
    _COMMAND = 'tracer'
    _HELP = "Run tracer on MrBayes step"

    def run(self):
        from .mr_bayes import run_tracer
        run_tracer(self.project.read_step(self.args.step, check_data_type=('mr_bayes', 'mr_bayes_s'), no_check=True),
                   self.args.exe_location)


class RunConsense(_RunOnMrBayes):
    _COMMAND = 'consense'
    _HELP = "Run consense on MrBayes step"

    def run(self):
        from .mr_bayes import run_consense
        run_consense(self.project.read_step(self.args.step, check_data_type=('mr_bayes', 'mr_bayes_s'), no_check=True),
                     self.args.exe_location)
