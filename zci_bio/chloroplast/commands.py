from step_project.base_commands import ProjectCommand, CreateStepFromStepCommand, CreateStepsFromStepCommand


class ChloroplastAnalyse(CreateStepFromStepCommand):
    _COMMAND = 'chloroplast_analyse'
    _HELP = "Analyse (and normalize) chloroplast genomes"
    _COMMAND_GROUP = 'Chloroplast'
    _STEP_BASE_NAME = 'ChloroAnalyse'
    _INPUT_STEP_DATA_TYPE = 'table'

    # @classmethod
    # def set_arguments(cls, parser):
    #     CreateStepFromStepCommand.set_arguments(parser)
    #     parser.add_argument('-r', '--run', action='store_true', help='Post data on mVISTA web page')

    def run(self, step_data):
        from .analyse import analyse_genomes
        return analyse_genomes_start(step_data, self._input_step())

    def finish(self, step_obj):
        from .analyse_genomes import analyse_genomes_proceed
        analyse_genomes_proceed(step_obj)


class ChloroplastAnnotation(ProjectCommand):
    _COMMAND = 'chloroplast_annotation'
    _HELP = "Displays simplified chloroplast annotation"
    _COMMAND_GROUP = 'Chloroplast'

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Input step')
        parser.add_argument('-n', '--num-genes', default=1, type=int, help='Number of genes to display on region ends')
        parser.add_argument('-t', '--feature-type', default='gene', help='Feature type: gene or CDS')
        parser.add_argument('-f', '--features', help='List of features to take. Format feature1;feature2;...')
        parser.add_argument('-s', '--sequences', help='List of sequences to print. Format seq1;seq2;...')

    def run(self):
        from .annotation import chloroplast_annotation
        a = self.args
        chloroplast_annotation(
            self.project.read_step(a.step, check_data_type='annotations'),
            num_genes=a.num_genes,
            feature_type=a.feature_type,
            features=a.features.split(';') if a.features else None,
            sequences=a.sequences.split(';') if a.sequences else None)


class ChloroplastOrientate(CreateStepsFromStepCommand):
    _COMMAND = 'chloroplast_orientate'
    _HELP = "Orientate chloroplast sequences in standard way"
    _COMMAND_GROUP = 'Chloroplast'
    _STEP_BASE_NAME = 'ChloroOrient'
    _INPUT_STEP_DATA_TYPE = 'annotations'

    @classmethod
    def set_arguments(cls, parser):
        CreateStepsFromStepCommand.set_arguments(parser)
        parser.add_argument('-a', '--annotation-command', help='Annotation command to use on repaired seqeunces')

    def common_db_identifier(self):
        from .orientate import CHLOROPLAST_ORIENTED_DB_NAME
        return self.sequence_db_identifier(CHLOROPLAST_ORIENTED_DB_NAME, 'sequences')

    def run(self, step_data):
        from .orientate import orientate_chloroplast
        return orientate_chloroplast(self, self.args, step_data, self._input_step(), self.get_common_db_object())


class ChloroplastIRsFind(CreateStepFromStepCommand):
    _COMMAND = 'chloroplast_irs_find'
    _HELP = "Find chloroplast IRs and other repeats"
    _COMMAND_GROUP = 'Chloroplast'
    _STEP_BASE_NAME = 'ChloroIRs'
    _INPUT_STEP_DATA_TYPE = ('sequences', 'annotations')
    # Note: mummer is fast enough so it is not needed to run it on server.
    #       Implementation is done in a way that process can be split in parts (init, run, finish)

    @classmethod
    def set_arguments(cls, parser):
        CreateStepFromStepCommand.set_arguments(parser)
        parser.add_argument(
            '-m', '--force-mummer-parse', action='store_true',
            help='Force parsing of Mummer output even if calculation is not needed!')

    def common_db_identifier(self):
        db = self._input_step().common_db_identifier()[1]
        return self.sequence_db_identifier(db, 'IRs_mummer')

    def run(self, step_data):
        from .inverted_repeats import create_irs_data
        return create_irs_data(step_data, self._input_step(), self.args, self.get_common_db_object())


class ChloroplastIRsShow(ProjectCommand):
    # For testing purpose
    _COMMAND = 'chloroplast_irs_show'
    _HELP = "Show found chloroplast IRs"
    _COMMAND_GROUP = 'Chloroplast'

    @classmethod
    def set_arguments(cls, parser):
        parser.add_argument('step', help='Step name')

    def run(self):
        from .inverted_repeats import show_irs_data
        return show_irs_data(self.project.read_step(self.args.step, check_data_type='annotations'))
