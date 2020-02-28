from step_project.base_commands import ProjectCommand, CreateStepsFromStepCommand


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

    def db_identifier(self):
        # Same as in input step
        return self._input_step().db_identifier()

    def sequence_db(self):
        from .orientate import CHLOROPLAST_ORIENTED_DB_NAME
        print('aaaaaa', CHLOROPLAST_ORIENTED_DB_NAME)
        return CHLOROPLAST_ORIENTED_DB_NAME

    def run(self, step_data):
        from .orientate import orientate_chloroplast
        return orientate_chloroplast(self, self.args, step_data, self._input_step(), self.get_common_db_object())
