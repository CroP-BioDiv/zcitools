from step_project.base_commands import ProjectCommand, CreateStepFromStepCommand


# Step 01: check chloroplast IRs
class ChloroplastAnalyse(ProjectCommand):
    _COMMAND = 'chloroplast_analyse'
    _HELP = "Analyse (and normalize) chloroplast genomes"
    _COMMAND_GROUP = 'Chloroplast'

    @staticmethod
    def set_arguments(parser):
        ProjectCommand.set_arguments(parser)
        parser.add_argument('step', help='Input sequences step')
        parser.add_argument('-o', '--output-file', default='chloroplast_analyse.xlsx', help='Output Excel file')

    def common_db_identifier(self):
        return self.sequence_db_identifier('base')

    def run(self):
        from .analyse import analyse_genomes_start
        a = self.args
        input_step = self.project.read_step(a.step, check_data_type='table')
        return analyse_genomes_start(input_step, a.output_file, self.get_common_db_object())


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


# Step 02: check orientation of chloroplast by one reference
class ChloroplastCheckOrientation(CreateStepFromStepCommand):
    _COMMAND = 'chloroplast_check_orientation'
    _HELP = "Check orientation of chloroplast sequences by reference genome"
    _COMMAND_GROUP = 'Chloroplast'
    _STEP_BASE_NAME = 'ChloroCheckOrient'
    _INPUT_STEP_DATA_TYPE = 'annotations'

    @classmethod
    def set_arguments(cls, parser):
        CreateStepFromStepCommand.set_arguments(parser)
        parser.add_argument('-r', '--referent-genome', default='TAN_', help='Referent genome. Whole or start')
        parser.add_argument('-l', '--length-to-check', type=int, default=1000, help='Length on sequence start to check')
        parser.add_argument('-p', '--force-parse', action='store_true', help='Force parsing of alignment files')
        parser.add_argument('-o', '--output-file', default='chloroplast_orientate.xlsx', help='Output Excel file')

    # def common_db_identifier(self):
    #     from .orientate import CHLOROPLAST_ORIENTED_DB_NAME
    #     return self.sequence_db_identifier(CHLOROPLAST_ORIENTED_DB_NAME, 'sequences')

    def run(self, step_data):
        from .orientate import orientate_chloroplast_start
        return orientate_chloroplast_start(step_data, self._input_step(), self.args)

    def finish(self, step_obj):
        from .orientate import orientate_chloroplast_finish
        orientate_chloroplast_finish(step_obj)  # , self.get_step_db_object(step_obj))


# # Step 03: orientate chloroplast sequences by Fast-Plast method
# class ChloroplastCheckOrientation(CreateStepFromStepCommand):
#     _COMMAND = 'chloroplast_orientate'
#     _HELP = "Orientate chloroplast sequences by Fast-Plat method"
#     _COMMAND_GROUP = 'Chloroplast'
#     _STEP_BASE_NAME = 'ChloroCheckOrient'
#     _INPUT_STEP_DATA_TYPE = 'annotations'

#     @classmethod
#     def set_arguments(cls, parser):
#         CreateStepFromStepCommand.set_arguments(parser)

#     # def common_db_identifier(self):
#     #     from .orientate import CHLOROPLAST_ORIENTED_DB_NAME
#     #     return self.sequence_db_identifier(CHLOROPLAST_ORIENTED_DB_NAME, 'sequences')

#     def run(self, step_data):
#         from .orientate import orientate_chloroplast_start
#         return orientate_chloroplast_start(step_data, self._input_step(), self.args)

#     def finish(self, step_obj):
#         from .orientate import orientate_chloroplast_finish
#         orientate_chloroplast_finish(step_obj)  # , self.get_step_db_object(step_obj))


# Info commands
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
