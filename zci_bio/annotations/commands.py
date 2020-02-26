from step_project.base_commands import ProjectCommand, NonProjectCommand, \
    CreateStepFromStepCommand, CreateStepFromStepCommand_CommonDB
from common_utils.exceptions import ZCItoolsValueError


class AnnotationDiff(ProjectCommand):
    _COMMAND = 'annotation_diff'
    _HELP = "Displays sequence distance marix based on annotation"
    _COMMAND_GROUP = 'Bio'

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Input step')

    def run(self):
        step = self.project.read_step(self.args.step, check_data_type='alignment')
        return step.diff_matrix()


# ToDo: ove prebaciti da rade s dummy annotation step objektom. Treba imati proxy feature objekt!
class AnnotationExtract(NonProjectCommand):
    _COMMAND = 'annotation_extract'
    _HELP = "Extract annotation parts into separate file(s)"
    _COMMAND_GROUP = 'Bio'

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('input_files', nargs='+', help='Input filenames')
        parser.add_argument('-i', '--input-format', help='Force input format')
        # parser.add_argument('-o', '--output-format', default='fasta', help='Output format')
        parser.add_argument('-d', '--output-directory', default='.', help='Output directory')
        parser.add_argument('-t', '--filter-type', default='CDS', help='Filter feature type')

    def run(self):
        import os.path
        from ..utils.import_methods import import_bio_seq_io
        from ..utils.helpers import get_bio_io_type, feature_qualifiers_to_desc
        from common_utils.file_utils import ensure_directory, write_fasta, basename_no_ext

        args = self.args
        od = args.output_directory
        SeqIO = import_bio_seq_io()
        ensure_directory(od)
        type_ = args.filter_type

        for i_filename in args.input_files:
            # Note: One sequence in one file!
            seq_rec = SeqIO.read(i_filename, get_bio_io_type(i_filename, args.input_format))
            # ToDo: filtrirati po necemu?
            # ToDo: sortirati po necemu?
            write_fasta(
                os.path.join(od, f"extract_{basename_no_ext(i_filename)}.fasta"),
                ((feature_qualifiers_to_desc(f), str(f.extract(seq_rec).seq))
                 for f in seq_rec.features if f.type == type_ and 'gene' in f.qualifiers))


class AnnotationGroup(AnnotationExtract):  # Works only with given file
    _COMMAND = 'annotation_group'
    _HELP = "Extract annotation parts and group them in files"

    def run(self):
        import os.path
        from collections import defaultdict
        from ..utils.import_methods import import_bio_seq_io
        from ..utils.helpers import get_bio_io_type
        from common_utils.file_utils import ensure_directory, write_fasta

        args = self.args
        od = args.output_directory
        SeqIO = import_bio_seq_io()
        ensure_directory(od)
        genes = defaultdict(dict)  # gene -> dict(species -> data)

        for i_filename in args.input_files:
            for seq in SeqIO.parse(i_filename, get_bio_io_type(i_filename, args.input_format)):
                name = seq.id
                split_on = name.index('_')
                gene = name[:split_on]
                species = name[(split_on + 1):]
                genes[gene][species] = seq.seq

        for gene, data in genes.items():
            write_fasta(os.path.join(od, f'{gene}.fasta'), sorted(data.items()))


# ---------------------------------------------------------
class GeSeqStep(CreateStepFromStepCommand_CommonDB):
    _COMMAND = 'ge_seq'
    _HELP = "Annotates chloroplast sequences with GeSeq"
    _STEP_BASE_NAME = 'GeSeq'
    _INPUT_STEP_DATA_TYPE = 'sequences'
    _COMMON_DB_IDENT = 'GeSeq'

    def run(self, step_data):
        from .ge_seq import create_ge_seq_data
        return create_ge_seq_data(step_data, self._input_step(), self.get_common_db_object())

    def finish(self, step_obj):
        from .ge_seq import finish_ge_seq_data
        finish_ge_seq_data(step_obj, self.get_common_db_object())


class CPGAVAS(CreateStepFromStepCommand_CommonDB):
    _COMMAND = 'cpgavas'
    _HELP = "Annotates chloroplast sequences with CPGAVAS"
    _STEP_BASE_NAME = 'CPGAVAS'
    _INPUT_STEP_DATA_TYPE = 'sequences'
    _COMMON_DB_IDENT = 'CPGAVAS'

    def run(self, step_data):
        from .cpgavas import create_cpgavas_data
        return create_cpgavas_data(step_data, self._input_step())

    def finish(self, step_obj):
        from .ge_seq import finish_cpgavas_data
        finish_cpgavas_data(step_obj)


# Presentations
class OGDRAW(CreateStepFromStepCommand):
    _COMMAND = 'ogdraw'
    _HELP = "Create OGDraw images of annotations"
    _STEP_BASE_NAME = 'OGDraw'
    _PRESENTATION = True
    _INPUT_STEP_DATA_TYPE = 'annotations'
    _IMAGE_FORMATS = ('svg', 'pdf', 'ps', 'png', 'jpg', 'tif', 'gif')

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Input sequences step')
        parser.add_argument('-f', '--image_format', default='ps', help='One of: svg, pdf, ps, png, jpg, tif, gif')
        parser.add_argument('-s', '--sequences', help="Filter only sequences, separate seq_idents by ';'.")
        dbs = CreateStepFromStepCommand.get_sequence_dbs()
        parser.add_argument('-d', '--database', default='base', help=f'Database to use: {", ".join(dbs)}')

    def db_identifier(self):
        step_command = self._input_step().get_command()  # Depends on annotation process
        return dict(static=True, data_identifier=self.get_sequence_db_ident(
            'OGDraw', step_command, self.args.image_format))

    def run(self, step_data):
        from .ogdraw import create_ogdraw
        img_f = self.args.image_format.lower()
        if img_f not in self._IMAGE_FORMATS:
            raise ZCItoolsValueError(f'Given format {img_f} is not supported!')
        return create_ogdraw(
            step_data, img_f, self._input_step(), self.get_common_db_object(), sequences=self.args.sequences)

    def finish(self, step_obj):
        from .ogdraw import finish_ogdraw
        finish_ogdraw(step_obj, self.get_common_db_object())


# class RSCU(_PresentationCommand):
#     _COMMAND = 'rscu'
#     _STEP_DATA_TYPE = 'annotations'
#     _CALCULATION_DIRECTORY = 'RSCU'
#     _HELP = "Calculates sequence RSCUs"

#     def run(self, step):
#         from .processing.sequence.cai import calcualte_rscu
#         calcualte_rscu(step, self._CALCULATION_DIRECTORY)
