from step_project.base_commands import ProjectCommand, NonProjectCommand, CreateStepCommand, CreateStepFromStepCommand
from step_project.common.common_db_commands import CommonDBCommand, CommonDB
from common_utils.exceptions import ZCItoolsValueError


class CheckAnnotations(NonProjectCommand):
    _COMMAND = 'check_annotations'
    _HELP = "Checks annotations for problems"
    _COMMAND_GROUP = 'Bio'

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('files', nargs='+', help='Input files or directories')
        parser.add_argument('-t', '--filter-type', nargs='*', help='Filter by feature type.')
        parser.add_argument('-o', '--output-filename', help='Output filename')

    def run(self):
        from ..utils.features import check_annotations
        a = self.args
        check_annotations(a.files, a.filter_type, a.output_filename)


class FeatureProperties(ProjectCommand):
    _COMMAND = 'feature_properties'
    _HELP = "Displays feature properties"
    _COMMAND_GROUP = 'Bio'

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Input annotation step')
        parser.add_argument('feature', help='Feature to display')

    def run(self):
        from ..utils.helpers import feature_location_desc
        step = self.project.read_step(self.args.step, check_data_type='annotations', update_mode=False, no_check=True)
        name = self.args.feature
        ret = []
        for seq_ident, seq in step._iterate_records():
            features = [f for f in seq.features if f.type == 'gene' and f.qualifiers['gene'][0] == name]
            if features:
                ret.append((seq_ident, len(seq), [(f.type, feature_location_desc(f.location)) for f in features]))
        for x in ret:
            print(*x)


class AnnotationsMerge(CreateStepCommand):
    _COMMAND = 'annotation_merge'
    _HELP = "Merge more annotations steps into one"
    _COMMAND_GROUP = 'Bio'

    @staticmethod
    def set_arguments(parser):
        CreateStepCommand.set_arguments(parser)
        parser.add_argument('steps', nargs='+', help='Input steps')

    def run(self, step_data):
        from .steps import AnnotationsStep
        from common_utils.file_utils import copy_file

        step = AnnotationsStep(self.project, step_data, remove_data=True)
        for s in self.args.steps:
            a_step = self.project.read_step(s, check_data_type='annotations')
            for seq_ident in a_step.all_sequences():
                step._sequences.add(seq_ident)
                copy_file(a_step.get_sequence_filename(seq_ident), step.step_file(seq_ident + '.gb'))
        step.save()
        return step


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
                 for f in seq_rec.features if location and f.type == type_ and 'gene' in f.qualifiers))


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
class GeSeq(CreateStepFromStepCommand):
    _COMMAND = 'ge_seq'
    _HELP = "Annotates chloroplast sequences with GeSeq"
    _STEP_BASE_NAME = 'GeSeq'
    _INPUT_STEP_DATA_TYPE = 'sequences'
    _COMMON_DB_IDENT = ('GeSeq',)

    @staticmethod
    def set_arguments(parser):
        CreateStepFromStepCommand.set_arguments(parser)
        parser.add_argument('-s', '--num-sequences-in-file', default=50, type=int,
                            help='Number of sequences in upload files')

    def run(self, step_data):
        from .ge_seq import create_ge_seq_data
        return create_ge_seq_data(step_data, self._input_step(), self.get_common_db_object(),
                                  self.args.num_sequences_in_file)

    def finish(self, step_obj):
        from .ge_seq import finish_ge_seq_data
        finish_ge_seq_data(step_obj, self.get_step_db_object(step_obj))


class GeSeqCDBClean(CommonDBCommand):
    _COMMAND = 'ge_seq_cache_clean'
    _HELP = "Remove outdated GeSeq records from CommonDB"

    def run(self):
        import io
        from ..utils.import_methods import import_bio_seq_io
        from ..utils.helpers import get_bio_io_type
        SeqIO = import_bio_seq_io()
        self._set_db_attrs(('GeSeq',))
        seq_db = CommonDB.get_zci_db(('sequences',))
        for rel_dir, rec in self._iterate_records():
            record_ident = self._db_ident(rel_dir, rec)
            if gs_ret := self._common_db.get_one_file(record_ident):
                if seq_ret := seq_db.get_one_file(record_ident):
                    _, g_filename, g_data = gs_ret
                    _, s_filename, s_data = seq_ret
                    g_seq = SeqIO.read(io.StringIO(g_data.decode('utf-8')), get_bio_io_type(g_filename))
                    s_seq = SeqIO.read(io.StringIO(s_data.decode('utf-8')), get_bio_io_type(s_filename))
                    if g_seq.seq != s_seq.seq:
                        print('ToDo: DELETE', f)


class CPGAVAS(CreateStepFromStepCommand):
    _COMMAND = 'cpgavas'
    _HELP = "Annotates chloroplast sequences with CPGAVAS"
    _STEP_BASE_NAME = 'CPGAVAS'
    _INPUT_STEP_DATA_TYPE = 'sequences'
    _COMMON_DB_IDENT = ('CPGAVAS',)

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
        CreateStepFromStepCommand.set_arguments(parser)
        parser.add_argument('-f', '--image_format', default='ps', help='One of: svg, pdf, ps, png, jpg, tif, gif')
        parser.add_argument('-s', '--sequences', help="Filter only sequences, separate seq_idents by ';'.")

    def common_db_identifier(self):
        step = self._input_step()
        step_command = step.get_command()  # Depends on annotation process
        return ('OGDraw', step_command, self.args.image_format)

    def run(self, step_data):
        from .ogdraw import create_ogdraw
        img_f = self.args.image_format.lower()
        if img_f not in self._IMAGE_FORMATS:
            raise ZCItoolsValueError(f'Given format {img_f} is not supported!')
        return create_ogdraw(
            step_data, img_f, self._input_step(), self.get_common_db_object(), sequences=self.args.sequences)

    def finish(self, step_obj):
        from .ogdraw import finish_ogdraw
        finish_ogdraw(step_obj, self.get_step_db_object(step_obj))


class OGDRAW_ExtractJPG(NonProjectCommand):
    _COMMAND = 'ogdraw_jpg'
    _HELP = "Extract OGDraw jpg image into CommonDB, from GeSeq annotation result file"

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Input sequences step')
        # ToDo: remove
        parser.add_argument('-S', '--sequence-db', help=f'Sequence database to use. For old projects!')

    def common_db_identifier(self):
        step = self.project.read_step(self.args.step, check_data_type='annotations', no_check=True)
        step_command = step.get_command()  # Depends on annotation process
        return ('OGDraw', step_command, 'jpg')

    def run(self):
        from .ogdraw import extract_jpg_to_common_db
        step = self.project.read_step(self.args.step, check_data_type='annotations', no_check=True)
        return extract_jpg_to_common_db(step, self.get_common_db_object())


# class RSCU(_PresentationCommand):
#     _COMMAND = 'rscu'
#     _STEP_DATA_TYPE = 'annotations'
#     _CALCULATION_DIRECTORY = 'RSCU'
#     _HELP = "Calculates sequence RSCUs"

#     def run(self, step):
#         from .processing.sequence.cai import calcualte_rscu
#         calcualte_rscu(step, self._CALCULATION_DIRECTORY)
