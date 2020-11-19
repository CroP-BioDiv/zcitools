from step_project.base_commands import CreateStepFromStepCommand
from common_utils.exceptions import ZCItoolsValueError


class AlignGenomes(CreateStepFromStepCommand):
    _COMMAND = 'align_genomes'
    _HELP = "Align sequences (using one of alignment programs)"
    _INPUT_STEP_DATA_TYPE = ('sequences', 'annotations')
    _ALIGNMENTS = (('w', 'whole'),
                   ('gs', 'genes single'), ('gc', 'genes concatenate'),
                   ('cs', 'CDS single'), ('cc', 'CDS concatenated'))
    _STEP_NAMES = dict(clustal_omega='Clustal', mafft='MAFFT', muscle='MUSCLE')

    @staticmethod
    def set_arguments(parser):
        from .common_methods import _align_programs
        CreateStepFromStepCommand.set_arguments(parser)
        parser.add_argument('-p', '--alignment_program', default='mafft',
                            choices=sorted(_align_programs.keys()),
                            help=f"Alignment program to use.")
        parser.add_argument(
            'alignments', nargs='+', help=f"To align: {', '.join(f'{c} ({d})' for c, d in AlignGenomes._ALIGNMENTS)}")
        parser.add_argument(
            '-w', '--whole-partition', choices=['gene', 'CDS'],
            help='For whole alignment which set partition type to use. Values: gene, CDS.')
        parser.add_argument('-r', '--run', action='store_true', help='Run Clustal Omega locale')

    def step_base_name(self):
        return self._format_step_name(self._STEP_NAMES[self.args.alignment_program])

    def run(self, step_data):
        from .common_methods import create_alignment_data
        align_params = [a.lower() for a in self.args.alignments]
        sup = set(a for a, _ in self._ALIGNMENTS)
        if not_sup_params := [a for a in align_params if a not in sup]:
            raise ZCItoolsValueError(f'No valid alignments set ({", ".join(sorted(not_sup_params))}).')
        return create_alignment_data(step_data, self._input_step(), align_params, self.args.whole_partition,
                                     self.args.run, self.args.alignment_program)

    def finish(self, step_obj):
        from .common_methods import finish_alignment_data
        # Check are needed files in zip, not something strange
        files = set(d['filename'].replace('sequences.fa', 'alignment.phy')
                    for d in read_yaml(step_obj.step_file('finish.yml')))
        finish_alignment_data(step_obj, files)


# -------------------------------------------------------------------
# Not in use!
class mVISTA(CreateStepFromStepCommand):
    _COMMAND = 'mvista'
    _HELP = "Align sequences with mVISTA programs"
    _STEP_BASE_NAME = 'mVISTA'
    _INPUT_STEP_DATA_TYPE = 'annotations'
    # Note: no global caching

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Input sequences step')
        parser.add_argument('-r', '--run', action='store_true', help='Post data on mVISTA web page')
        parser.add_argument('-e', '--email', help='email address (needed for posting, default in settings)')

    def run(self, step_data):
        from .mvista import create_mvista_data
        from common_utils.file_utils import get_settings
        email = self.args.email or get_settings()['email']
        if self.args.run and not email:
            raise ZCItoolsValueError('Email address is needed to post mVISTA data!')
        return create_mvista_data(step_data, self._input_step(), self.args.run, email)

    def finish(self, step_obj):
        from .mvista import finish_mvista_data
        finish_mvista_data(step_obj)
