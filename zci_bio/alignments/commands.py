from step_project.base_commands import CreateStepFromStepCommand
from common_utils.exceptions import ZCItoolsValueError
from common_utils.file_utils import get_settings


class ClustalO(CreateStepFromStepCommand):
    _COMMAND = 'clustal'
    _HELP = "Align sequences using Clustal Omega"
    _STEP_BASE_NAME = 'Clustal'
    _INPUT_STEP_DATA_TYPE = 'annotations'
    _ALIGNMENTS = (('w', 'whole'),
                   ('gs', 'genes single'), ('gc', 'genes concatenate'),
                   ('cs', 'CDS single'), ('cc', 'CDS concatenated'))
    _SUPPORTED_ALIGNS = set(a for a, _ in _ALIGNMENTS)
    # Note: no global caching

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Input sequences step')
        parser.add_argument(
            'alignments', nargs='+', help=f"To align: {', '.join(f'{c} ({d})' for c, d in ClustalO._ALIGNMENTS)}")
        parser.add_argument('-r', '--run', action='store_true', help='Run Clustal Omega locale')

    def run(self, step_data):
        from .clustal_omega import create_clustal_data
        align_params = [a.lower() for a in self.args.alignments]
        not_sup_params = [a for a in align_params if a not in self._SUPPORTED_ALIGNS]
        if not_sup_params:
            raise ZCItoolsValueError(f'No valid alignments set ({", ".join(sorted(not_sup_params))}).')
        return create_clustal_data(step_data, self._input_step(), align_params, self.args.run)

    def finish(self, step_obj):
        from .clustal_omega import finish_clustal_data
        finish_clustal_data(step_obj)


class MAFFT(ClustalO):
    _COMMAND = 'mafft'
    _HELP = "Align sequences using MAFFT"
    _STEP_BASE_NAME = 'MAFFT'

    def run(self, step_data):
        from .mafft import create_mafft_data
        align_params = [a.lower() for a in self.args.alignments]
        not_sup_params = [a for a in align_params if a not in self._SUPPORTED_ALIGNS]
        if not_sup_params:
            raise ZCItoolsValueError(f'No valid alignments set ({", ".join(sorted(not_sup_params))}).')
        return create_mafft_data(step_data, self._input_step(), align_params, self.args.run)

    def finish(self, step_obj):
        from .mafft import finish_mafft_data
        finish_mafft_data(step_obj)


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
        run = self.args.run
        email = self.args.email or get_settings()['email']
        if run and not email:
            raise ZCItoolsValueError('Email address is needed to post mVISTA data!')
        return create_mvista_data(step_data, self._input_step(), self.get_cache_object(), run, email)

    def finish(self, step_obj):
        from .mvista import finish_mvista_data
        finish_mvista_data(step_obj, self.get_cache_object())
