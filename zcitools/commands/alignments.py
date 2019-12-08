from .base import _CreateStepFromStepCommand
from ..utils.exceptions import ZCItoolsValueError

# Check base.py for description


class ClustalO(_CreateStepFromStepCommand):
    _COMMAND = 'clustal'
    _HELP = "Align sequences with Clustal Omega"
    _STEP_BASE_NAME = 'Clustal'
    _INPUT_STEP_DATA_TYPE = 'annotations'
    _ALIGNMENTS = (('w', 'whole'),
                   ('gs', 'genes single'), ('gc', 'genes concatenate'),
                   ('cs', 'CDS single'), ('cc', 'CDS concatenated'))
    # Note: no global caching

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Input sequences step')
        parser.add_argument(
            'alignments', nargs='+', help=f"To align: {', '.join(f'{c} ({d})' for c, d in ClustalO._ALIGNMENTS)}")
        parser.add_argument('-r', '--run', action='store_true', help='Run Clustal Omega locale')

    def run(self, step_data):
        from ..processing.alignment.clustal_omega import create_clustal_data
        sup_aligns = [a for a, _ in self._ALIGNMENTS]
        align_params = [a.lower() for a in self.args.alignments]
        align_params = [a for a in align_params if a in sup_aligns]
        if not align_params:
            raise ZCItoolsValueError('No valid alignments set ({self.args.alignments}).')
        return create_clustal_data(step_data, self._input_step(), self.get_cache_object(), align_params, self.args.run)

    def finish(self, step_obj):
        from ..processing.alignment.clustal_omega import finish_clustal_data
        finish_clustal_data(step_obj, self.get_cache_object())
