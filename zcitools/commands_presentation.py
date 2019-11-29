# Check command_classes.py for description


class _PresentationCommand:
    _COMMAND = None
    _COMMAND_TYPE = 'presentation'
    _STEP_DATA_TYPE = None         # Type of step commnad works on
    # Only one of these should be set
    _CALCULATION_DIRECTORY = None  # Name of subdirectory result
    _CALCULATION_FILENAME = None   # Name of file result

    def __init__(self, args):
        self.args = args

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Step name')
        parser.add_argument(
            '-f', '--force', action='store_true', help='Force recalculation (removes existing directory data)')

    def run(self, step):
        raise NotImplementedError(f'Method {self.__class__.__name__}.run(step) is not implemented!')


class RSCU(_PresentationCommand):
    _COMMAND = 'rscu'
    _STEP_DATA_TYPE = 'annotations'
    _CALCULATION_DIRECTORY = 'RSCU'
    _HELP = "Calculates sequence RSCUs"

    def run(self, step):
        from .processing.sequence.cai import calcualte_rscu
        calcualte_rscu(step, self._CALCULATION_DIRECTORY)


class OGDRAW(_PresentationCommand):
    _COMMAND = 'ogdraw'
    _STEP_DATA_TYPE = 'annotations'
    _CALCULATION_DIRECTORY = 'OGDraw'
    _HELP = "Create OGDraw images of annotations"

    def run(self, step):
        from .processing.annotation.ogdraw import calculate_ogdraw
        calculate_ogdraw(step, self._CALCULATION_DIRECTORY)
