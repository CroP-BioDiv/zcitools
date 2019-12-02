from .base import _CreateStepFromStepCommand

# Check base.py for description


class GeSeqStep(_CreateStepFromStepCommand):
    _COMMAND = 'ge_seq'
    _HELP = "Annotates chloroplast sequences with GeSeq"
    _STEP_BASE_NAME = 'GeSeq'
    _INPUT_STEP_DATA_TYPE = 'sequences'

    def cache_identifier(self):
        return dict(static=True, data_identifier=['NCBI'])

    def run(self, step_data):
        from ..processing.annotation.ge_seq import create_ge_seq_data
        return create_ge_seq_data(step_data, self._input_step())

    def finish(self, step_obj):
        from ..processing.annotation.ge_seq import finish_ge_seq_data
        finish_ge_seq_data(step_obj)


class CPGAVAS(_CreateStepFromStepCommand):
    _COMMAND = 'cpgavas'
    _HELP = "Annotates chloroplast sequences with CPGAVAS"
    _STEP_BASE_NAME = 'CPGAVAS'
    _INPUT_STEP_DATA_TYPE = 'sequences'

    def cache_identifier(self):
        return dict(static=True, data_identifier=['CPGAVAS'])

    def run(self, step_data):
        from ..processing.annotation.cpgavas import create_cpgavas_data
        return create_cpgavas_data(step_data, self._input_step())

    def finish(self, step_obj):
        from ..processing.annotation.ge_seq import finish_cpgavas_data
        finish_cpgavas_data(step_obj)


# Presentations
class OGDRAW(_CreateStepFromStepCommand):
    _COMMAND = 'ogdraw'
    _HELP = "Create OGDraw images of annotations"
    _STEP_BASE_NAME = 'OGDraw'
    _PRESENTATION = True
    _INPUT_STEP_DATA_TYPE = 'annotations'

    def cache_identifier(self):
        return dict(static=True, data_identifier=['OGDraw'])

    def run(self, step_data):
        from ..processing.annotation.ogdraw import calculate_ogdraw
        return calculate_ogdraw(step_data, self._input_step())

    def finish(self, step_obj):
        from ..processing.annotation.ogdraw import finish_ogdraw
        finish_ogdraw(step_obj)


# class RSCU(_PresentationCommand):
#     _COMMAND = 'rscu'
#     _STEP_DATA_TYPE = 'annotations'
#     _CALCULATION_DIRECTORY = 'RSCU'
#     _HELP = "Calculates sequence RSCUs"

#     def run(self, step):
#         from .processing.sequence.cai import calcualte_rscu
#         calcualte_rscu(step, self._CALCULATION_DIRECTORY)
