from zcitools.base.commands import CreateStepFromStepCommand
from common_utils.exceptions import ZCItoolsValueError

# Check base.py for description


class GeSeqStep(CreateStepFromStepCommand):
    _COMMAND = 'ge_seq'
    _HELP = "Annotates chloroplast sequences with GeSeq"
    _STEP_BASE_NAME = 'GeSeq'
    _INPUT_STEP_DATA_TYPE = 'sequences'

    def cache_identifier(self):
        return dict(static=True, data_identifier=['GeSeq'])

    def run(self, step_data):
        from ..processing.annotation.ge_seq import create_ge_seq_data
        return create_ge_seq_data(step_data, self._input_step(), self.get_cache_object())

    def finish(self, step_obj):
        from ..processing.annotation.ge_seq import finish_ge_seq_data
        finish_ge_seq_data(step_obj, self.get_cache_object())


class CPGAVAS(CreateStepFromStepCommand):
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

    def cache_identifier(self):
        step_commnad = self._input_step().get_command()  # Depends on annotation process
        return dict(static=True, data_identifier=['OGDraw', step_commnad, self.args.image_format])

    def run(self, step_data):
        from ..processing.annotation.ogdraw import create_ogdraw
        img_f = self.args.image_format.lower()
        if img_f not in self._IMAGE_FORMATS:
            raise ZCItoolsValueError(f'Given format {img_f} is not supported!')
        return create_ogdraw(step_data, img_f, self._input_step(), self.get_cache_object())

    def finish(self, step_obj):
        from ..processing.annotation.ogdraw import finish_ogdraw
        finish_ogdraw(step_obj, self.get_cache_object())


# class RSCU(_PresentationCommand):
#     _COMMAND = 'rscu'
#     _STEP_DATA_TYPE = 'annotations'
#     _CALCULATION_DIRECTORY = 'RSCU'
#     _HELP = "Calculates sequence RSCUs"

#     def run(self, step):
#         from .processing.sequence.cai import calcualte_rscu
#         calcualte_rscu(step, self._CALCULATION_DIRECTORY)
