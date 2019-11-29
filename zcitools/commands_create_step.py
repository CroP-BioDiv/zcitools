import os.path
from .steps import read_step

# Check command_classes.py for description


class _CreateStepCommand:
    _COMMAND_TYPE = 'new_step'
    _COMMAND = None
    _STEP_BASE_NAME = None

    def __init__(self, args):
        self.args = args

    def prev_steps(self):
        return [os.path.normpath(p) for p in self._prev_steps()]

    def _prev_steps(self):
        return []

    def cache_identifier(self):
        return None

    def run(self, step_data):
        raise NotImplementedError(f'Method {self.__class__.__name__}.run(step_data) is not implemented!')

    def step_base_name(self):
        return self._STEP_BASE_NAME or self._COMMAND

    def finish(self, step_obj):
        pass


class TableStep(_CreateStepCommand):
    _COMMAND = 'table'
    _HELP = """
Creates table step. Mandatory argument is input data file.
Additional arguments specify how to interpret input data.
"""

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('filename', help='Input data')
        parser.add_argument('-f', '--format', help='Input data format (if needed). Values: excel, csv, txt')
        parser.add_argument('-c', '--columns', help='Columns. Format name1,type1:name2,type2:...')

    def run(self, step_data):
        from .processing.input_file import create_table_step
        args = self.args
        return create_table_step(step_data, args.filename, data_format=args.format, columns=args.columns)


class NCBIStep(_CreateStepCommand):
    _COMMAND = 'ncbi'
    _HELP = "Creates sequences step. Mandatory argument is a table step."
    _STEP_BASE_NAME = 'NCBI'

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Input table step')
        parser.add_argument(
            '-f', '--force-download', action='store_true', help='Download even if step already contains data.')

    def _prev_steps(self):
        return [self.args.step]

    def cache_identifier(self):
        return dict(static=True, data_identifier=['NCBI'])

    def run(self, step_data):
        from .processing.sequence.ncbi import download_ncbi
        step = read_step(self.args.step, check_data_type='table')
        return download_ncbi(step_data, step, force_download=self.args.force_download)


# Annotations
class GeSeqStep(_CreateStepCommand):
    _COMMAND = 'ge_seq'
    _HELP = "Annotates chloroplast sequences with GeSeq"
    _STEP_BASE_NAME = 'GeSeq'

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Input sequences step')

    def _prev_steps(self):
        return [self.args.step]

    def run(self, step_data):
        from .processing.annotation.ge_seq import create_ge_seq_data
        step = read_step(self.args.step, check_data_type='sequences')
        return create_ge_seq_data(step_data, step)

    def finish(self, step_obj):
        from .processing.annotation.ge_seq import finish_ge_seq_data
        finish_ge_seq_data(step_obj)


class GeSeqStep(GeSeqStep):
    _COMMAND = 'cpgavas'
    _HELP = "Annotates chloroplast sequences with CPGAVAS"
    _STEP_BASE_NAME = 'CPGAVAS'

    def run(self, step_data):
        from .processing.annotation.cpgavas import create_cpgavas_data
        step = read_step(self.args.step, check_data_type='sequences')
        return create_cpgavas_data(step_data, step)

    def finish(self, step_obj):
        from .processing.annotation.ge_seq import finish_cpgavas_data
        finish_cpgavas_data(step_obj)
