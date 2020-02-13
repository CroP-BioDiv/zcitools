from step_project.base_commands import NonProjectCommand, CreateStepCommand


class FetchSequencesStep(CreateStepCommand):
    _COMMAND = 'fetch_seqs'
    _HELP = "Creates sequences step. Mandatory argument is a table step."
    _STEP_BASE_NAME = 'seqs'

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('step', help='Input table step')

    def _prev_steps(self):
        return [self.args.step]

    def cache_identifier(self):
        return dict(static=True, data_identifier=['sequences'])

    def run(self, step_data):
        from .fetch import fetch_sequences
        step = self.project.read_step(self.args.step, check_data_type='table')
        return fetch_sequences(step_data, step, self.get_cache_object())


class SequenceReadsStep(NonProjectCommand):
    _COMMAND = 'sequence_reads'
    _HELP = "Browse given directory and find default sequences reads data."
    _COMMAND_GROUP = 'Bio'

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('-d', '--directory', default='.', help='Directory to browse. Default working directory.')
        parser.add_argument(
            '-o', '--output-filename', default=False, nargs='?', help='Save data. Default <dir>/sequence_reads.yml')
        #
        parser.add_argument('-p', '--platform', default='Illumina', help='Sequencing platform')
        parser.add_argument('-l', '--read-length', type=int, help='Read length')
        parser.add_argument('-i', '--insert-length', type=int, help='Insert length')

    def run(self):
        import os.path
        from ..utils.sequence_reads import SequenceReads
        args = self.args
        seq_reads = SequenceReads.from_directory(
            args.directory, platform=args.platform, read_length=args.read_length, insert_length=args.insert_length)
        print('Sequence reads data:')
        seq_reads.print_data()
        if args.output_filename is not False:  # False means that -o is not set. None means -o without value
            seq_reads.write_data(args.output_filename or os.path.join(args.directory, 'sequence_reads.yml'))


class ConvertSequence(NonProjectCommand):
    _COMMAND = 'convert_sequence'
    _HELP = "Converts sequence(s) data to given format"

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('input_filename', help='Input sequence file')
        parser.add_argument('output_filename', help='Output sequence file')
        parser.add_argument(
            '-i', '--input-format', help='Input sequence format, if not given deduced by file extention.')
        parser.add_argument(
            '-o', '--output-format', help='Output sequence format, if not given deduced by file extention.')

    def run(self):
        from ..utils.helpers import convert_sequence_file
        a = self.args
        convert_sequence_file(
            a.input_filename, a.output_filename, input_format=a.input_format, output_format=a.output_format)
