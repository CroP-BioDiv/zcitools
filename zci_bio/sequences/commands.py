from step_project.base_commands import ProjectCommand, CreateStepCommand


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


class SequenceReadsStep(ProjectCommand):
    _COMMAND = 'sequence_reads'
    _HELP = "Browse given directory and find default sequences reads data."

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('-d', '--directory', default='.', help='Directory to browse. Default working directory.')
        parser.add_argument(
            '-o', '--output-filename', default=False, nargs='?', help='Save data. Default <dir>/sequence_reads.yml')
        #
        parser.add_argument('-i', '--instrument', help='Stores instrument to reads')
        parser.add_argument('-l', '--read-length', type=int, help='Stores read length to reads')
        parser.add_argument('-g', '--gap-length', type=int, help='Stores gap length to reads')

    def run(self):
        import os.path
        from ..utils.sequence_reads import SequenceReads
        args = self.args
        seq_reads = SequenceReads.from_directory(
            args.directory, instrument=args.instrument, read_length=args.read_length, gap_length=args.gap_length)
        print('Sequence reads data:')
        seq_reads.print_data()
        if args.output_filename is not False:  # False means that -o is not set. None means -o without value
            seq_reads.write_data(args.output_filename or os.path.join(args.directory, 'sequence_reads.yml'))
