import numpy
from step_project.base_commands import NonProjectCommand, CreateStepFromStepCommand
from common_utils.exceptions import ZCItoolsValueError


class FetchSequencesStep(CreateStepFromStepCommand):
    _COMMAND = 'fetch_seqs'
    _HELP = "Creates sequences step. Mandatory argument is a table step."
    _STEP_BASE_NAME = 'seqs'
    _INPUT_STEP_DATA_TYPE = 'table'
    _COMMON_DB_IDENT = ('sequences',)

    @staticmethod
    def set_arguments(parser):
        CreateStepFromStepCommand.set_arguments(parser)
        parser.add_argument('-c', '--column-name', help='Table column to use for sequence identifiers')

    def run(self, step_data):
        from .fetch import fetch_sequences
        return fetch_sequences(
            step_data, self._input_step(), self.get_common_db_object(), column_name=self.args.column_name)


class UpdateCahcedSequencesStep(NonProjectCommand):
    _COMMAND = 'cdb_update_seqs'
    _COMMAND_GROUP = 'NCBI'
    _HELP = "Update cached sequences"

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('-s', '--start-seq-ident', help='Start with seq_ident')
        parser.add_argument('-e', '--end-seq-ident', help='End with seq_ident')

    def run(self):
        from .fetch import update_cached_sequences
        return update_cached_sequences(self.get_idents_common_db_object(tuple()),
                                       self.args.start_seq_ident, self.args.end_seq_ident)


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
        parser.add_argument('-s', '--as-single-reads', action='store_true',
                            help='Paired end data represent as not paired end.')

    def run(self):
        import os.path
        from ..utils.sequence_reads import SequenceReads
        args = self.args
        seq_reads = SequenceReads.from_directory(args.directory, args)
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


class ChangeSequence(NonProjectCommand):
    _COMMAND = 'change_sequence'
    _HELP = "Change sequence(s) data to given (in some way) and save result"

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('method', help='Change method: revert, translate')
        parser.add_argument('input_filename', help='Input sequence file')
        parser.add_argument('output_filename', help='Output sequence file')
        parser.add_argument(
            '-i', '--input-format', help='Input sequence format, if not given deduced by file extention.')
        parser.add_argument(
            '-o', '--output-format', help='Output sequence format, if not given deduced by file extention.')
        parser.add_argument('-p', '--position', type=int, help='Translate position')

    def run(self):
        from ..utils.helpers import change_sequence_data
        a = self.args
        method = a.method.lower()
        if method not in ('revert', 'translate'):
            raise ZCItoolsValueError(f"Not known method {a.method}!")
        change_sequence_data(
            method, a.input_filename, a.output_filename,
            input_format=a.input_format, output_format=a.output_format, position=a.position)


class SequencesEqual(NonProjectCommand):
    _COMMAND = 'sequences_equal'
    _HELP = """Check are given sets of sequences equal.
Sequence(s) can be one sequence file, directory with sequence files,
or step directory of type sequences, annotations or alignements"""

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('left', help='Left hand side sequence(s)')
        parser.add_argument('right', help='Right hand side sequence(s)')
        parser.add_argument('-e', '--extensions', help='Files with only given extensions. Format ext1;ext2')

    def run(self):
        from ..utils.helpers import read_raw_sequences_from_all
        extensions = self.args.extensions.split(';') if self.args.extensions else None
        errors = []
        left = dict(read_raw_sequences_from_all(self.args.left, extensions=extensions))
        # right = dict(read_raw_sequences_from_all(self.args.right, extensions=extensions))
        # print(sorted(right.keys()))
        not_in_right = set(left.keys())
        for seq_ident, r_seq in read_raw_sequences_from_all(self.args.right, extensions=extensions):
            l_seq = left.get(seq_ident)
            if l_seq:
                not_in_right.discard(seq_ident)
                if r_seq != l_seq:
                    errors.append(f'Differs {seq_ident}, left {len(l_seq)}, right {len(r_seq)}!')
            else:
                errors.append(f'Not in left {seq_ident}!')
        if not_in_right:
            errors.append(f'Not in right {", ".join(sorted(not_in_right))}!')
        if errors:
            print('ERRORS:')
            for e in errors:
                print(e)


class SequencesSumLength(NonProjectCommand):
    # This one is better: https://github.com/MikeTrizna/assembly_stats
    _COMMAND = 'sequences_summary'
    _HELP = "Sum length of sequences"

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('filenames', nargs='+', help='Sequence filename(s)')

    def _stat_from_lengths(self, lengths, gc, ns):
        sorted_lengths = numpy.array(sorted(lengths))
        sum_l = int(numpy.sum(sorted_lengths))

        csum = numpy.cumsum(sorted_lengths)
        nx = int(sum_l // 2)
        csumn = min(csum[csum >= nx])
        l_50 = int(numpy.where(csum == csumn)[0])
        n_50 = int(sorted_lengths[l_50])

        return dict(sum_length=sum_l,
                    max_length=int(sorted_lengths[-1]),
                    min_length=int(sorted_lengths[0]),
                    num_parts=len(lengths),
                    median=int(numpy.median(sorted_lengths)),
                    mean=numpy.mean(sorted_lengths),
                    L50=l_50, N50=n_50,
                    gc=100 * (gc / (sum_l - ns)))

    def run(self):
        from ..utils.helpers import read_sequences
        lengths = []
        gc = 0
        ns = 0
        for f in self.args.filenames:
            f_lengths = []
            f_gc = 0
            f_ns = 0
            for seq in read_sequences(f):
                f_lengths.append(len(seq))
                f_gc += seq.seq.count('G') + seq.seq.count('C')
                f_ns += seq.seq.count('N')
            lengths.extend(f_lengths)
            gc += f_gc
            ns += f_ns
            #
            d = self._stat_from_lengths(f_lengths, f_gc, f_ns)
            print(f"  File {f}: {d['num_parts']} sequences, of length {d['sum_length']}, " +
                  f"max length {d['max_length']}, N50 {d['N50']}, GC {d['gc']:.2f}%.")
        #
        d = self._stat_from_lengths(lengths, gc, ns)
        print(f"All: {d['num_parts']} sequences, of length {d['sum_length']}, " +
              f"max length {d['max_length']}, N50 {d['N50']}, GC {d['gc']:.2f}%.")
