import os.path
from step_project.base_step import Step, StepCollection
from common_utils.exceptions import ZCItoolsValueError
from common_utils.file_utils import read_fasta_identifiers, write_lines_in_file
from common_utils.misc import sets_equal
from common_utils.terminal_layout import TableByColumns
from ..utils.import_methods import import_bio_seq_io, import_bio_align_io


class AlignmentStep(Step):
    """
Stores an alignment between 2 or more sequences.
Data is stored:
 - desription.yml     : list of sequence identifiers to align and type of these sequences.
 - sequences.fa       : sequences to align
 - alignments.phy|fa  : alignment
 # For phylogenetic analysis, storing information of sequence parts/features
 - <seq_ident>._parts : sequence partitions
        If alignment is on only one feature that these files are not stored.
        File contains lines of format: <end index> <description>.
        Indices are zero indexed and means that indices with 'last_end <= indices < end' are in part
        Description doesn't have to be unique inside a file.
        This should be real partition. Covers whole sequence
        Note: extension ._parts is used since MrBayes has output file with extension .parts
 - <seq_ident>.gb : annotated sequence
        Same sequnece with filtered features
"""
    # ToDo: other alignment formats?
    _STEP_TYPE = 'alignment'
    _SEQUENCE_TYPES = ('gene', 'genes', 'whole')

    # Init object
    def _init_data(self, type_description):
        self._sequences = set()  # seq_ident
        self._seq_type = None    # None or string from _SEQUENCE_TYPES
        if type_description:
            self._sequences.update(type_description['sequences'])
            self._seq_type = type_description['sequence_type']

    def _check_data(self):
        # Check does alignment file exist
        if all(not os.path.isfile(self.step_file(a)) for a in ('alignment.phy', 'alignment.fa')):
            raise ZCItoolsValueError(f'No alignment file for step {self.directory}!')

        # Check are all sequences in sequence file
        exist_seq_idents = set(read_fasta_identifiers(self.step_file('sequences.fa')))
        sets_equal(self._sequences, exist_seq_idents, 'sequence', step=self.directory)

        # Check are all sequences in alignment file
        # ToDo: ...

    # Set data
    def set_sequences(self, seqs):
        self._sequences.update(seqs)

    def seq_sequence_type(self, t):
        if t is not None and t not in self._SEQUENCE_TYPES:
            raise ZCItoolsValueError(f'Wrong alignment sequence type ({t})!')
        self._seq_type = t

    # Save/load data
    def save(self, create=True, completed=True):
        # Store description.yml
        self.save_description(dict(sequences=sorted(self._sequences), sequence_type=self._seq_type),
                              create=create, completed=completed)

    def store_partition(self, seq_ident, partition):
        assert seq_ident in self._sequences, seq_ident
        write_lines_in_file(self.step_file(f'{seq_ident}._parts'), [f'{n} {desc}' for n, desc in partition])

    # Retrieve data methods
    def all_sequences(self):
        return self._sequences

    def get_sequence_type(self):
        return self._seq_type

    def is_short(self):
        return self._seq_type == 'gene'

    def get_alignment_obj(self):
        for f, t in (('alignment.phy', 'phylip-relaxed'), ('alignment.fa', 'fasta')):
            f = self.step_file(f)
            if os.path.isfile(f):
                with open(f) as in_f:
                    return import_bio_align_io().read(in_f, t)

    def get_phylip_file(self):
        f = self.step_file('alignment.phy')
        if not os.path.isfile(f):
            fa_f = self.step_file('alignment.fa')
            AlignIO = import_bio_align_io()
            p_format = 'phylip-relaxed' if any(len(s) > 10 for s in all_sequences) else 'phylip'
            AlignIO.convert(fa_f, 'fasta', f, p_format)
        return f

    def get_nexus_file(self):
        f = self.step_file('alignment.nex')
        if os.path.isfile(f):
            return f
        if p_f := self.get_phylip_file():
            import_bio_align_io().convert(p_f, 'phylip-relaxed', f, 'nexus', molecule_type='DNA')
            return f

    def make_summary_data(self):
        align = self.get_alignment_obj()
        # ToDo: some more data
        return dict(alignment_length=len(align[0]))

    # Partitions (used for phylogeny analysis)
    def get_partition_from_file(self, seq_ident):
        # Returns None or list of integers
        try:
            ret = []
            with open(self.step_file(f'{seq_ident}._parts'), 'r') as r:
                for line in r.readlines():
                    n, d = line.strip().split(' ', maxsplit=1)
                    ret.append((int(n), d))
            return ret
        except OSError:
            return

    def get_alignment_map_indices(self, with_reverse=False):
        from .alignment_map_indices import AlignmentMapIndices
        return AlignmentMapIndices(filename=self.step_file('alignment.phy'), with_reverse=with_reverse)

    # Show data
    def show_data(self, params=None):
        print('Alignment', self.directory, self._seq_type)

    def diff_matrix(self):
        seqs = dict((seq.name, seq.seq) for seq in self.get_alignment_obj())
        seq_idents = sorted(seqs.keys())

        num_s = len(seq_idents)
        em = [''] * num_s
        columns = [['', ''] + seq_idents[1:]]
        for i, s1 in enumerate(seq_idents):
            print(seqs[s1][:100])
            col = [s1, sum(int(c != '-') for c in seqs[s1])] + ['-'] * i
            columns.append(col)
            for s2 in seq_idents[i + 1:]:
                col.append(sum(int(c1 != c2) for c1, c2 in zip(seqs[s1], seqs[s2])))

        print(TableByColumns(columns))


class AlignmentsStep(StepCollection):
    _STEP_TYPE = 'alignments'
    _SUBSTEP_CLASS = AlignmentStep


# ToDo:
class mVISTAStep(AlignmentStep):
    """
...
"""
    _STEP_TYPE = 'mvista'
