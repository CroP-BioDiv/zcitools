import os.path
from step_project.base_step import Step, StepCollection
from common_utils.exceptions import ZCItoolsValueError
from common_utils.file_utils import read_fasta_identifiers, write_lines_in_file
from common_utils.misc import sets_equal
from ..utils.import_methods import import_bio_seq_io, import_bio_align_io, import_bio_alphabet


class AlignmentStep(Step):
    """
Stores an alignment between 2 or more sequences.
Data is stored:
 - desription.yml    : list of sequence identifiers to align and type of these sequences.
 - sequences.fa      : sequences to align
 - alignments.phy    : alignment
 - <seq_ident>.parts : sequence partitions. If alignment is on only one feature that these files are not stored.
                       File has lines of format: <end index> <description>.
                       Note: this should be real partition. Covers whole sequence
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
        if not os.path.isfile(self.step_file('alignment.phy')):
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

    def store_partitions(self, seq_ident, partitions):
        assert seq_ident in self._sequences, seq_ident
        write_lines_in_file(self.step_file(f'{seq_ident}.parts'), [f'{n} {desc}' for n, desc in partitions])

    # Retrieve data methods
    def all_sequences(self):
        return self._sequences

    def get_sequence_type(self):
        return self._seq_type

    def is_short(self):
        return self._seq_type == 'gene'

    def is_composite(self):  # Contains more genes
        return self._seq_type == 'genes'

    def get_phylip_file(self):
        f = self.step_file('alignment.phy')
        if os.path.isfile(f):
            return f
        # ToDo: when other formats are in, change this

    def get_nexus_file(self):
        f = self.step_file('alignment.nex')
        if os.path.isfile(f):
            return f
        p_f = self.get_phylip_file()
        if p_f:
            AlignIO = import_bio_align_io()
            Alphabet = import_bio_alphabet()
            AlignIO.convert(p_f, 'phylip', f, 'nexus', alphabet=Alphabet.generic_dna)
            return f

    #
    def iterate_partitions(self):
        for seq_ident in self._sequences:
            p = self.get_partitions(seq_ident)
            if p:
                yield seq_ident, p

    def get_partitions(self, seq_ident):
        # Returns None or list of integers
        f = self.step_file(f'{seq_ident}.parts')
        if os.path.isfile(f):
            ret = []
            with open(f, 'r') as r:
                for line in r.readlines():
                    n, d = line.strip().split(' ', maxsplit=1)
                    ret.append((int(n), d))
            return ret

    def get_alignment_map_indices(self):
        from .alignment_map_indices import AlignmentMapIndices
        return AlignmentMapIndices(filename=self.step_file('alignment.phy'))

    # Show data
    def show_data(self, params=None):
        print('Alignment', self.directory, self._seq_type)


class AlignmentsStep(StepCollection):
    _STEP_TYPE = 'alignments'
    _SUBSTEP_CLASS = AlignmentStep


# ToDo:
class mVISTAStep(AlignmentStep):
    """
...
"""
    _STEP_TYPE = 'mvista'
