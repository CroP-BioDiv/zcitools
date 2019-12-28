import os.path
from step_project.base.step import Step, StepCollection
from common_utils.exceptions import ZCItoolsValueError
from common_utils.file_utils import read_fasta_identifiers
from common_utils.misc import sets_equal
from ..utils.import_methods import import_bio_align_io, import_bio_alphabet


class AlignmentStep(Step):
    """
Stores an alignment between 2 or more sequences.
Data is stored:
 - desription.yml stores list of sequence identifiers to align and type of these sequences.
 - sequences.fa stores sequences to align
 - alignments.phy stores alignment
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

    # Retrieve data methods
    def all_sequences(self):
        return self._sequences

    def get_sequence_type(self):
        return self._seq_type

    def is_short(self):
        return self._seq_type == 'gene'

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
