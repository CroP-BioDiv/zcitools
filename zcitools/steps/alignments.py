import os.path
from .step import Step
from ..utils.exceptions import ZCItoolsValueError
from ..utils.file_utils import read_fasta_identifiers
from ..utils.helpers import sets_equal
from ..utils.show import print_ls_like_list


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
    def save(self, create=True, needs_editing=False):
        # Store description.yml
        self.save_description(dict(sequences=sorted(self._sequences), sequence_type=self._seq_type),
                              create=create, needs_editing=needs_editing)

    # Retrieve data methods
    def all_sequences(self):
        return self._sequences

    # Show data
    def show_data(self, params=None):
        print('Alignment', self.directory, self._seq_type)


class AlignmentsStep(Step):
    """
Stores multiple aligments between different sets of sequences.
Each alignment is step of type AligmentStep.
Alignment is identified by AligmentStep's directory name.
Note: list of alignments is not stored in description.yml.
"""
    _STEP_TYPE = 'alignments'

    # Init object
    def _init_data(self, type_description):
        pass

    def _check_data(self):
        # Check AlignmentStep objects
        for a in self.all_alignments():
            self.read_substep(a)._check_data()

    # Save/load data
    def save(self, create=True, needs_editing=False):
        # Store description.yml
        self.save_description(dict(), create=create, needs_editing=needs_editing)

    # Retrieve data methods
    def all_alignments(self):
        return self.step_subdirectories()

    # Show data
    def show_data(self, params=None):
        print_ls_like_list('Alignments', self.all_alignments(), sort=True, min_rows_to_split=20)
