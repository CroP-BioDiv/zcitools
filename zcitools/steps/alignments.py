from .step import Step


class AlignmentStep(Step):
    """Stores an aligment between 2 or more sequences."""
    _STEP_TYPE = 'alignment'

    # Init object
    def _init_data(self, type_description):
        self._sequences = set()  # seq_ident
        if type_description:
            self._sequences.update(type_description['sequences'])

    # Set data
    def set_sequences(self, seqs):
        self._sequences.update(seqs)

    # Save/load data
    def save(self, create=True, needs_editing=False):
        # Store description.yml
        self.save_description(dict(sequences=sorted(self._sequences)), create=create, needs_editing=needs_editing)

    # Retrieve data methods
    def all_sequences(self):
        return self._sequences


class AlignmentsStep(Step):
    """
Stores multiple aligments between different sets of sequences.
Each alignment is step of type AligmentStep.
Alignment is identified by string.
"""
    _STEP_TYPE = 'alignments'

    # Init object
    def _init_data(self, type_description):
        pass
        # self._alignments = set()  # alignment_ident
        # if type_description:
        #     self._alignments.update(type_description['alignments'])

    # # Set data
    # def add_alignment(self, alignment):
    #     self._alignments.update(alignment)

    # Save/load data
    def save(self, create=True, needs_editing=False):
        # Store description.yml
        # self.save_description(dict(alignments=sorted(self._alignments)), create=create, needs_editing=needs_editing)
        self.save_description(dict(), create=create, needs_editing=needs_editing)

    # Retrieve data methods
    def all_alignmets(self):
        return self._alignmets
