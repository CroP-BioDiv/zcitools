from .step import Step


class AnnotationsStep(Step):
    """
Stores list of (DNA) sequences with there annotations.
List of sequence identifier are stored in description.yml.
Each annotation (with sequence) is stored in separate file?
"""
    ""
    _STEP_TYPE = 'annotations'
    # _KNOWN_EXTENSIONS = frozenset(['.gb', '.fa'])

    # Init object
    def _init_data(self, type_description):
        self._sequences = set()  # seq_ident
        if type_description:
            self._sequences.update(type_description['sequences'])
        #     existing_seqs = self._find_existing_seqs()
        #     for seq_ident in type_description['sequences']:
        #         self._sequences[seq_ident] = existing_seqs.get(seq_ident, [])

    def _check_data(self):
        pass

    # Set data
    def set_sequences(self, seqs):
        self._sequences.update(seqs)

    # Save/load data
    def save(self):
        # Store description.yml
        self.save_description(dict(sequences=sorted(self._sequences)))
