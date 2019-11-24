from .step import Step


class AnnotationsStep(Step):
    """
Stores list of (DNA) sequences with there annotations.
List of sequence identifier are stored in description.yml.
Annotations are stored:
 - in file annotations.gb, for whole sequnece set,
 - or in files <seq_ident>.gb for each sequence separately.
"""
    _STEP_TYPE = 'annotations'
    _ALL_FILENAME = 'annotations.gb'

    # Init object
    def _init_data(self, type_description):
        self._sequences = set()  # seq_ident
        if type_description:
            self._sequences.update(type_description['sequences'])
        #     existing_seqs = self._find_existing_seqs()
        #     for seq_ident in type_description['sequences']:
        #         self._sequences[seq_ident] = existing_seqs.get(seq_ident, [])

    def _check_data(self):
        # if self.step_file('annotations.gb')
        pass

    # Set data
    def set_sequences(self, seqs):
        self._sequences.update(seqs)

    # Save/load data
    def get_all_annotation_filename(self):
        return self.step_file(self._ALL_FILENAME)

    def save(self, needs_editing=False):
        # Store description.yml
        self.save_description(dict(sequences=sorted(self._sequences)), needs_editing=needs_editing)
