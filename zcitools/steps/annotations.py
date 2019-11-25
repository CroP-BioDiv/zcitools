import os.path
from .step import Step
from ..utils.import_methods import import_bio_seq_io
from ..utils.terminal_layout import StringColumns


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

    def _check_data(self):
        exist_seq_idents = set(seq_ident for seq_ident, _ in self._iterate_records())
        # Are all sequences presented
        not_exist = self._sequences - exist_seq_idents
        if not_exist:
            raise ZCItoolsValueError(f"Sequence data not presented for: {', '.join(sorted(not_exist))}")

        # Is there more sequences
        more_data = exist_seq_idents - self._sequences
        if more_data:
            raise ZCItoolsValueError(f"Data exists for not listed sequence(s): {', '.join(sorted(more_data))}")

    # Set data
    def set_sequences(self, seqs):
        self._sequences.update(seqs)

    # Save/load data
    def get_all_annotation_filename(self):
        return self.step_file(self._ALL_FILENAME)

    def save(self, needs_editing=False):
        # Store description.yml
        self.save_description(dict(sequences=sorted(self._sequences)), needs_editing=needs_editing)

    # Retrieve data methods
    def _iterate_records(self):
        SeqIO = import_bio_seq_io()

        all_f = self.get_all_annotation_filename()
        if os.path.isfile(all_f):
            with open(all_f, 'r') as in_s:
                for seq_record in SeqIO.parse(in_s, 'genbank'):
                    yield seq_record.id, seq_record
        else:
            raise 'Not implemented!!!'

    # Show data
    def show_data(self, format=None):
        header = ['seq_ident', 'Record ID', 'Length', 'Num features']
        rows = [[seq_ident, seq_record.id, len(seq_record.seq), len(seq_record.features)]
                for seq_ident, seq_record in self._iterate_records()]
        print(StringColumns(rows, header=header))
