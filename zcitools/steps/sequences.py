import os.path
from collections import defaultdict
from .step import Step
from ..utils.exceptions import ZCItoolsValueError
from ..utils.import_methods import import_bio_seq_io
from ..utils.terminal_layout import StringColumns
from ..utils.helpers import write_fasta


class SequencesStep(Step):
    """
Stores list of (DNA) sequences.
List of sequence identifier are stored in description.yml.
Each sequence can be stored in one or more files in different formats.
"""
    _STEP_TYPE = 'sequences'
    _TYPE_LOAD_PRIORITY = (('.gb', 'genbank'), ('.fa', 'fasta'))
    _KNOWN_EXTENSIONS = tuple(e for e, _ in _TYPE_LOAD_PRIORITY)
    _EXT_2_TYPE = dict(_TYPE_LOAD_PRIORITY)
    _TYPE_2_EXT = dict(x[::-1] for x in _TYPE_LOAD_PRIORITY)

    # Init object
    def _init_data(self, type_description):
        self._sequences = dict()  # seq_ident -> list of files
        if type_description:
            existing_seqs = self._find_existing_seqs()
            for seq_ident in type_description['sequences']:
                self._sequences[seq_ident] = existing_seqs.get(seq_ident, [])

    def _check_data(self):
        existing_seqs = self._find_existing_seqs()
        exist_seq_idents = set(existing_seqs.keys())
        needed_seq_idents = set(self._sequences.keys())
        # Are all sequences presented
        not_exist = needed_seq_idents - exist_seq_idents
        if not_exist:
            raise ZCItoolsValueError(f"Sequence data not presented for: {', '.join(sorted(not_exist))}")

        # Is there more sequences
        more_data = exist_seq_idents - needed_seq_idents
        if more_data:
            raise ZCItoolsValueError(f"Data exists for not listed sequence(s): {', '.join(sorted(more_data))}")

    def _find_existing_seqs(self):
        existing_seqs = defaultdict(list)
        for f in self.step_files(not_cached=True):
            for e in self._KNOWN_EXTENSIONS:
                if f.endswith(e):
                    existing_seqs[f[:-len(e)]].append(f)
                    break
        return existing_seqs

    # Set data
    def add_sequence_file(self, f):
        # Filename is relative inside step directory
        seq_ident, ext = os.path.splitext(f)
        if ext not in self._KNOWN_EXTENSIONS:
            raise ZCItoolsValueError(f"Extension '{ext}' is not known sequence format! {f}")

        if seq_ident in self._sequences:
            # Remove other (old) files of same sequence
            for old_f in self._sequences[seq_ident]:
                if old_f != f:
                    silent_remove_file(self.step_file(old_f))
        else:
            self._sequences[seq_ident] = [f]
        #
        self.remove_cache_files()

    # Save/load data
    def save(self, create=True, needs_editing=False):
        # Store description.yml
        self.save_description(dict(sequences=sorted(self._sequences)), create=create, needs_editing=needs_editing)
        # Data files are handled with add_sequence_file() method

    # Retrieve data methods
    def sequence_exists(self, ident):
        return bool(self._sequences.get(ident))

    def all_sequences(self):
        return self._sequences.keys()

    def _iterate_records(self):
        # Iterate through all sequences, returns Bio.SeqRecord objects.
        for seq_ident, files in sorted(self._sequences.items()):
            seq_record = self._read_record(seq_ident, files=files)
            if seq_record is not None:
                yield seq_ident, seq_record

    def _read_record(self, seq_ident, files=None):
        if files is None:
            files = self._sequences[seq_ident]
        for ext, st in self._TYPE_LOAD_PRIORITY:
            f = seq_ident + ext
            if f in files:
                with open(self.step_file(f), 'r') as in_s:
                    return import_bio_seq_io().read(in_s, st)

    # def get_all_seqs_fa(self):
    #     f = self.cache_file('all_seqs.fa')
    #     if not os.path.isfile(f):
    #         # Note: sequences are named by initial seq_ident (not seq_record.id)!
    #         write_fasta(f, ((seq_ident, seq_record.seq) for seq_ident, seq_record in self._iterate_records()))
    #     return f

    def concatenate_seqs_fa(self, filename, seq_idents):
        write_fasta(filename, ((seq_ident, self._read_record(seq_ident).seq) for seq_ident in seq_idents))

    def get_sequence(self, seq_ident):
        # Returns sequence as a string
        seq_record = self._read_record(seq_ident)
        return str(seq_record.seq)

    def get_sequence_file(self, seq_ident, file_type):
        # Returns sequence of given file type
        # First check does it exist
        ext = self._TYPE_2_EXT[file_type]
        for f in self._sequences[seq_ident]:
            if f.endswith(ext):
                return self.step_file(f)

        # If not, than make convertion
        seq_record = self._read_record(seq_ident)
        seq_f = self.step_file(seq_ident + ext)
        with open(seq_f, 'w') as output_handle:
            count = import_bio_seq_io().write(seq_record, output_handle, file_type)
        self._sequences[seq_ident].append(seq_ident + ext)
        return seq_f

    # Show data
    def show_data(self, params=None):
        header = ['seq_ident', 'Record ID', 'Length', 'Num features']
        rows = [[seq_ident, seq_record.id, len(seq_record.seq), len(seq_record.features)]
                for seq_ident, seq_record in self._iterate_records()]
        print(StringColumns(rows, header=header))
