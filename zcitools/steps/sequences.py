import os.path
from collections import defaultdict
from .step import Step
from ..utils.exceptions import ZCItoolsValueError
from ..utils.import_methods import import_bio_seq_io


class SequencesStep(Step):
    """
Stores list of (DNA) sequences.
List of sequence identifier are stored in description.yml.
Each sequence can be stored in one or more files in different formats.
"""
    _STEP_TYPE = 'sequences'
    _KNOWN_EXTENSIONS = ('.gb', '.fa')
    _SeqIO_TYPES = ('genbank', 'fasta')

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
            si = seq_ident + '.'
            for old_f in self._sequences[seq_ident]:
                if old_f != f:
                    silent_remove_file(self.step_file(old_f))
        else:
            self._sequences[seq_ident] = [f]
        #
        self.remove_cache_files()

    # Save/load data
    def save(self):
        # Store description.yml
        self.save_description(dict(sequences=sorted(self._sequences)))
        # Data files are handled with add_sequence_file() method

    # Retrieve data methods
    def sequence_exists(self, ident):
        return bool(self._sequences.get(ident))

    def all_sequences(self):
        return self._sequences.keys()

    # Cach files are prfixed with '_c_'
    def _iterate_seq_records(self):
        SeqIO = import_bio_seq_io()

        # Iterate through all sequences, returns Bio.SeqRecord objects.
        for seq_ident, files in sorted(self._sequences.items()):
            print(seq_ident, files)
            for ext, st in zip(self._KNOWN_EXTENSIONS, self._SeqIO_TYPES):
                f = seq_ident + ext
                if f in files:
                    print('  eva', f)
                    with open(self.step_file(f), 'r') as in_s:
                        yield from SeqIO.parse(in_s, st)
                    break

    def get_all_seqs_fa(self):
        f = self.step_file('_c_all_seqs.fa')
        if not os.path.isfile(f):
            with open(f, 'w') as fa:
                for seq_record in self._iterate_seq_records():
                    fa.write(f">{seq_record.id}\n{seq_record.seq}\n")
                    # fa.write(f">{seq_record.id} {seq_record.description}\n{seq_record.seq}\n")
        return f
