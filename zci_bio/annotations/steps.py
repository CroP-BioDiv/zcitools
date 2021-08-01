import os.path
from collections import defaultdict, Counter
from step_project.base_step import Step
from common_utils.misc import sets_equal
from common_utils.show import print_table
from common_utils.file_utils import silent_remove_file
from ..utils.import_methods import import_bio_seq_io
from ..utils.helpers import feature_qualifiers_to_desc, feature_location_desc, concatenate_sequences, fix_sequence


class AnnotationsStep(Step):
    """
Stores list of (DNA) sequences with there annotations.
List of sequence identifier are stored in description.yml.
Annotations are stored:
 - in file annotations.gb, for whole sequnece set,
 - or in files <seq_ident>.gb for each sequence separately.
"""
    _STEP_TYPE = 'annotations'

    # Init object
    def _init_data(self, type_description):
        self._sequences = set()  # seq_ident
        if type_description:
            self._sequences.update(type_description['sequences'])
        self._cached_seq_recs = dict()

    def _check_data(self):
        exist_seq_idents = set(seq_ident for seq_ident, _ in self._iterate_records())
        sets_equal(self._sequences, exist_seq_idents, 'sequence', step=self.directory)

    # Set data
    def set_sequences(self, seqs):
        self._sequences.update(seqs)

    def add_sequence_file(self, f):
        # Filename is relative inside step directory
        seq_ident, ext = os.path.splitext(f)
        if ext != '.gb':
            raise ZCItoolsValueError(f"Extension '{ext}' is not known sequence format! {f}")
        self._sequences.add(seq_ident)

    # Save/load data
    def save(self, create=True, completed=True, additional_data=None):
        # Store description.yml
        data = dict(sequences=sorted(self._sequences))
        if additional_data:
            data.update(additional_data)
        self.save_description(data, create=create, completed=completed)

    # Retrieve data methods
    def all_sequences(self):
        return self._sequences

    def has_sequences(self):
        return bool(self._sequences)

    def get_sequence_filename(self, seq_ident):
        return self.step_file(seq_ident + '.gb')

    def get_sequence_record(self, seq_ident, cache=True):
        if not (seq_rec := self._cached_seq_recs.get(seq_ident)):
            with open(self.get_sequence_filename(seq_ident), 'r') as in_s:
                seq_record = import_bio_seq_io().read(in_s, 'genbank')
                assert seq_ident == seq_record.id.split('.', 1)[0], (seq_ident, seq_record.id)
                seq_rec = fix_sequence(seq_record)
                if cache:
                    self._cached_seq_recs[seq_ident] = seq_rec
        return seq_rec

    def get_sequence(self, seq_ident):
        return self.get_sequence_record(seq_ident).seq

    def _iterate_records(self, filter_seqs=None):
        for seq_ident in sorted(filter_seqs or self._sequences):
            yield seq_ident, self.get_sequence_record(seq_ident)

    def concatenate_seqs_genbank(self, filename, seq_idents, filter_seqs=None):
        concatenate_sequences(filename, [self.get_sequence_filename(si) for si in filter_seqs or self._sequences])

    def can_be_completed(self):
        # Note: Work for GeSeq!
        return any(f.endswith('.zip') for f in os.listdir(self.directory))

    #
    def extract_shared_features(self, feature_type, filter_seqs=None):
        # Returns tuple (set of same gene, dict ((seq_ident, gene) -> seq part))
        data = self._get_feature_desc(feature_type, filter_seqs=filter_seqs)
        if len(data) <= 1:
            print('Not enough data to find shared features!')
            return

        same_genes = self._intersect_features(data)
        if not same_genes:
            print('There is no shared features in annotations!')
            return

        parts = defaultdict(str)  # (seq_ident, feature name) -> sequence
        for seq_ident in data.keys():
            seq_record = self.get_sequence_record(seq_ident)
            for f in seq_record.features:
                if f.type == feature_type:
                    name = feature_qualifiers_to_desc(f)
                    if name in same_genes:
                        parts[(seq_ident, name)] += str(f.extract(seq_record).seq)
        #
        return same_genes, parts

    def get_features(self, feature_type, filter_seqs=None):
        # Iterate through sequences and there genes
        for seq_ident, seq_record in self._iterate_records(filter_seqs=filter_seqs):
            yield seq_ident, seq_record, (f for f in seq_record.features if f.type == feature_type)

    def _get_feature_desc(self, feature_type, filter_seqs=None):
        # Returns dict seq_ident -> dict (gene name -> num occurences)
        return dict((seq_ident, Counter(feature_qualifiers_to_desc(f) for f in features))
                    for seq_ident, _, features in self.get_features(feature_type, filter_seqs=filter_seqs))

    # Show data
    def show_data(self, params=None):
        # If listed, filter only these sequences
        filter_seqs = None
        if params:
            filter_seqs = self._sequences & set(params)
            params = [p for p in params if p not in filter_seqs]  # Remove processed params

        cmd = params[0].lower() if params else 'by_type'  # Default print
        if params:
            params = params[1:]

        if cmd == 'by_type':
            all_types = set()
            data = dict()  # seq_ident -> dict(length=int, features=int, <type>=str)
            for seq_ident, seq_record in self._iterate_records(filter_seqs=filter_seqs):
                type_2_list = defaultdict(list)
                for f in seq_record.features:
                    if f.type != 'source':
                        type_2_list[f.type].append(feature_qualifiers_to_desc(f))
                all_types.update(type_2_list.keys())
                data[seq_ident] = dict(length=len(seq_record.seq), features=len(seq_record.features))
                data[seq_ident].update((t, f"{len(fs)}/{len(set(fs))}" if t != 'repeat_region' else str(len(fs)))
                                       for t, fs in type_2_list.items())

            all_types = sorted(all_types)
            print_table(['seq_ident', 'Length', 'Features'] + all_types,
                        sorted([seq_ident, d['length'], d['features']] + [d.get(t, '') for t in all_types]
                               for seq_ident, d in sorted(data.items())))

        # Genes
        elif cmd == 'genes':
            self._all_features(self._get_feature_desc('gene', filter_seqs=filter_seqs))
        elif cmd == 'repeated_genes':
            self._repeated_features(self._get_feature_desc('gene', filter_seqs=filter_seqs))
        elif cmd == 'shared_genes':
            self._shared_features(self._get_feature_desc('gene', filter_seqs=filter_seqs))

        # CDSs
        elif cmd == 'cds':
            self._all_features(self._get_feature_desc('CDS', filter_seqs=filter_seqs))
        elif cmd == 'repeated_cds':
            self._repeated_features(self._get_feature_desc('CDS', filter_seqs=filter_seqs))
        elif cmd == 'shared_cds':
            self._shared_features(self._get_feature_desc('CDS', filter_seqs=filter_seqs))

        # IR
        elif cmd == 'ir':
            # seq_ident -> list of feature locations
            data = dict((
                seq_ident,
                [feature_location_desc(f.location)
                 for f in seq_record.features if f.type == 'repeat_region'])
                for seq_ident, seq_record in self._iterate_records(filter_seqs=filter_seqs))

            for seq_ident, locations in sorted(data.items()):
                print(f"{seq_ident} ({len(locations)}): {', '.join(map(str, sorted(locations)))}")
        else:
            print(f'Wrong show command ({cmd})!')

    def _all_features(self, data):
        # data is dict seq_ident -> dict (gene -> num occurences)
        for seq_ident, genes in sorted(data.items()):
            print(f"{seq_ident} ({len(genes)}): {', '.join(sorted(genes.keys()))}")

    def _repeated_features(self, data):
        # data is dict seq_ident -> dict (gene -> num occurences)
        for seq_ident, genes in sorted(data.items()):
            num_2_genes = defaultdict(list)
            num_repeted = 0
            for g, num in genes.items():
                if num > 1:
                    num_2_genes[num].append(g)
                    num_repeted += 1
            if num_2_genes:
                print(f"{seq_ident} ({num_repeted}):")
                for num, genes in sorted(num_2_genes.items(), reverse=True):
                    print(f"    {num}: {', '.join(sorted(genes))}")

    def _shared_features(self, data):
        # data is dict seq_ident -> dict (gene -> num occurences)
        if len(data) > 1:
            same_features = self._intersect_features(data)
            print('Genes not shared by all sequences:')
            for seq_ident, genes in sorted(data.items()):
                rest_features = set(genes.keys()) - same_features
                print(f"    {seq_ident} ({len(rest_features)}): {', '.join(sorted(rest_features))}")

            print(f"Shared ({len(same_features)}): {', '.join(sorted(same_features))}")
        else:
            print('Not enough data to find same genes!')

    def _intersect_features(self, data):
        # data is dict seq_ident -> dict (gene -> num occurences)
        return set.intersection(*(set(g.keys()) for g in data.values()))
