import os.path
from collections import defaultdict
from .step import Step
from ..utils.import_methods import import_bio_seq_io
from ..utils.show import print_table
from ..utils.helpers import feature_qualifiers_to_desc, feature_location_desc, concatenate_sequences


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
    def save(self, create=True, needs_editing=False):
        # Store description.yml
        self.save_description(dict(sequences=sorted(self._sequences)), create=create, needs_editing=needs_editing)

    # Retrieve data methods
    def all_sequences(self):
        return self._sequences

    def concatenate_seqs_genbank(self, filename, seq_idents):
        concatenate_sequences(filename, [f for _, f in self._iterate_files(filter_seqs=seq_idents)])

    def _iterate_files(self, filter_seqs=None):
        for seq_ident in sorted(self._sequences):
            if not filter_seqs or seq_ident in filter_seqs:
                yield seq_ident, self.step_file(seq_ident + '.gb')

    def _iterate_records(self, filter_seqs=None):
        SeqIO = import_bio_seq_io()
        for seq_ident, seq_filename in self._iterate_files(filter_seqs=filter_seqs):
            with open(seq_filename, 'r') as in_s:
                for seq_record in SeqIO.parse(in_s, 'genbank'):
                    assert seq_ident == seq_record.id, (seq_ident, seq_record.id)
                    yield seq_record.id, seq_record

    #
    def _extract_features(self, iter_features):
        for seq_ident, seq_record, features in iter_features:
            # ToDo: how to sort them. For now by name, it is possible to sort by location
            features = sorted(features, key=feature_qualifiers_to_desc)
            seq = ''.join(str(f.extract(seq_record).seq) for f in features)
            # if len(seq) % 3:
            #     print(seq_ident, len(seq))
            #     # print([len(s) for s in (f.extract(seq_record).seq for f in features) if len(s) % 3])
            #     for f in features:
            #         if 'translation_input' in f.qualifiers:
            #             t = ''.join(f.qualifiers['translation_input']).replace(' ', '')
            #             print(feature_qualifiers_to_desc(f), len(f.extract(seq_record).seq), len(t))
            #             print(t)
            #             # print(f.qualifiers['translation_input'])
            yield seq_ident, seq

    def get_genes(self, filter_seqs=None):
        # Iterate through sequences and there genes
        for seq_ident, seq_record in self._iterate_records(filter_seqs=filter_seqs):
            yield seq_ident, seq_record, (f for f in seq_record.features if f.type == 'gene')

    def _get_genes_desc(self, filter_seqs=None):
        # Returns dict seq_ident -> set of gene names
        return dict((seq_ident, set(feature_qualifiers_to_desc(f) for f in features))
                    for seq_ident, _, features in self.get_genes(filter_seqs=filter_seqs))

    def extract_all_genes(self, filter_seqs=None):
        return self._extract_features(self.get_genes(filter_seqs=filter_seqs))

    def get_cds(self, filter_seqs=None):
        # Iterate through sequences and there genes
        for seq_ident, seq_record in self._iterate_records(filter_seqs=filter_seqs):
            yield seq_ident, seq_record, (f for f in seq_record.features if f.type == 'CDS')

    def _get_cds_desc(self, filter_seqs=None):
        # Returns dict seq_ident -> set of gene names
        return dict((seq_ident, set(feature_qualifiers_to_desc(f) for f in features))
                    for seq_ident, _, features in self.get_cds(filter_seqs=filter_seqs))

    def extract_all_cds(self, filter_seqs=None):
        return self._extract_features(self.get_cds(filter_seqs=filter_seqs))

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
            data = dict()  # seq_ident -> dict(length=int, features=int, <type>=num)
            for seq_ident, seq_record in self._iterate_records(filter_seqs=filter_seqs):
                d = defaultdict(int)
                genes = set()
                for f in seq_record.features:
                    if f.type != 'source':
                        d[f.type] += 1
                        if f.type == 'gene':
                            genes.add(feature_qualifiers_to_desc(f))
                if genes:
                    d['gene_unique'] = len(genes)
                all_types.update(d.keys())
                d['length'] = len(seq_record.seq)
                d['features'] = len(seq_record.features)
                data[seq_ident] = d

            all_types = sorted(all_types)
            print_table(['seq_ident', 'Length', 'Features'] + all_types,
                        sorted([seq_ident, d['length'], d['features']] + [d.get(t, 0) for t in all_types]
                               for seq_ident, d in sorted(data.items())))

        elif cmd == 'genes':
            self._all_features(self._get_genes_desc(filter_seqs=filter_seqs).items())

        elif cmd == 'shared_genes':
            self._shared_features(self._get_genes_desc(filter_seqs=filter_seqs))

        elif cmd == 'cds':
            self._all_features(self._get_cds_desc(filter_seqs=filter_seqs).items())

        elif cmd == 'shared_cds':
            self._shared_features(self._get_cds_desc(filter_seqs=filter_seqs))

        elif cmd == 'ir':
            # seq_ident -> list of feature locations
            data = dict((seq_ident,
                         [feature_location_desc(f.location) for f in seq_record.features if f.type == 'repeat_region'])
                        for seq_ident, seq_record in self._iterate_records(filter_seqs=filter_seqs))

            for seq_ident, locations in sorted(data.items()):
                print(f"{seq_ident} ({len(locations)}): {', '.join(map(str, sorted(locations)))}")
        else:
            print(f'Wrong show command ({cmd})!')

    def _all_features(self, data):
        for seq_ident, genes in sorted(data):
            print(f"{seq_ident} ({len(genes)}): {', '.join(sorted(genes))}")

    def _shared_features(self, data):
        if len(data) > 1:
            same_genes = set.intersection(*data.values())
            print('Genes not shared by all sequences:')
            for seq_ident, genes in sorted(data.items()):
                rest_genes = genes - same_genes
                print(f"    {seq_ident} ({len(rest_genes)}): {', '.join(sorted(rest_genes))}")

            print(f"Shared ({len(same_genes)}): {', '.join(sorted(same_genes))}")
        else:
            print('Not enough data to find same ganes!')
