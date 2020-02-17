import os.path
from collections import defaultdict, Counter
from step_project.base_step import Step
from common_utils.misc import sets_equal
from common_utils.show import print_table
from common_utils.cache import cache_args
from common_utils.terminal_layout import StringColumns
from ..utils.import_methods import import_bio_seq_io
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
        sets_equal(self._sequences, exist_seq_idents, 'sequence', step=self.directory)

    # Set data
    def set_sequences(self, seqs):
        self._sequences.update(seqs)

    # Save/load data
    def save(self, create=True, completed=True):
        # Store description.yml
        self.save_description(dict(sequences=sorted(self._sequences)), create=create, completed=completed)

    # Retrieve data methods
    def all_sequences(self):
        return self._sequences

    @cache_args
    def get_sequence_record(self, seq_ident):
        with open(self.step_file(seq_ident + '.gb'), 'r') as in_s:
            seq_record = import_bio_seq_io().read(in_s, 'genbank')
            assert seq_ident == seq_record.id, (seq_ident, seq_record.id)
        return seq_record

    def get_sequence(self, seq_ident):
        return self.get_sequence_record(seq_ident).seq

    def _iterate_records(self, filter_seqs=None):
        for seq_ident in filter_seqs or self._sequences:
            yield seq_ident, self.get_sequence_record(seq_ident)

    def concatenate_seqs_genbank(self, filename, seq_idents):
        concatenate_sequences(filename, [self.step_file(si + '.gb') for si in filter_seqs or self._sequences])

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
            data = dict((seq_ident,
                         [feature_location_desc(f.location) for f in seq_record.features if f.type == 'repeat_region'])
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

    #
    def chloroplast_annotation(self, num_genes=1, feature_type='gene', features=None, sequences=None):
        sequences = (self._sequences & set(sequences)) if sequences else self._sequences
        rows = []
        for seq_ident in sorted(sequences):
            seq = self.get_sequence_record(seq_ident)
            rep_regs = [f for f in seq.features if f.type == 'repeat_region']
            if len(rep_regs) != 2:
                print(f"Sequence {seq_ident} doesn't have inverted repeats!")
                continue
            ira, irb = rep_regs
            if len(seq) != irb.location.end:
                print(f"  warning ({seq_ident}): sequence's IRB doesn't end on sequence end!")

            locs = [(1, ira.location.start),
                    (ira.location.start, ira.location.end),
                    (ira.location.end, irb.location.start)]
            regions = [[] for _ in range(4)]
            borders = [[] for _ in range(4)]
            for f in sorted((f for f in seq.features if f.type == feature_type and f.location),
                            key=lambda x: x.location.start):
                ls = f.location.start
                le = f.location.end
                if le < ls:
                    borders[0].append(f)
                else:
                    for bi, (i1, i2) in enumerate(locs):
                        if le <= i2:
                            regions[bi].append(f)
                            break
                        elif ls < i2:
                            borders[bi + 1].append(f)
                            break
                    else:
                        regions[-1].append(f)
            #
            header = ['|', 'LSC', '  |', 'IRa', '  |', 'SSC', '  |', 'IRb']
            indices = ['', '', str(ira.location.start), '', str(ira.location.end), '', str(irb.location.start), '']
            row = []

            for i in range(4):
                row.append(borders[i][0].qualifiers['gene'][0] if borders[i] else '')
                rs = regions[i]
                if len(rs) <= 2 * num_genes:
                    row.append(' '.join(r.qualifiers['gene'][0] for r in rs))
                else:
                    row.append(' '.join(r.qualifiers['gene'][0] for r in rs[:num_genes]) + ' ... ' +
                               ' '.join(r.qualifiers['gene'][0] for r in rs[-num_genes:]))

            #
            rows.append([seq_ident] + [''] * 7)
            rows.extend([indices, header, row])
        #
        print(StringColumns(rows))
