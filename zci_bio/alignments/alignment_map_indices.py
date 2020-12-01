from collections import defaultdict
from ..utils.helpers import read_alignment, feature_qualifiers_to_desc


class IndexInAlignment:
    # Maps original index into alined
    # Note: indices are zero indexed!
    def __init__(self, seq):
        self.indices = [i for i, c in enumerate(seq.seq) if c != '-']

    def __getitem__(self, i):
        return self.indices[i]

    def __len__(self):
        return len(self.indices)

    def all_in_feature(self, feature):
        return set.union(*(set(range(self.indices[p.start], self.indices[p.end - 1] + 1))
                           for p in feature.location.parts))


class AlignmentMapIndices:
    def __init__(self, alignment=None, filename=None):
        assert alignment or filename
        self.alignment = alignment or read_alignment(filename)
        self.seq_indices = dict((seq.name, IndexInAlignment(seq)) for seq in self.alignment)

    def map_partitions(self, data):
        num_ps = len(next(iter(orig_parts.values())))
        # All have to have same number of parts
        assert all(len(p) == num_ps for p in data.values()), [(s, len(p)) for s, p in data.items()]
        # All have to be indexed
        assert not (not_in := [s for s in data.keys() if s not in self.seq_indices]), not_in


    def create_raxml_partitions_from_parts(self, orig_parts, partitions_filename):
        # orig_parts is iterable of tuples (seq_ident, partitions)
        orig_parts = dict(orig_parts)  # seq_ident -> list of tuples (index, description)
        if not orig_parts:
            return

        num_genes = len(next(iter(orig_parts.values())))
        assert all(len(p) == num_genes for p in orig_parts.values()), \
            [(s, len(p)) for s, p in orig_parts.items()]
        assert all(s in self.seq_indices for s in orig_parts.keys()), \
            [s for s in orig_parts.keys() if s not in self.seq_indices]

        #
        with open(partitions_filename, 'w') as output:
            last_start = 1  # RAxML indexes form 1

            for i in range(num_genes):                       # For each gene
                for seq_ident, parts in orig_parts.items():  # Check between all sequences
                    aligned = self.seq_indices[seq_ident]    # Is there start/end pair aligned on neighbouring positions
                    end, gene = parts[i]
                    if end >= len(aligned):
                        output.write(f'DNA, {gene} = {last_start}-{len(self.alignment[0])}\n')
                        assert i == (num_genes - 1), (i, num_genes)
                        break
                    p = aligned[end - 1]  # In current part
                    n = aligned[end]      # In next part
                    if n == p + 1:
                        output.write(f'DNA, {gene} = {last_start}-{n}\n')
                        last_start = n + 1  # Indices are zero indexed!
                        break
                else:
                    raise ValueError(parts[i][1])

    def create_raxml_partitions_from_features(self, orig_features, partitions_filename):
        # orig_features is iterable of tuples (seq_ident, seq_record with filtered features)
        def _group_features(features):
            ret = defaultdict(list)
            for f in seq_record.features:
                ret[feature_qualifiers_to_desc(f)].append(f)
            return ret

        orig_features = dict((seq_ident, _group_features(seq_record.features))
                             for seq_ident, seq_record in orig_features)
        if not orig_features:
            return

        # Find shared features
        shared_features = set.intersection(*(set(g.keys()) for fs in orig_features.values()))

        final_features = []  # Tuples (from, to)
        for f_name in shared_features:
            # ToDo: Slower but safer way with set of feature's indices
            # Calculate feature indices in alignment
            f_locs = dict((seq_ident, [aligned.all_in_feature(f) for f in orig_features[seq_ident][f_name]])
                          for seq_ident, aligned in self.seq_indices.items())

            # Find intersection region in alignment
            indices_in = set.intersection(*(set.union(*locs) for locs in f_locs.values()))
            if not indices_in:
                continue

            # Find which original features contain intersections
            if indices_in:
                all_in = set.union(*(set.union(*(f for f in seq_features if not indices_in.isdisjoin(f)))
                                     for seq_features in f_locs.values()))

                final_features.extend((s, e, f_name) for s, e in _in_ranges(sorted(all_in)))
        final_features.sort()


def _in_ranges(sorted_list):
    # Converts list of sorted indices into ranges, inclusive
    last_start = last_end = sorted_list[0]
    for i in sorted[1:]:
        if i != last_end + 1:
            yield last_start, last_end
            last_start = i
        last_end = i
    yield last_start, last_end
