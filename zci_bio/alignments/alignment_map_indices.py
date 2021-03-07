from collections import defaultdict
from ..utils.helpers import read_alignment, feature_qualifiers_to_desc


class IndexInAlignment:
    # Maps original index into alined
    # Note: indices are zero indexed!
    def __init__(self, seq, with_reverse=False):
        self.seq_name = seq.name
        self.seq_length = len(seq)
        self.indices = [i for i, c in enumerate(seq.seq) if c != '-']
        if with_reverse:
            self.reverse_indices = []
            idx = 0
            for c in seq.seq:
                self.reverse_indices.append(idx)
                if c != '-':
                    idx += 1
            self.reverse_indices.append(self.seq_length)  # For end+

    def __getitem__(self, i):
        assert 0 <= i < len(self.indices), (i, len(self.indices))
        return self.indices[i]

    def __len__(self):
        return len(self.indices)

    def reverse(self, i):
        return self.reverse_indices[i]

    def reverse_length(self, s, e):
        return self.reverse_indices[e] - self.reverse_indices[s]

    def all_in_feature(self, feature):
        return set.union(*(set(range(self.indices[p.start], self.indices[p.end - 1] + 1))
                           for p in feature.location.parts))


class AlignmentMapIndices:
    def __init__(self, alignment=None, filename=None, with_reverse=False):
        assert alignment or filename
        self.alignment = alignment or read_alignment(filename)
        self.seq_indices = dict((seq.name, IndexInAlignment(seq, with_reverse=with_reverse)) for seq in self.alignment)
        self.alignment_length = len(self.alignment[0])

    def all_sequences(self):
        return self.seq_indices.keys()

    def has_sequence(self, seq_ident):
        return seq_ident in self.seq_indices

    def get_alignment(self, seq_ident):
        return self.seq_indices.get(seq_ident)
