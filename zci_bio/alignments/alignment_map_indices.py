from ..utils.helpers import read_alignment


class IndexInAlignment:
    # Maps original index into alined
    # Note: indices are zero indexed!
    def __init__(self, seq):
        self.indices = [i for i, c in enumerate(seq.seq) if c != '-']

    def __getitem__(self, i):
        return self.indices[i]

    def __len__(self):
        return len(self.indices)


class AlignmentMapIndices:
    def __init__(self, alignment=None, filename=None):
        assert alignment or filename
        self.alignment = alignment or read_alignment(filename)
        self.seq_indices = dict((seq.name, IndexInAlignment(seq)) for seq in self.alignment)

    def create_raxml_partitions(self, orig_parts, partitions_filename):
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
            last_start = 1

            for i in range(num_genes):                       # For each gene
                for seq_ident, parts in orig_parts.items():  # Check between all sequences
                    aligned = self.seq_indices[seq_ident]    # Is there start/end pair aligned on neighbouring positions
                    end, gene = parts[i]
                    if end >= len(aligned):
                        output.write(f'DNA, {gene} = {last_start}-{len(aligned)}\n')
                        assert i == (num_genes - 1), (i, num_genes)
                        break
                    p = aligned[end - 1]  # In current part
                    n = aligned[end]      # In next part
                    if n == p + 1:
                        output.write(f'DNA, {gene} = {last_start}-{p}\n')
                        last_start = n + 1  # Indices are zero indexed!
                        break
                else:
                    raise ValueError(parts[i][1])
