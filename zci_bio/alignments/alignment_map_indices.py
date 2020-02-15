from ..utils.helpers import read_alignment


class IndexInAlignment:
    # Maps original index into alined
    def __init__(self, seq):
        self.indices = [i + 1 for i, c in enumerate(seq.seq) if c != '-']

    def __getitem__(self, i):
        return self.indices[i]


class AlignmentMapIndices:
    def __init__(self, alignment=None, filename=None):
        assert alignment or filename
        self.alignment = alignment or read_alignment(filename)
        self.seq_indices = dict((seq.name, IndexInAlignment(seq)) for seq in self.alignment)

    def create_raxml_partitions(self, orig_parts, filename):
        # orig_parts is iterable of tuples (seq_ident, partitions)
        # Naci sto je sve anotirano, sortirati po redu
        orig_parts = dict((s, k) for s, k in orig_parts)  # seq_ident -> list of tuples (index, description)
        if not orig_parts:
            return

        num_genes = len(next(iter(orig_parts.values())))
        assert all(len(p) == num_genes for p in orig_parts.values()), \
            [(s, len(p)) for s, p in orig_parts.items()]
        assert all(s in self.seq_indices for s in orig_parts.keys()), \
            [s for s in orig_parts.keys() if s not in self.seq_indices]

        #
        with open(filename, 'w') as output:
            last_start = 1

            for i in range(num_genes - 1):                   # For each gene
                for seq_ident, parts in orig_parts.items():  # Check between all sequences
                    aligned = self.seq_indices[seq_ident]    # Is there start/end pair aligned on neighbouring positions
                    end = parts[i][0]
                    p = aligned[end - 1]
                    n = aligned[end]
                    if n == p + 1:
                        output.write(f'DNA, {parts[i][1]} = {last_start}-{p}\n')
                        last_start = n
                        break
                else:
                    raise ValueError(parts[i][1])

            # Add last one
            gene = next(iter(orig_parts.values()))[-1][1]
            output.write(f'DNA, {gene} = {last_start}-{len(self.alignment[0])}\n')
