from collections import defaultdict
from itertools import groupby
from common_utils.exceptions import ZCItoolsValueError
from ..utils.features import Feature, find_disjunct_genes


class Partitions:
    # Describes genome partition used for phylogenetic calculation
    # Partition is list of tuples (name, end index).
    #  - Indexing is zero based, and end index is used checked with <.
    #  - Last partition has to have end index equal to length of alignment.
    def __init__(self, make_partitions=True, annotations_step=None):
        self.make_partitions = make_partitions
        self.annotations_step = annotations_step

    def _partition_from_part_files(self, align_step):
        assert False, 'Not tested yet!!!'
        # Collect data
        sequences = align_step.all_sequences()  # Of type set
        partitions = dict((seq_ident, p) for seq_ident in sequences if (p := align_step.get_partition_from_file(seq_ident)))
        if (not_in := (sequences - set(partitions.keys()))):
            if len(not_in) == len(sequences):
                raise ZCItoolsValueError('None of aligned sequences has partiotion file!')
            else:
                raise ZCItoolsValueError(f"Sequence(s) without partition file! {', '.join(sorte(not_in))}")

        # Check partitions
        ami = align_step.get_alignment_map_indices()
        num_genes = len(next(iter(partitions.values())))
        # All have to have same number of parts
        assert all(len(p) == num_genes for p in partitions.values()), [(s, len(p)) for s, p in partitions.items()]
        # All have to be indexed
        assert not (not_in := (sequences - set(ami.all_sequence()))), not_in

        # Create list of parts
        alignment_length = ami.alignment_length
        for i in range(num_genes):                       # For each gene
            for seq_ident, parts in partitions.items():  # Check between all sequences
                aligned = ami.seq_indices[seq_ident]     # Is there start/end pair aligned on neighbouring positions
                end, gene = parts[i]
                if end >= alignment_length:
                    partition.append((gene, alignment_length))
                    assert i == (num_genes - 1), (i, num_genes)
                    break
                p = aligned[end - 1]  # Index in current part
                n = aligned[end]      # Index in next part
                if n == p + 1:        # If positions are neighbouring than we have boundary
                    partition.append((gene, n))
                    break
            else:
                raise ValueError(parts[i][1])
        return partition

    def _partition_from_part_annotations(self, align_step):
        if not (annotations := self.annotations_step):
            raise ZCItoolsValueError('Annotations step is not specified!')

        sequences = align_step.all_sequences()  # Of type set
        if (not_in := (sequences - annotations.all_sequences())):
            if len(not_in) == len(sequences):
                raise ZCItoolsValueError('None of aligned sequences are present in annotations step!')
            else:
                raise ZCItoolsValueError(f"Sequence(s) not present in annotations step! {', '.join(sorte(not_in))}")

        #
        ami = align_step.get_alignment_map_indices(with_reverse=True)
        seq_recs = dict((seq_ident, annotations.get_sequence_record(seq_ident)) for seq_ident in sequences)
        data = [(seq_ident, find_disjunct_genes_fix_name(seq_rec), ami.get_alignment(seq_ident))
                for seq_ident, seq_rec in seq_recs.items()]
        # Take (some) sequence with the most genes
        max_ident, max_seq, max_align = max(data, key=lambda x: len(x[1]))
        partition = []
        name_counter = defaultdict(int)

        def _add_p(name, end_idx):
            partition.append((f'{name}_{name_counter[name]}', end_idx))
            name_counter[name] += 1

        for max_gene in max_seq:
            # Naci mu range u alignmentu
            half = max_align[(max_gene.real_start + max_gene.real_end) // 2]  # Position in alignment
            ends = []
            not_in_seqs = []
            wrong_gene = False
            for seq_ident, genes, align in data:
                seq_idx = align.reverse(half)
                for g in genes:
                    s, e = g.ends()
                    if s <= seq_idx < e:
                        if g.name == max_gene.name:
                            ends.append((align[s], align[e - 1] + 1))
                        else:
                            wrong_gene = True
                        break
                else:
                    # Location not in any of sequence genes
                    not_in_seqs.append((seq_ident, align))
                if wrong_gene:
                    break

            if wrong_gene:
                continue

            # Consolidate data
            start = min(s for s, _ in ends)  # ? To think about
            end = max(e for _, e in ends)

            if not_in_seqs:
                # If alignment is 'appropriate' than maybe annotation is missing
                c_length = 0.8 * len(max_gene)  # Length of not gene
                if any((align.reverse_length(start, end) < c_length) for _, align in not_in_seqs):
                    # print(max_gene.name, len(max_gene), c_length, [(s, align.reverse_length(start, end)) for s, align in not_in_seqs])
                    continue

            if (not partition and start > 0) or (partition and partition[-1][-1] < start):
                _add_p('nc', start)
            _add_p(max_gene.name[0], end)

        if partition and partition[-1][-1] < ami.alignment_length:
            _add_p('nc', ami.alignment_length)

        return partition

    def create_partitions(self, align_step):
        if (st := align_step.get_sequence_type()) == 'genes':
            return self._partition_from_part_files(align_step)
        if st == 'whole':
            return self._partition_from_part_annotations(align_step)
        assert st == 'gene', f'Wrong alignment sequence type: {st}'

    def create_raxml_partitions(self, align_step, partitions_filename):
        if self.make_partitions and (partition := self.create_partitions(align_step)):
            with open(partitions_filename, 'w') as output:
                output.write('\n'.join(f'DNA, {g} = {f}-{t}' for g, f, t in self._from_to_1(align_step, partition)))
            return partitions_filename

    def create_mrbayes_partitions(self, align_step):
        if self.make_partitions and (partition := self.create_partitions(align_step)):
            p = '\n'.join(f'charset {g} = {f}-{t};' for g, f, t in self._from_to_1(align_step, partition))
            return f'{p}\npartition partition_1 = {len(partition)}: {", ".join(g for g, _ in partition)};'
        return ''

    #
    def _from_to_1(self, align_step, partition):
        last_start = 1
        for name, end in partition:
            yield name, last_start, end
            last_start = end + 1


# In some helper file?
def find_disjunct_genes_fix_name(seq_rec):
    # Gene is identified with a name and a strand
    genes = find_disjunct_genes(seq_rec)
    for g in genes:
        g.name = (g.name, g.strand)
    return genes


def longest_common_subsequence_genes(lists_of_genes):
    # Find shared genes
    shared_genes = set(g.name for g in lists_of_genes[0])
    for gs in lists_of_genes[1:]:
        shared_genes &= set(g.name for g in gs)
    if not shared_genes:
        return

    # Filter genes to shared
    # Note: [g for g in gs if g in shared_genes] is needed since gene can occure on more location
    lists_of_genes = [[g for g in gs if g.name in shared_genes] for gs in lists_of_genes]
    return longest_common_subsequence(lists_of_genes, lambda x: x.name)


def longest_common_subsequence(list_of_seqs, extract_str):
    if not list_of_seqs or any(not l for l in list_of_seqs):
        return []

    if extract_str is None:
        extract_str = lambda x: x

    num_seqs = len(list_of_seqs)
    seq_lens = [len(l) for l in list_of_seqs]
    idx_to_check = [0] * num_seqs
    lsc = [[] for _ in range(list_of_seqs)]

    return _lsc_rec(list_of_seqs, seq_lens, 0, idx_to_check, lsc)


def _lsc_rec(list_of_seqs, seq_lens, found_length, idx_to_check, lsc):
    # Check is it possible to find subsequence of length longer than found_length
    max_additional = min(seq_lens[i] - idx_to_check[i] for i in range(len(lsc)))
    if max_additional <= 0 or (len(lsc[0]) + max_additional) <= found_length:
        return

    # Find next chars in sequences
    next_pos = defaultdict(list)
    for l_idx, (seq, idx) in enumerate(zip(list_of_seqs, idx_to_check)):
        next_pos[extract_str(seq[idx])].append(l_idx)

    # If all characters are the same, than add it in lsc
    if len(next_pos) == 1:
        # Found common char
        next_idx_to_check = [i + 1 for i in idx_to_check]
        next_lsc = [ls + [list_of_seqs[idx]] for ls, idx in zip(lsc, idx_to_check)]
        return _lsc_rec(list_of_seqs, seq_lens, found_length, next_idx_to_check, next_lsc)

    # Try all chars, prefer ones with more occurences
    best_lsc = lsc
    for _, possibilities in groupby(next_pos.values(), key=len):
        for use_idxs in sorted(possibilities, key=lambda idxs: min(seq_lens[i] - idx_to_check[i] for i in idxs)):
            next_idx_to_check = list(idx_to_check)
            for idx in use_idxs:
                next_idx_to_check[idx] += 1
            next_lsc = lsc

            if new_lsc := _lsc_rec(list_of_seqs, seq_lens, found_length, next_idx_to_check, lsc):
                best_lsc = new_lsc
                assert found_length < len(best_lsc[0]), (found_length < len(best_lsc[0]))
                found_length = len(best_lsc[0])

    return best_lsc
