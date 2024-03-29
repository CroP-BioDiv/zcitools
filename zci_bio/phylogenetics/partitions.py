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
        # Collect data
        sequences = align_step.all_sequences()  # Of type set
        partitions = dict((seq_ident, p) for seq_ident in sequences
                          if (p := align_step.get_partition_from_file(seq_ident)))
        if (not_in := (sequences - set(partitions.keys()))):
            if len(not_in) == len(sequences):
                raise ZCItoolsValueError('None of aligned sequences has partiotion file!')
            else:
                raise ZCItoolsValueError(f"Sequence(s) without partition file! {', '.join(sorte(not_in))}")

        # Check partitions
        ami = align_step.get_alignment_map_indices()
        some_parts = next(iter(partitions.values()))
        num_genes = len(some_parts)
        # All have to have same number of parts
        assert all(len(p) == num_genes for p in partitions.values()), [(s, len(p)) for s, p in partitions.items()]
        # All have to be indexed
        assert not (not_in := (sequences - set(ami.all_sequences()))), not_in

        # Create list of parts
        alignment_length = ami.alignment_length
        partition = []
        last_end = 0
        for i in range(num_genes - 1):                   # For each gene, except the last
            for seq_ident, parts in partitions.items():  # Check between all sequences
                aligned = ami.seq_indices[seq_ident]     # Is there start/end pair aligned on neighbouring positions
                end, gene = parts[i]
                if end >= alignment_length:
                    partition.append((gene, [(last_end, alignment_length)]))
                    assert i == (num_genes - 1), (i, num_genes)
                    break
                p = aligned[end - 1]  # Index in current part
                n = aligned[end]      # Index in next part
                if n == p + 1:        # If positions are neighbouring than we have boundary
                    partition.append((gene, [(last_end, n)]))
                    last_end = n
                    break
            else:
                raise ValueError(parts[i][1])

        # Last partition
        _, gene = some_parts[num_genes - 1]
        partition.append((gene, [(last_end, alignment_length)]))
        return partition

    def _partition_from_part_annotations(self, align_step, min_nc_length=10):
        if not (annotations := self.annotations_step):
            raise ZCItoolsValueError('Annotations step is not specified!')

        sequences = align_step.all_sequences()  # Of type set
        if (not_in := (sequences - annotations.all_sequences())):
            if len(not_in) == len(sequences):
                raise ZCItoolsValueError('None of aligned sequences are present in annotations step!')
            else:
                raise ZCItoolsValueError(f"Sequence(s) not present in annotations step! {', '.join(sorted(not_in))}")

        #
        ami = align_step.get_alignment_map_indices(with_reverse=True)
        if not_in := [s for s in ami.all_sequences() if s not in sequences]:
            raise ZCItoolsValueError(f'Not all aligned sequences have index: {", ".join(sorted(not_in))}')

        seq_recs = dict((seq_ident, annotations.get_sequence_record(seq_ident)) for seq_ident in sequences)
        data = [(seq_ident, find_disjunct_genes_fix_name(seq_rec), ami.get_alignment(seq_ident))
                for seq_ident, seq_rec in seq_recs.items()]
        # Take (some) sequence with the most genes
        max_ident, max_seq, max_align = max(data, key=lambda x: len(x[1]))
        nc_parts = []    # Non gene partition. List of tuple (start_idx, end_idx).
        gene_parts = []  # Gene partiotions. List of tuples (name, start_idx, end_idx)
        name_counter = defaultdict(int)

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

            # Check for minimal non coding length
            last_end = gene_parts[-1][-1] if gene_parts else 0
            if start < last_end + min_nc_length:
                start = last_end
            else:
                nc_parts.append((last_end, start))

            name = max_gene.name[0]
            gene_parts.append((f'{name}_{name_counter[name]}', start, end))
            name_counter[name] += 1

        if not gene_parts:  # No partitions at all
            return

        last_part = gene_parts[-1]
        last_end = last_part[-1]
        if ami.alignment_length < last_end + min_nc_length:
            gene_parts[-1] = (last_part[0], last_part[1], ami.alignment_length)
        else:
            nc_parts.append((last_end, ami.alignment_length))

        return [(n, [(s, e)])for n, s, e in gene_parts] + [('nc', nc_parts)]

    def create_partitions(self, align_step, from_one=True, min_nc_length=10):
        # Returns list of tuples (name, list of pairs of indices [from index, to index])
        if (st := align_step.get_sequence_type()) == 'genes':
            partition = self._partition_from_part_files(align_step)
        elif st == 'whole':
            partition = self._partition_from_part_annotations(align_step, min_nc_length=min_nc_length)
        else:
            assert st == 'gene', f'Wrong alignment sequence type: {st}'
            return
        # Check partition
        assert all(isinstance(x[0], str) for x in partition), [x for x in partition if not isinstance(x[0], str)]
        assert all(isinstance(x[1], (list, tuple)) and bool(x) for x in partition), [x for x in partition if not (isinstance(x[1], (list, tuple)) and bool(x))]
        assert all(all(isinstance(y, (list, tuple)) and len(y) == 2 for y in x[1]) for x in partition), partition

        # Reindex
        if from_one:
            partition = [[x, [(s + 1, e) for s, e in indices]] for x, indices in partition]
        else:
            partition = [[x, [(s, e - 1) for s, e in indices]] for x, indices in partition]
        return partition

    def create_raxml_partitions(self, align_step, partitions_filename):
        if self.make_partitions and (partition := self.create_partitions(align_step)):
            with open(partitions_filename, 'w') as output:
                # Note: RAxML doesn't like partitions of length 1 like: "DNA, gene = n-n"
                output.write('\n'.join(f'DNA, {g} = {", ".join(f"{s}-{e}" for s, e in idxs)}' for g, idxs in partition))
            return partitions_filename

    def create_mrbayes_partitions(self, align_step, nexus_file, ident='    '):
        if self.make_partitions and (partition := self.create_partitions(align_step, min_nc_length=10)):
            def _fix_n(n):
                n = n.replace('.', '_')
                n = n.replace('-', '_')
                return n

            partition = [(_fix_n(g), idxs) for g, idxs in partition]
            names = [g for g, _ in partition]

            # Note: MrBayes doesn't like partition with more parts
            # It starts calculation but can hangs after some time :-/
            # Because of that alignment is repartioned.
            if True:  # all(len(idxs) == 1 for g, idxs in partition):
                p = '\n'.join(f'{ident}charset {g} = {" ".join(f"{s}-{e}" for s, e in idxs)};' for g, idxs in partition)
            else:
                # Open nexus_file
                from ..utils.import_methods import import_bio_align_io
                AlignIO = import_bio_align_io()
                alignment = AlignIO.read(nexus_file, 'nexus')

                s, e = partition[0][1][0]
                new_a = alignment[:, s:e]
                for s, e in partition[0][1][1:]:
                    new_a += alignment[:, s:e]
                lengths = [sum((e - s) for s, e in partition[0][1])]
                # Repartition
                for g, idxs in partition[1:]:
                    for s, e in idxs:
                        new_a += alignment[:, s:e]
                    lengths.append(sum((e - s) for s, e in idxs))

                # Save
                with open(nexus_file, "w") as handle:
                    # Check Bio/AlignIO/NexusIO.py, method _classify_mol_type_for_nexus()
                    # ToDo: how to know what type to set?
                    for x in new_a:
                        if 'molecule_type' not in x.annotations:
                            x.annotations['molecule_type'] = 'DNA'
                    AlignIO.write(new_a, handle, 'nexus')

                #
                p = ''
                last_end = 1
                for g, _l in zip(names, lengths):
                    p += f'{ident}charset {g} = {last_end}-{last_end + _l - 1};\n'
                    last_end += _l

            return f"""
{p}
{ident}partition partition_1 = {len(names)}: {",".join(names)};
{ident}set partition = partition_1;
"""
        return ''


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
