from itertools import chain
from common_utils.cache import cache
from zci_bio.utils.ncbi_taxonomy import get_ncbi_taxonomy, order_ranks


class StatByTaxonomy:
    # Stores statistics (counters) for given taxons and group arguments
    _STAT_ATTRS = []
    _EXCEL_COLUMNS = []

    def __init__(self, taxids, ranks=None, names=None):
        self.taxid_2_stat = dict()  # Taxid -> dict (parent, nums)

        # Find all taxids, species and parent clades
        ncbi_taxonomy = get_ncbi_taxonomy()
        all_taxids = set(taxids)
        gl = ncbi_taxonomy.get_lineage
        self.taxid_2_lineage = dict((t, gl(t)) for t in all_taxids)
        all_taxids.update(chain.from_iterable(self.taxid_2_lineage.values()))
        self.taxid_2_rank = ncbi_taxonomy.get_rank(all_taxids)
        self.taxid_2_name = ncbi_taxonomy.get_taxid_translator(all_taxids)

        # Extract only taxid of interest
        self.group_taxids = set()
        if ranks:
            self.group_taxids.update(t for t, r in self.taxid_2_rank.items() if r in ranks)
        if names:
            for name in names:
                if n_taxid := ncbi_taxonomy.get_exact_name_translation(name):
                    self.group_taxids.add(n_taxid)

    def stat_attrs(self):
        return self._STAT_ATTRS

    def excel_columns(self):
        return self._EXCEL_COLUMNS

    @cache
    def name_2_taxid(self):
        return dict((n, t) for t, n in self.taxid_2_name.items())

    def add(self, taxid, *args, **kwargs):
        parents = [t for t in self.taxid_2_lineage[taxid] if t in self.group_taxids]
        if not parents:
            print(f'Warning: stat for taxid {taxid} is not set, because there is no parent to add to!')
            return

        # Add stat to parent nodes
        for taxid, parent in zip(parents, [None] + parents[:-1]):
            if node := self.taxid_2_stat.get(taxid):
                assert node['parent'] == parent, (taxid, node['parent'], parent)
            else:
                node = dict(taxid=taxid, parent=parent, num=0)
                self._init_new_node(node)
                self.taxid_2_stat[taxid] = node
            #
            node['num'] += 1
            self._add(node, *args, **kwargs)

    def _init_new_node(self, node):
        node.update((a, 0) for a in self.stat_attrs())

    def _sum(self, stats):
        # Note: contains num_groups attr!
        assert stats
        ss = dict(taxid=None, parent=stats[0]['parent'], num=0, num_groups=len(stats))
        self._init_new_node(ss)
        for a in ['num'] + self.stat_attrs():
            for s in stats:
                ss[a] += s[a]
        return ss

    #
    def get(self, parent):
        # ToDo: more efficient
        return [(s, self.get(s['taxid'])) for s in self.taxid_2_stat.values() if s['parent'] == parent]

    def _depth(self, stats):
        return (max(self._depth(sub_s) for _, sub_s in stats) + 1) if stats else 0

    def to_table(self, minimum_sequences):
        stats = self.get(None)
        max_depth = self._depth(stats)
        return max_depth, self._to_table(stats, minimum_sequences, max_depth, 0)

    def _to_table(self, stats, minimum_sequences, max_depth, depth):
        # Collaps groups with sequences less than given minimum
        # Note: Single leaf is not collapsed. Single inner node is collapsed.
        if minimum_sequences and \
           (rest := [s for s, _ in stats if s['num'] < minimum_sequences]) \
           and (len(rest) > 1 or depth < max_depth - 1):
            stats = sorted(((s, sub_s) for s, sub_s in stats if s['num'] >= minimum_sequences),
                           key=lambda x: self.taxid_2_name[x[0]['taxid']])
            stats.append((self._sum(rest), []))
        else:
            stats = sorted(stats, key=lambda x: self.taxid_2_name[x[0]['taxid']])

        rows = []
        for s, sub_s in stats:
            rows.append(self._to_row(s, max_depth, depth))
            if sub_s:
                rows += self._to_table(sub_s, minimum_sequences, max_depth, depth + 1)

        # Add sum row
        if depth == 0 and len(stats) > 1:
            rows.append(self._to_sum_row(self._sum([s for s, _ in stats]), max_depth))

        #
        return rows

    def _to_row_label(self, s, max_depth, depth, label):
        row = [''] * max_depth + [s['num']] + [s[a] for a in self.stat_attrs()]
        row[depth] = label
        return row

    def _to_row(self, s, max_depth, depth):
        return self._to_row_label(
            s, max_depth, depth, (self.taxid_2_name[s['taxid']] if s['taxid'] else f"... ({s['num_groups']})"))

    def _to_sum_row(self, s, max_depth):
        return self._to_row_label(s, max_depth, 0, 'all')

    def export_excel(self, filename, minimum_sequences):
        from common_utils.value_data_types import rows_2_excel
        md, rows = self.to_table(minimum_sequences)
        columns = [''] * md + ['Num sequences'] + self.excel_columns()
        rows_2_excel(filename, columns, rows, index=False)

    def export_sheet(self, sheet_name, minimum_sequences):
        md, rows = self.to_table(minimum_sequences)
        columns = [''] * md + ['Num sequences'] + self.excel_columns()
        return (sheet_name, columns, rows)

    #
    def get_node_from_name(self, taxa_name):
        if not (taxid := self.name_2_taxid().get(taxa_name)):
            print(f'No taxon with name {taxa_name}!')
            return
        if taxid not in self.group_taxids:
            print(f'Taxon {taxa_name} is not used for grouping!')
            return
        if not (node := self.taxid_2_stat.get(taxid)):
            print(f'Taxon {taxa_name} has no stat data!')
            return
        return node


#
class _Node:
    def __init__(self, is_leaf, taxid, parent_node, rank, depth, name):
        self.is_leaf = is_leaf
        self.taxid = taxid
        self.parent_node = parent_node
        self.rank = rank
        self.depth = depth
        self.name = name
        self.max_depth = depth  # In subtree
        if is_leaf:
            self.objects = []
        else:
            self.children = []  # List of nodes

    def add_child(self, node):
        self.children.append(node)
        self.set_depth(node.depth)

    def set_depth(self, depth):
        if depth > self.max_depth:
            self.max_depth = depth
            if self.parent_node:
                self.parent_node.set_depth(depth)

    @staticmethod
    def dummy_node(is_leaf=False, depth=0, name='-'):
        return _Node(is_leaf, 0, None, None, depth, name)


_DUMMY_NODE_INNER = _Node.dummy_node(name='-')
_DUMMY_NODE_INNER_S = _Node.dummy_node(name='')


class GroupByTaxonomy:
    # Stores data (list of objects) for given taxons and group arguments
    # Used for grouping rows or statistics on data grouped by taxonomic structure.
    # Note: more general than class StatByTaxonomy. Previous class is here for backward compatibility
    def __init__(self, objects, ranks=None, names=None, taxid_method=None, taxid_attr=None):
        if taxid_method:
            object_taxids = [taxid_method(o) for o in objects]
        elif taxid_attr:
            object_taxids = [o[taxid_attr] for o in objects]
        else:
            assert False, '?'
            # taxid_method = lambda o: o[taxid_attr]

        # Find all taxids, species and parent clades
        all_taxids = set(object_taxids)
        ncbi_taxonomy = get_ncbi_taxonomy()
        gl = ncbi_taxonomy.get_lineage
        taxid_2_lineage = dict((t, gl(t)) for t in all_taxids)
        all_taxids.update(chain.from_iterable(taxid_2_lineage.values()))
        taxid_2_rank = ncbi_taxonomy.get_rank(all_taxids)
        self.taxid_2_name = ncbi_taxonomy.get_taxid_translator(all_taxids)

        # Extract only taxid of interest
        group_taxids = set()
        self.ranks = order_ranks(ranks) if ranks else None
        if ranks:
            group_taxids.update(t for t, r in taxid_2_rank.items() if r in ranks)
        if names:
            for name in names:
                if n_taxid := ncbi_taxonomy.get_exact_name_translation(name):
                    group_taxids.add(n_taxid)

        # Set objects
        self.nodes = dict()
        for taxid, obj in zip(object_taxids, objects):
            parent_taxids = [t for t in taxid_2_lineage[taxid] if t in group_taxids]
            if not parent_taxids:
                print(f'Warning: stat for taxid {taxid} is not set, because there is no parent to add to!')
                continue

            # Add nodes
            # Note: node creation is from higher to lower taxonomy
            for depth, (taxid, parent_taxid) in enumerate(zip(parent_taxids, [None] + parent_taxids[:-1])):
                node = self._add_node(taxid == parent_taxids[-1], taxid, parent_taxid, taxid_2_rank[taxid], depth)
            node.objects.append(obj)  # Last node is a leaf
        #
        # Finish!!!!

    def _add_node(self, is_leaf, taxid, parent_taxid, rank, depth):
        if node := self.nodes.get(taxid):
            assert node.is_leaf == is_leaf, (node.is_leaf, is_leaf, taxid, parent_taxid, depth)
            assert node.taxid == taxid, (node.taxid, is_leaf, taxid, parent_taxid, depth)
            assert node.depth == depth, (node.depth, is_leaf, taxid, parent_taxid, depth)
            return node
        parent = self.nodes[parent_taxid] if parent_taxid is not None else None
        self.nodes[taxid] = node = _Node(is_leaf, taxid, parent, rank, depth, self.taxid_2_name[taxid])
        if parent:
            parent.add_child(node)
        return node

    def grouped_columns(self):
        # Assumes that names are higher nodes and ranks are lower nodes.
        max_depth = max(n.depth for n in self.nodes.values())
        num_names = max_depth + 1 - len(self.ranks)
        assert num_names >= 0, (max_depth, len(self.ranks))
        if num_names == 0:
            return self.ranks
        if num_names == 1:
            return ['Clade'] + self.ranks
        return [f'Clade_{n}' for n in range(1, num_names + 1)] + self.ranks

    def sorted_objects(self, objects_sort=None):
        # Iterate through objects in taxonomical order
        for n in sorted((n for n in self.nodes.values() if n.depth == 0), key=lambda n: (-n.max_depth, n.name)):
            yield from self._sorted_objects(n, objects_sort)

    def _sorted_objects(self, node, objects_sort):
        if node.is_leaf:
            yield from (sorted(n.objects, key=objects_sort) if objects_sort else n.objects)
        else:
            for n in sorted(node.children, key=lambda n: n.name):
                yield from self._sorted_objects(n, objects_sort)

    def sorted_nodes_objects(self, objects_sort=None, return_names=False, compact_names=False, lowest_rank=None):
        # Iterate through objects in taxonomical order.
        # Returns pairs (list of nodes, list of objects)
        max_depth = max(n.depth for n in self.nodes.values()) + 1  # +1 since depth 0 means the heighest node
        s_nodes = sorted((n for n in self.nodes.values() if n.depth == 0), key=lambda n: (-n.max_depth, n.name))
        if lowest_rank == -1:
            all_objs = []
            for n in s_nodes:
                for _, objs in self._sorted_nodes_objects(
                        [n], max_depth, objects_sort, return_names, compact_names, None):
                    all_objs.extend(objs)
            yield [], all_objs
        else:
            for n in s_nodes:
                yield from self._sorted_nodes_objects(
                    [n], max_depth, objects_sort, return_names, compact_names, lowest_rank)

    def _sorted_nodes_objects(self, prev_nodes, max_depth, objects_sort, return_names, compact_names, lowest_rank):
        node = prev_nodes[-1]
        if node.is_leaf:
            if return_names:
                prev_nodes = [n.name for n in prev_nodes]
            if (num_dummies := max(0, max_depth - len(prev_nodes))):
                prev_nodes = ['-' if return_names else _DUMMY_NODE_INNER] * num_dummies + prev_nodes
            yield (prev_nodes, (sorted(node.objects, key=objects_sort) if objects_sort else node.objects))
        else:
            if (lowest_rank is not None) and lowest_rank in (node.rank, node.depth):
                all_objs = []
                for n in sorted(node.children, key=lambda n: n.name):
                    for nodes, objs in self._sorted_nodes_objects(prev_nodes + [n], max_depth, objects_sort,
                                                                  return_names, compact_names, lowest_rank):
                        all_objs.extend(objs)
                if return_names:
                    prev_nodes = [n.name for n in prev_nodes]
                yield prev_nodes, all_objs
            else:
                next_prev_names = ([_DUMMY_NODE_INNER_S] * len(prev_nodes)) if compact_names else prev_nodes
                for n in sorted(node.children, key=lambda n: n.name):
                    yield from self._sorted_nodes_objects(prev_nodes + [n], max_depth, objects_sort,
                                                          return_names, compact_names, lowest_rank)
                    prev_nodes = next_prev_names
