import itertools
from functools import wraps
from .import_methods import import_ete3_Tree
from common_utils.cache import cache
from common_utils.exceptions import ZCItoolsValueError


def _check_same_nodes(et1, et2):
    # Check
    n1 = set(et1.iter_leaf_names())
    n2 = set(et2.iter_leaf_names())
    not_in_2 = n1 - n2
    not_in_1 = n2 - n1
    if not_in_2 or not_in_1:
        if not_in_2:
            print("Nodes missing in the second tree:", ', '.join(sorted(not_in_2)))
        if not_in_1:
            print("Nodes missing in the first tree:", ', '.join(sorted(not_in_1)))
        raise ZCItoolsValueError("Trees don't have same set of nodes.")


def _distance(func):
    @wraps(func)
    def func_wrapper(self, t2, rooted=False):
        assert isinstance(t2, PhylogeneticTree), type(t2)
        if rooted:
            et1 = self.rooted_tree()
            et2 = t2.rooted_tree()
        else:
            et1 = self.unrooted_tree()
            et2 = t2.unrooted_tree()
        _check_same_nodes(et1, et2)
        return func(self, et1, et2, rooted)
    return func_wrapper


class PhylogeneticTree:
    # Tree based on Ete3
    def __init__(self, newick_filename, outgroup, rename_nodes=None):
        self.newick_filename = newick_filename
        self.outgroup = outgroup
        self.rename_nodes = rename_nodes  # Callable that returns new name

    def _load_tree(self):
        tree = import_ete3_Tree()(self.newick_filename)
        if self.rename_nodes:
            for node in tree.traverse():
                node.name = self.rename_nodes(node.name)
        return tree

    @cache
    def unrooted_tree(self):
        return self._load_tree()

    @cache
    def rooted_tree(self):
        tree = self._load_tree()
        tree.set_outgroup(self.outgroup)
        return tree

    def _branch_splits(self, tree):
        # Return iterable of pairs (node of a branch split, tuple of sorted node names containing outgroup)
        # Note: for implementation check ete3 method TreeNode.iter_edges()
        orient_node = tree.get_leaves_by_name(self.outgroup)
        assert len(orient_node) == 1, (self.outgroup, orient_node, self.newick_filename)
        orient_node = orient_node[0]
        cached_content = tree.get_cached_content()
        all_leaves = cached_content[tree]
        for n, side1 in cached_content.items():
            if orient_node in side1:
                yield (n, tuple(sorted(x.name for x in side1)))
            else:
                yield (n, tuple(sorted(x.name for x in (all_leaves - side1))))

    def _branch_lengths(self, et1, et2):
        d1 = dict((ns, n.dist) for n, ns in self._branch_splits(et1))
        d2 = dict((ns, n.dist) for n, ns in self._branch_splits(et2))
        all_ns = set(d1.keys())
        all_ns.update(d2.keys())
        return ((d1.get(ns, 0), d2.get(ns, 0)) for ns in all_ns)

    def _edge_stats(self, tree):
        dists = [node.dist for node in tree.traverse()]
        sum_length = sum(dists)
        num_edges = len(dists)
        return dict(min_length=min(dists),
                    max_length=max(dists),
                    sum_length=sum_length,
                    average_length=sum_length / num_edges,
                    num_edges=num_edges)

    # Distances
    @_distance
    def distance_robinson_foulds(self, et1, et2, rooted):
        # Keys: rf, max_rf, ref_edges_in_source, source_edges_in_ref, effective_tree_size,
        #       norm_rf, treeko_dist, source_subtrees, common_edges, source_edges, ref_edges
        r = et1.compare(et2, unrooted=not rooted)
        return [r['rf'], r['max_rf']]  # These values are all we need

    @_distance
    def distance_branche_score(self, et1, et2, rooted):
        if rooted:
            assert False, 'ToDo...'
        else:
            stat_1 = self._edge_stats(et1)
            stat_2 = self._edge_stats(et2)
            return dict(bs=sum(abs(a - b) for a, b in self._branch_lengths(et1, et2)),
                        stat_1=stat_1,
                        stat_2=stat_2)

    def distance_kendall_colijn(self, t2, lambda_factor=0.5):
        # Kendall-Colijn is on rooted trees
        _check_same_nodes(self.rooted_tree(), t2.rooted_tree())
        # dλ(Ta, Tb) = ||vλ(Ta) - vλ(Tb)||
        kc1 = self.kendall_colijn_lambda(lambda_factor)
        kc2 = t2.kendall_colijn_lambda(lambda_factor)
        return sum((a - b)**2 for a, b in zip(kc1, kc2))**0.5

    def distance_kendall_colijn_topo(self, t2):
        # Kendall-Colijn with lambda = 0
        ms1, _ = self.kendall_colijn_vectors()
        ms2, _ = t2.kendall_colijn_vectors()
        return sum((a - b)**2 for a, b in zip(ms1, ms2))**0.5

    def kendall_colijn_lambda(self, _l):
        # vλ(T) = (1−λ)m(T) + λM(T)⁠
        ms, Ms = self.kendall_colijn_vectors()
        l_1 = 1 - _l
        return [(l_1 * m + _l * M) for m, M in zip(ms, Ms)]

    @cache
    def kendall_colijn_vectors(self):
        # m(T) = (m_1_2, m_1_3, ..., m_k-1_k, 1, ..., 1)
        # M(T) = (M_1_2, M_1_3, ..., M_k-1_k, p_1, ..., p_k)
        et = self.rooted_tree()
        leaves = sorted(et.iter_leaves(), key=lambda n: n.name)
        ms = []
        Ms = []
        for n1, n2 in itertools.combinations(leaves, 2):
            m = M = 0
            current = et.get_common_ancestor(n1, n2)
            while current != et:
                m += 1
                M += current.dist
                current = current.up
            ms.append(m)
            Ms.append(M)
        ms.extend([1] * len(leaves))
        Ms.extend(n.dist for n in leaves)
        assert len(ms) == len(Ms), (len(ms), len(Ms))
        return ms, Ms
