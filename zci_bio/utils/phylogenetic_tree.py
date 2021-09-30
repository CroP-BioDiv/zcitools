import itertools
from functools import wraps
from .import_methods import import_ete3_Tree
from common_utils.cache import cache, cache_args
from common_utils.exceptions import ZCItoolsValueError


def _is_leaf_fn(node):
    return node.is_leaf() and bool(node.name)


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
        self._have_same_nodes(t2, int(rooted))
        return func(self, et1, et2, rooted)
    return func_wrapper


def branch_splits(ete_tree):  # ETE3
    # Returns list of tuples (node, set_of_names_1, set_of_names_2)
    cached_content = ete_tree.get_cached_content()
    all_leaves = cached_content[ete_tree]
    num_l = len(all_leaves)
    return [(n, set(x.name for x in side1), set(x.name for x in (all_leaves - side1)))
            for n, side1 in cached_content.items()
            if 0 != len(side1) != num_l]


class PhylogeneticTree:
    # Tree based on Ete3
    def __init__(self, newick_filename, outgroup, rename_nodes=None, newick_format=0):
        self.newick_filename = newick_filename
        self.outgroup = outgroup
        self.num_leaves = [0, 0]        # unrooted, rooted
        self.node_names = [None, None]  # unrooted, rooted
        self.rename_nodes = rename_nodes  # Callable that returns new name
        self.newick_format = newick_format

    def _load_tree(self, tree_idx):
        tree = import_ete3_Tree()(self.newick_filename, format=self.newick_format)
        if self.rename_nodes:
            for node in tree.traverse():
                node.name = self.rename_nodes(node.name)
        #
        self.node_names[tree_idx] = set(tree.iter_leaf_names(is_leaf_fn=_is_leaf_fn))
        self.num_leaves[tree_idx] = len(self.node_names[tree_idx])
        #
        return tree

    @cache
    def unrooted_tree(self):
        return self._load_tree(0)

    @cache
    def rooted_tree(self):
        tree = self._load_tree(1)
        tree.set_outgroup(self.outgroup)
        return tree

    @cache
    def unrooted_num_edges(self):
        return len(self._branch_splits())

    #
    def _have_same_nodes(self, t2, tree_idx):
        assert self.node_names[tree_idx], tree_idx
        assert t2.node_names[tree_idx], tree_idx
        not_in_2 = self.node_names[tree_idx] - t2.node_names[tree_idx]
        not_in_1 = t2.node_names[tree_idx] - self.node_names[tree_idx]
        if not_in_2 or not_in_1:
            if not_in_2:
                print("Nodes missing in the second tree:", ', '.join(sorted(not_in_2)))
            if not_in_1:
                print("Nodes missing in the first tree:", ', '.join(sorted(not_in_1)))
            raise ZCItoolsValueError("Trees don't have same set of nodes.")

    #
    @cache
    def _branch_splits(self):
        tree = self.unrooted_tree()
        orient_node = tree.get_leaves_by_name(self.outgroup)
        assert len(orient_node) == 1, (self.outgroup, orient_node, self.newick_filename)
        orient_node = orient_node[0]
        cached_content = tree.get_cached_content()
        all_leaves = cached_content[tree]
        num_l = len(all_leaves)
        ret = [(n, side1 if orient_node in side1 else (all_leaves - side1)) for n, side1 in cached_content.items()]
        return [(n, tuple(sorted(x.name for x in s))) for n, s in ret if len(s) != num_l]

    def _branch_lengths(self, t2):
        d1 = dict((ns, n.dist) for n, ns in self._branch_splits())
        d2 = dict((ns, n.dist) for n, ns in t2._branch_splits())
        all_ns = set(d1.keys())
        all_ns.update(d2.keys())
        return ((d1.get(ns, 0), d2.get(ns, 0)) for ns in all_ns)

    @cache_args
    def _edge_stats(self, rooted):
        tree = self.rooted_tree() if rooted else self.unrooted_tree()
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
        return dict(rf=r['rf'], max_rf=r['max_rf'])  # These values are all we need

    def distance_branche_score(self, t2, rooted=False):
        assert not rooted, 'For now!'
        return dict(bs=sum(abs(a - b) for a, b in self._branch_lengths(t2)),
                    stat_1=self._edge_stats(rooted),
                    stat_2=t2._edge_stats(rooted))

    def distance_kendall_colijn(self, t2, lambda_factor=0.5):
        # Kendall-Colijn is on rooted trees
        # dλ(Ta, Tb) = ||vλ(Ta) - vλ(Tb)||
        kc1 = self.kendall_colijn_lambda(lambda_factor)
        kc2 = t2.kendall_colijn_lambda(lambda_factor)
        self._have_same_nodes(t2, 1)
        # Returns tuple (distance, num_ms, num_leaves)
        num_l = self.num_leaves[1]  # On rooted
        return dict(l_1=sum(abs(a - b) for a, b in zip(kc1, kc2)),
                    l_2=sum((a - b)**2 for a, b in zip(kc1, kc2))**0.5,
                    num_ms=num_l * (num_l - 1) // 2,
                    num_leaves=num_l)

    def distance_kendall_colijn_topo(self, t2):
        # Kendall-Colijn with lambda = 0
        return self.distance_kendall_colijn(t2, lambda_factor=0)

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
        leaves = sorted(et.iter_leaves(is_leaf_fn=_is_leaf_fn), key=lambda n: n.name)
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
