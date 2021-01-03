from functools import wraps
from .import_methods import import_ete3_Tree
from common_utils.cache import cache
from common_utils.exceptions import ZCItoolsValueError


def _distance(func):
    @wraps(func)
    def func_wrapper(self, t2, rooted):
        assert isinstance(t2, PhylogeneticTree), type(t2)
        if rooted:
            et1 = self.rooted_tree()
            et2 = t2.rooted_tree()
        else:
            et1 = self.unrooted_tree()
            et2 = t2.unrooted_tree()
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
        #
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
        assert len(orient_node) == 1
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

    # Distances
    @_distance
    def distance_robinson_foulds(self, et1, et2, rooted):
        return et1.compare(et2, unrooted=not rooted)

    @_distance
    def distance_branche_score(self, et1, et2, rooted):
        if rooted:
            assert False, 'ToDo...'
        else:
            return sum(abs(a - b) for a, b in self._branch_lengths(et1, et2))

    @_distance
    def distance_kendall_colijn(self, et1, et2, rooted):
        pass
