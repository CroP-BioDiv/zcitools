import re
from ..utils.import_methods import import_bio_phylo


def read_tree(filename, format='newick'):
    return import_bio_phylo().read(filename, format)


def read_trees(filename, format='newick'):
    return list(import_bio_phylo().parse(filename, format))
