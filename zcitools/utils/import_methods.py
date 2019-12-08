from functools import wraps

"""
Method for importing non standard python libraries.
Idea of this methods is to inform user how to import needed library if it is not installed.

Methods receive *args as arguments and return:
 - library if args is empty,
 - list of library attributes if args are given.
"""


# Wrapper for import protocol
def _import_method(error_msg):
    def wrap(func):
        def wrapped_f(*args):
            try:
                lib = func()
            except ImportError:
                print()
                print('=' * 80)
                print(error_msg)
                print('=' * 80)
                print()
                raise

            if args:
                if len(args) == 1:
                    return getattr(lib, args[0])  # No unpacking!
                return [getattr(lib, a) for a in args]
            return lib

        return wrapped_f
    return wrap


# Biopython
_missing_bio = """
Biopython library (https://biopython.org/) is missing.

For installation instruction check web page:
https://biopython.org/wiki/Download

Short: pip install biopython
"""


@_import_method(_missing_bio)
def import_bio_seq_io():
    from Bio import SeqIO
    return SeqIO


@_import_method(_missing_bio)
def import_bio_align_io():
    from Bio import AlignIO
    return AlignIO


@_import_method(_missing_bio)
def import_bio_alphabet():
    from Bio import Alphabet
    return Alphabet


@_import_method(_missing_bio)
def import_bio_entrez():
    from Bio import Entrez
    return Entrez


_missing_cai = """
Codon Adaptation Index (CAI) is missing.

For installation instruction check web page:
https://github.com/Benjamin-Lee/CodonAdaptationIndex

Short: pip install cai
"""


@_import_method(_missing_cai)
def import_CAI():
    import CAI
    return CAI


_missing_pygraphviz = """
PyGraphviz library is missing.

For installation instruction check web page:
https://pygraphviz.github.io/documentation/latest/install.html

Short: pip install pygraphviz
"""


@_import_method(_missing_pygraphviz)
def import_pygraphviz():
    import pygraphviz
    return pygraphviz
