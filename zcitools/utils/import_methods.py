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
                return [getattr(lib, a) for a in args]
            return lib

        return wrapped_f
    return wrap


# Biopython
_missing_bio = """
Biopython library (https://biopython.org/) is missing.

For installation instruction check web page:
https://biopython.org/wiki/Download
"""

# @_import_method(_missing_bio)
# def import_bio():
#     import Bio
#     return Bio


@_import_method(_missing_bio)
def import_bio_seq_io():
    from Bio import SeqIO
    return SeqIO


@_import_method(_missing_bio)
def import_bio_align_io():
    from Bio import AlignIO
    return AlignIO


@_import_method(_missing_bio)
def import_bio_entrez():
    from Bio import Entrez
    return Entrez
