from common_utils.import_method import import_method

# Biopython
_missing_bio = """
Biopython library (https://biopython.org/) is missing.

For installation instruction check web page:
https://biopython.org/wiki/Download

Short: pip install biopython
"""


@import_method(_missing_bio)
def import_bio_seq_io():
    from Bio import SeqIO
    return SeqIO


@import_method(_missing_bio)
def import_bio_align_io():
    from Bio import AlignIO
    return AlignIO


@import_method(_missing_bio)
def import_bio_phylo():
    from Bio import Phylo
    return Phylo


@import_method(_missing_bio)
def import_bio_entrez():
    from Bio import Entrez
    # Set email from settings file
    from common_utils.file_utils import get_settings
    email = get_settings()['email']
    if email:
        Entrez.email = email
    return Entrez


@import_method("""
ETE Toolkit library (http://etetoolkit.org/docs/latest/index.html) is missing.

For installation instruction check web page:
http://etetoolkit.org/download/

Short: pip install ete3
""")
def import_ete3_NCBITaxa():
    from ete3 import NCBITaxa
    return NCBITaxa


@import_method("""
GFF parsing library (https://biopython.org/wiki/GFF_Parsing) is missing.

For installation instruction check web page:
https://github.com/chapmanb/bcbb/tree/master/gff

Short: pip install bcbio-gff
""")
def import_bcbio_gff():
    from BCBio import GFF
    return GFF


@import_method("""
Codon Adaptation Index (CAI) is missing.

For installation instruction check web page:
https://github.com/Benjamin-Lee/CodonAdaptationIndex

Short: pip install cai
""")
def import_CAI():
    import CAI
    return CAI
