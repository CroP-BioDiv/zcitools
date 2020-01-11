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
def import_bio_alphabet():
    from Bio import Alphabet
    return Alphabet


@import_method(_missing_bio)
def import_bio_entrez():
    from Bio import Entrez
    # Set email from settings file
    from common_utils.file_utils import get_settings
    email = get_settings()['email']
    if email:
        Entrez.email = email
    return Entrez


_missing_bcbio = """
GFF parsing library (https://biopython.org/wiki/GFF_Parsing) is missing.

For installation instruction check web page:
https://github.com/chapmanb/bcbb/tree/master/gff

Short: pip install bcbio-gff
"""


@import_method(_missing_bcbio)
def import_bcbio_gff():
    from BCBio import GFF
    return GFF


_missing_cai = """
Codon Adaptation Index (CAI) is missing.

For installation instruction check web page:
https://github.com/Benjamin-Lee/CodonAdaptationIndex

Short: pip install cai
"""


@import_method(_missing_cai)
def import_CAI():
    import CAI
    return CAI
