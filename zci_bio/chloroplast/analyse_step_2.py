from .analyse_step_3 import find_irs_by_similar
from ..utils.ncbi_taxonomy import get_ncbi_taxonomy


# *Problem*
# Sequence is offseted to the left, since IRs were found short on LSC side than LSC is offseted 'left'.
#
# *Possible solution*
# Gene trnH-GUG is (usually) located at LSC start.
# Idea is, if:
#  - parts are normally detected,
#  - position of trnH-GUG is located right of LSC start,
# than it is possible that IRs were detected shorter than they should be (maybe very small error).
# It is worth cheking for more 'reliable' IRs positions.
#
# Note: if IRs were located shorter on SSC side,
# that can cause problem only if IRs are wrongly oriented and fixing reverts them
#
def evaluate_credibility(seq_descs):
    for seq_ident, d in seq_descs.items():
        if d._parts and abs(d.part_trnH_GUG or 0) > 100:  # Contains parts and trnH-GUG located far enough
            if res := find_irs_by_similar(seq_descs, seq_ident, d):
                d.set_took_part(*res, 'trnH-GUG')
