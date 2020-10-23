from .analyse_step_3 import find_best_irs_by_similar
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
            if res := find_best_irs_by_similar(seq_descs, seq_ident, d):
                # Only LSC is not good for sure.
                # ToDo: it is even possible to extend SSC part, but it is not important (at all)
                ira, irb, took_ident = res
                ssc = d._parts['ssc']
                ira = (ira[0], ssc.real_start)
                irb = (ssc.real_end, irb[1])
                d.set_took_part(ira, irb, took_ident, 'trnH-GUG')


# def _is_more_left(seq_length, seq_start, a, b):
#     # Returns True if index a < index b in circular genome
#     d = a - b
#     if d == 0:
#         return False
#     if d > 0:  # ?
#         return d > (seq_length // 2)  # ?
#     return d > (seq_length // 2)  # ?
