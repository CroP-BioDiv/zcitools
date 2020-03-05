from common_utils.exceptions import ZCItoolsValueError
from ..utils.features import Feature, Partition


def find_chloroplast_partition(seq_ident, seq):
    # Returns None or Partition object with parts named: lsc, ira, ssc, irb.
    rep_regs = [f for f in seq.features if f.type == 'repeat_region']
    if not rep_regs:
        return None
    if len(rep_regs) != 2:
        raise ZCItoolsValueError(f"Sequence {seq_ident} doesn't have inverted repeats!")
    #
    l_seq = len(seq)
    partition = Partition(
        [Feature(l_seq, name='ira', feature=rep_regs[0]), Feature(l_seq, name='irb', feature=rep_regs[1])],
        fill=True)
    n_parts = partition.not_named_parts()
    assert len(n_parts) == 2, len(n_parts)
    ssc_ind = int(len(n_parts[0]) > len(n_parts[1]))
    n_parts[1 - ssc_ind].name = 'lsc'
    n_parts[ssc_ind].name = 'ssc'
    return partition
