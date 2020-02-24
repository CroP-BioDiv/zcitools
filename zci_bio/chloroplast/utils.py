from common_utils.exceptions import ZCItoolsValueError
from ..utils.features import Feature, Partition


def find_regions(seq_ident, seq):
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


# def _irb_start_end(irb):
#     loc = irb.location
#     if loc.__class__.__name__ == 'CompoundLocation':
#         assert len(loc.parts) == 2, (len(loc.parts), irb)
#         for l in loc.parts:
#             if l.start != 0:
#                 return l.start, l.end
#         assert False, irb
#     return loc.start, loc.end


# def find_irs_locations(seq_ident, seq):
#     ira, irb = find_irs_features(seq_ident, seq)
#     if not ira:
#         return

#     ira_start, ira_end = ira.location.start, ira.location.end
#     irb_start, irb_end = _irb_start_end(irb)
#     if len(seq) != irb_end:
#         print(f"  warning ({seq_ident}): sequence's IRB doesn't end on sequence end! ({len(seq)}, {irb_end})")
#     return ira_start, ira_end, irb_start, irb_end


# def find_irs_features(seq_ident, seq):
#     rep_regs = [f for f in seq.features if f.type == 'repeat_region']
#     if not rep_regs:
#         return None, None
#     if len(rep_regs) != 2:
#         raise ZCItoolsValueError(f"Sequence {seq_ident} doesn't have inverted repeats!")
#     ira, irb = rep_regs
#     return ira, irb


# def split_features_in_regions(seq_ident, seq, ira_start, ira_end, irb_start, irb_end, feature_type='gene'):
#     locs = [(1, ira_start), (ira_start, ira_end), (ira_end, irb_start)]
#     regions = [[] for _ in range(4)]
#     borders = [[] for _ in range(4)]
#     features_s = sorted((f for f in seq.features if f.type == feature_type and f.location),
#                         key=lambda x: x.location.start)
#     for f in features_s:
#         ls = f.location.start
#         le = f.location.end
#         if le < ls:
#             borders[0].append(f)
#         else:
#             for bi, (i1, i2) in enumerate(locs):
#                 if le <= i2:
#                     regions[bi].append(f)
#                     break
#                 elif ls < i2:
#                     borders[bi + 1].append(f)
#                     break
#             else:
#                 regions[-1].append(f)
#     # Checks
#     if any(len(b) > 1 for b in borders):
#         print(f"  warning ({seq_ident}): more features on border!", [len(b) for b in borders])
#     return features_s, regions, borders
