from common_utils.terminal_layout import StringColumns
# from Bio.SeqFeature import FeatureLocation, CompoundLocation


def _irb_start_end(irb):
    loc = irb.location
    if loc.__class__.__name__ == 'CompoundLocation':
        assert len(loc.parts) == 2, (len(loc.parts), irb)
        for l in loc.parts:
            if l.start != 0:
                return l.start, l.end
        assert False, irb
    return loc.start, loc.end


def chloroplast_annotation(annotations, num_genes=1, feature_type='gene', features=None, sequences=None):
    sequences = (annotations._sequences & set(sequences)) if sequences else annotations._sequences
    rows = []
    have_first_col = False
    for seq_ident in sorted(sequences):
        seq = annotations.get_sequence_record(seq_ident)
        rep_regs = [f for f in seq.features if f.type == 'repeat_region']
        if len(rep_regs) != 2:
            print(f"Sequence {seq_ident} doesn't have inverted repeats!")
            continue
        ira, irb = rep_regs
        ira_start, ira_end = ira.location.start, ira.location.end
        irb_start, irb_end = _irb_start_end(irb)
        if len(seq) != irb_end:
            print(f"  warning ({seq_ident}): sequence's IRB doesn't end on sequence end!",
                  len(seq), irb_end)

        locs = [(1, ira_start), (ira_start, ira_end), (ira_end, irb_start)]
        regions = [[] for _ in range(4)]
        borders = [[] for _ in range(4)]
        features_s = sorted((f for f in seq.features if f.type == feature_type and f.location),
                            key=lambda x: x.location.start)
        for f in features_s:
            ls = f.location.start
            le = f.location.end
            if le < ls:
                borders[0].append(f)
            else:
                for bi, (i1, i2) in enumerate(locs):
                    if le <= i2:
                        regions[bi].append(f)
                        break
                    elif ls < i2:
                        borders[bi + 1].append(f)
                        break
                else:
                    regions[-1].append(f)
        # Checks
        if any(len(b) > 1 for b in borders):
            print(f"  warning ({seq_ident}): more features on border!", [len(b) for b in borders])
        if borders[0]:
            have_first_col = True
        #
        header = ['|', f'LSC ({len(regions[0])})', '  |', f'IRa ({len(regions[1])})', '  |',
                  f'SSC ({len(regions[2])})', '  |', f'IRb  ({len(regions[3])})']
        indices = ['', '', str(ira_start), '', str(ira_end), '', str(irb_start), '']
        row = []

        for i in range(4):
            row.append(borders[i][0].qualifiers['gene'][0] if borders[i] else '')
            rs = regions[i]
            if len(rs) <= num_genes + 1:  # 2 * num_genes:
                # All
                row.append(','.join(r.qualifiers['gene'][0] for r in rs))
            else:
                row.append(','.join(r.qualifiers['gene'][0] for r in rs[:num_genes]) + ' - ' +
                           ','.join(r.qualifiers['gene'][0] for r in rs[-1:]))

        #
        rows.append([seq_ident, str(len(seq)), f'({len(features_s)})'] + [''] * 5)
        rows.extend([indices, header, row])
    #
    if not have_first_col:
        rows = [r[1:] if i % 3 else r[:-1] for i, r in enumerate(rows)]
    print(StringColumns(rows))
