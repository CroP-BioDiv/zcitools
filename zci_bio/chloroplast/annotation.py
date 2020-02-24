from common_utils.terminal_layout import StringColumns
from .utils import find_regions  # find_irs_locations, split_features_in_regions
from ..utils.features import Feature  # find_irs_locations, split_features_in_regions


def chloroplast_annotation(annotations, num_genes=1, feature_type='gene', features=None, sequences=None):
    sequences = (annotations._sequences & set(sequences)) if sequences else annotations._sequences
    rows = []
    with_first_col = False
    for seq_ident in sorted(sequences):
        print(seq_ident)
        seq = annotations.get_sequence_record(seq_ident)
        partition = find_regions(seq_ident, seq)
        if not partition:
            continue

        ira = partition.get_part_by_name('ira')
        irb = partition.get_part_by_name('irb')
        l_seq = len(seq)
        in_parts = partition.put_features_in_parts(
            Feature(l_seq, feature=f) for f in seq.features if f.type == feature_type and f.location)
        with_first_col |= ('irb-lsc' in in_parts)
        #
        header = ['|', f'LSC ({len(in_parts.get("lsc", []))})', '  |', f'IRa ({len(in_parts.get("ira", []))})', '  |',
                  f'SSC ({len(in_parts.get("ssc", []))})', '  |', f'IRb  ({len(in_parts.get("irb", []))})']
        indices = ['', '', str(ira.real_start), '', str(ira.real_end), '', str(irb.real_start), '']
        row = []

        p_names = list(partition.get_part_names())
        for prev, p_name in zip(p_names[-1:] + p_names[:-1], p_names):
            prev_border = in_parts.get(f'{prev}-{p_name}')
            row.append(prev_border[0].name if prev_border else '')
            rs = in_parts.get(p_name, [])
            if len(rs) <= num_genes + 1:  # 2 * num_genes:
                # All
                row.append(','.join(r.name for r in rs))
            else:
                row.append(','.join(r.name for r in rs[:num_genes]) + ' - ' +
                           ','.join(r.name for r in rs[-1:]))

        #
        rows.append([seq_ident, str(len(seq)), f'({sum(len(fs) for fs in in_parts.values())})'] + [''] * 5)
        rows.extend([indices, header, row, [''] * len(row)])

    #
    if not with_first_col:
        rows = [r[1:] if i % 5 else r[:-1] for i, r in enumerate(rows)]
    print(StringColumns(rows))
