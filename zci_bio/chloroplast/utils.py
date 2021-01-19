from Bio.SeqFeature import FeatureLocation, CompoundLocation
from common_utils.exceptions import ZCItoolsValueError
from ..utils.features import Feature, Partition
from ..utils.helpers import feature_qualifiers_to_desc


def find_chloroplast_irs(seq, check_length=True):
    # Finds the longest pair of inverted repeats
    _ir = ('inverted',)
    rep_regs = [f for f in seq.features
                if f.type == 'repeat_region' and f.qualifiers.get('rpt_type', _ir)[0] == 'inverted']

    # Repair features with location of type <~start..ira.start>
    half_length = len(seq) // 2
    for f in rep_regs:
        if len(f) >= half_length:
            loc = f.location
            f.location = CompoundLocation([FeatureLocation(loc.end, len(seq), strand=1),
                                           FeatureLocation(0, loc.start + 1, strand=1)])

    if len(rep_regs) >= 2:
        max_len = max(map(len, rep_regs)) - 10  # Some tolerance :-)
        max_regs = [f for f in rep_regs if len(f) >= max_len]
        if len(max_regs) != 2 and not check_length:
            # Backup
            max_regs = sorted(rep_regs, key=len)[-2:]
            if len(max_regs[0]) - len(max_regs[1]) > 1000:  # This is too much of difference
                return None
        if len(max_regs) == 2:
            ira, irb = max_regs
            diff_1 = (irb.location.parts[0].start - ira.location.parts[-1].end) % len(seq)
            diff_2 = (ira.location.parts[0].start - irb.location.parts[-1].end) % len(seq)
            return (ira, irb) if diff_1 < diff_2 else (irb, ira)


def irb_start(irb):
    return int(irb.location.parts[0].start)


def find_chloroplast_partition(seq, check_length=True):
    # Returns None or Partition object with parts named: lsc, ira, ssc, irb.
    irs = find_chloroplast_irs(seq, check_length=check_length)
    if irs:
        ira, irb = irs
        return create_chloroplast_partition(len(seq), ira, irb)


def create_chloroplast_partition(l_seq, ira, irb, in_interval=False):
    if in_interval:
        ps = [Feature(l_seq, name='ira', interval=ira), Feature(l_seq, name='irb', interval=irb)]
    else:
        ps = [Feature(l_seq, name='ira', feature=ira), Feature(l_seq, name='irb', feature=irb)]

    partition = Partition(ps, fill=True)
    n_parts = partition.not_named_parts()
    assert len(n_parts) == 2, len(n_parts)
    ssc_ind = int(len(n_parts[0]) > len(n_parts[1]))
    n_parts[1 - ssc_ind].name = 'lsc'
    n_parts[ssc_ind].name = 'ssc'
    return partition


def create_chloroplast_partition_all(l_seq, starts):
    assert len(starts) == 4, starts
    return Partition([Feature(l_seq, name=n, interval=(s, e))
                      for n, s, e in zip(('lsc', 'ira', 'ssc', 'irb'), starts, starts[1:] + starts[:1])])


def find_referent_genome(seq_idents, referent_seq_ident):
    if referent_seq_ident in seq_idents:
        return referent_seq_ident
    refs = [seq_ident for seq_ident in seq_idents if seq_ident.startswith(referent_seq_ident)]
    if not refs:
        raise ZCItoolsValueError(f'No referent genome which name starts with {referent_seq_ident}!')
    elif len(refs) > 1:
        raise ZCItoolsValueError(f'More genomes which name starts with {referent_seq_ident}!')
    return refs[0]


def chloroplast_parts_orientation(seq_rec, partition, genes=None):
    # Check chloroplast sequence part orientation.
    # Default orientation is same as one uses in Fast-Plast. Check:
    #  - source file orientate_plastome_v.2.0.pl
    #    (https://github.com/mrmckain/Fast-Plast/blob/master/bin/orientate_plastome_v.2.0.pl)
    #  - explanation https://github.com/mrmckain/Fast-Plast/issues/22#issuecomment-389302606
    # Consitent with Wikipedia image:
    #  - https://en.wikipedia.org/wiki/Chloroplast_DNA#/media/File:Plastomap_of_Arabidopsis_thaliana.svg

    l_seq = len(seq_rec)
    if genes is None:
        genes = [f for f in seq_rec.features if f.type == 'gene']
    in_parts = partition.put_features_in_parts(Feature(l_seq, feature=f) for f in genes)

    lsc_count = sum(f.feature.strand for f in in_parts.get('lsc', []) if any(x in f.name for x in ('rpl', 'rps')))
    ssc_count = sum(f.feature.strand for f in in_parts.get('ssc', []))
    ira_count = sum(f.feature.strand for f in in_parts.get('ira', []) if 'rrn' in f.name)

    return dict(lsc=(lsc_count <= 0),
                ssc=(ssc_count <= 0),
                ira=(ira_count >= 0))


# -----------------
def cycle_distance_min(a, b, cycle_len):
    d = abs(a - b)
    return d if d < (cycle_len / 2) else (cycle_len - d)


def cycle_distance_lt(a, b, cycle_len):
    # In positive direction. Takes that a <= b in cycle direction.
    ret = (b - a) % cycle_len
    return (ret - cycle_len) if ret > cycle_len // 2 else ret  # Negative is clearer than large number


def rotate_by_offset(seq_rec, offset, keep_offset=None, reverse=False):
    # Returns None if there is no need for rotation or new SeqRecord
    if offset is None:
        return seq_rec.reverse_complement() if reverse else None

    offset %= len(seq_rec.seq)

    # Check if there is need to rotate
    if not offset or (keep_offset and cycle_distance_min(0, offset, len(seq_rec.seq)) <= keep_offset):
        return seq_rec.reverse_complement() if reverse else None

    new_seq = seq_rec[offset:] + seq_rec[:offset]
    return new_seq.reverse_complement() if reverse else new_seq


def rotate_to_offset(seq_rec, parts, keep_offset=None):
    return rotate_by_offset(seq_rec, parts['lsc'].real_start, keep_offset=keep_offset)


# Orientate parts
def orient_chloroplast_parts_by_data(seq_rec, orientation, starts=None, partition=None):
    if not partition:
        if starts:
            partition = create_chloroplast_partition_all(len(seq_rec.seq), starts)
        else:
            raise ZCItoolsValueError(f'No partition data to orient chloroplast parts for sequence {seq_rec.name}!')

    parts = partition.extract(seq_rec)  # dict name -> Seq object
    if 'lsc' in orientation:  # LSC
        parts['lsc'] = parts['lsc'].reverse_complement()
    if 'ssc' in orientation:  # SSC
        parts['ssc'] = parts['ssc'].reverse_complement()
    if 'ira' in orientation:  # IRs
        parts['ira'], parts['irb'] = parts['irb'].reverse_complement(), parts['ira'].reverse_complement()

    new_seq = parts['lsc'] + parts['ira'] + parts['ssc'] + parts['irb']
    assert len(seq_rec.seq) == len(new_seq.seq), \
        (seq_rec.name, len(seq_rec.seq), len(new_seq.seq), starts,
            [(n, len(p)) for n, p in parts.items()],
            [(n, len(p)) for n, p in partition.extract(seq_rec).items()])
    return new_seq


def orient_chloroplast_parts(seq_rec):
    partition = find_chloroplast_partition(seq)
    orientation = chloroplast_parts_orientation(seq_rec, partition)
    orientation = [p for p, is_in in orientation.items() if is_in]
    return orient_chloroplast_parts_by_data(seq_rec, orientation, partition=partition)


# Orientate by gene trnF-GAA
def trnF_GAA_start(seq_rec, partition):
    # Returns offset of trnF-GAA gene regarding start of LSC region
    l_seq = len(seq_rec.seq)
    all_genes = [f for f in seq_rec.features if f.type == 'gene' and feature_qualifiers_to_desc(f) == 'trnF-GAA']
    if not all_genes:
        # print(f'Warning: no trnF-GAA found in sequence {seq_rec.name}!')
        return

    if partition and (lsc := partition['lsc']):
        orientation = chloroplast_parts_orientation(seq_rec, partition)
        pos_oriented = orientation['lsc']
        start_p = lsc.real_start if pos_oriented else lsc.real_end

        # Gene trnF-GAA is on strand 1
        for ts in ([t for t in all_genes if (t.strand > 0) == pos_oriented],
                   [t for t in all_genes if (t.strand > 0) != pos_oriented]):
            if ts:
                t = min(ts, key=lambda t: cycle_distance_min(
                    start_p, (t.location.start if pos_oriented else t.location.end), l_seq))
                return cycle_distance_lt(start_p, t.location.start if pos_oriented else t.location.end, l_seq)

    # Offset to 0!
    return min(t.location.start for t in all_genes)


def orient_by_trnF_GAA_by_data(seq_rec, offset, orientation, partition=None, starts=None, keep_offset=None):
    if starts or partition:
        if starts:
            lsc_start = starts[0]
        elif partition:
            lsc_start = partition['lsc'].real_start
        lsc_has_offset = (abs(lsc_start) > (keep_offset or 0))
        #
        if orientation or lsc_has_offset:
            if not partition:
                partition = create_chloroplast_partition_all(len(seq_rec.seq), starts)
            if partition:
                seq_rec = orient_chloroplast_parts_by_data(seq_rec, orientation, partition=partition)
            else:
                print(f"Warning: sequence {seq_rec.name} has LSC offset but doesn't have partitions!")
    return (rotate_by_offset(seq_rec, offset, keep_offset=keep_offset) or seq_rec)


# Orientate by gene trnH-GUG
def trnH_GUG_start(seq_rec, partition):
    # Possibilities and return value:
    #  * no gene at all     : None
    #  * genome has parts   : dict(lsc_offset=<num>, strategy=<str>)
    #  * doesn't have parts : dict(zero_offset=<num>, reverse=<bool>)
    l_seq = len(seq_rec.seq)
    all_trnhs = [Feature(l_seq, feature=f) for f in seq_rec.features
                 if f.type == 'gene' and f.qualifiers['gene'][0] == 'trnH-GUG']
    if not all_trnhs:
        print(f'Warning: no trnH-GUG found in sequence {seq_rec.name}!')
        return

    if partition and (lsc := partition['lsc']):
        # Take that orientation of LSC that is calculated on more genes is more reliable then strand of only one gene!
        orientation = chloroplast_parts_orientation(seq_rec, partition)
        pos_oriented = orientation['lsc']

        # Feature data really doesn't have to be nice :-/
        # Here are strategies of finding gene location sorted by priority, with strategy tag:
        #  - (S) gene is located in LSC, near LSC start or intersects LSC start, and on right strand
        #  - (B) gene is located in IRb, near LSC start, on right strand
        #  - (A) gene is located in IRa, near LSC end, on right strand
        # In E and F cases, we assume that if there is gane in IRa,
        # there should be gene on IRb on symmetrical location, and we use that location.
        # If strategy didn't find gene on the right strand,
        # than opposite strand is also checked and 'O' is appended on strategy tag.
        # Note: trnH-GUG is short gene ~70b!

        s, e = lsc.ends()
        start_p = s if pos_oriented else e
        half_len = len(lsc) // 2
        cdm = cycle_distance_min
        trnh = None

        # (S)
        if trnhs := [t for t in all_trnhs if lsc.intersects(t) and cdm(start_p, t.real_start, l_seq) < half_len]:
            for strategy, ts in [('S', [t for t in trnhs if (t.strand < 0) == pos_oriented]),
                                 ('SO', [t for t in trnhs if (t.strand < 0) != pos_oriented])]:
                if ts:
                    trnh = min(ts, key=lambda t: cdm(start_p, (t.real_start if pos_oriented else t.real_end), l_seq))
                    return dict(strategy=strategy,
                                lsc_offset=cycle_distance_lt(start_p,
                                                             trnh.real_start if pos_oriented else trnh.real_end,
                                                             l_seq))

        # (B)
        ir_oriented = orientation['ira']
        half_len = len(partition['irb']) // 2
        if trnhs := [t for t in all_trnhs if cdm(start_p, t.real_start, l_seq) < half_len]:
            for strategy, ts in [('B', [t for t in trnhs if (t.strand < 0) == ir_oriented]),
                                 ('BO', [t for t in trnhs if (t.strand < 0) != ir_oriented])]:
                if ts:
                    trnh = min(ts, key=lambda t: cdm(start_p, (t.real_start if ir_oriented else t.real_end), l_seq))
                    return dict(strategy=strategy,
                                lsc_offset=cycle_distance_lt(start_p,
                                                             trnh.real_start if ir_oriented else trnh.real_end,
                                                             l_seq))

        # (A)
        end_p = e if pos_oriented else s
        if trnhs := [t for t in all_trnhs if cdm(end_p, t.real_start, l_seq) < half_len]:
            # Note: IRa is reverse complement of IRb, so strand is opposite
            # ?(t.real_start if ir_oriented else t.real_end) or opposite?
            return dict(strategy='A', lsc_offset=0)
            for strategy, ts in [('A', [t for t in trnhs if (t.strand < 0) != ir_oriented]),
                                 ('AO', [t for t in trnhs if (t.strand < 0) == ir_oriented])]:
                if ts:
                    trnh = min(ts, key=lambda t: cdm(end_p, (t.real_end if ir_oriented else t.real_start), l_seq))
                    # Map offset from LSC end to LSC start
                    offset = -cycle_distance_lt(end_p, trnh.real_end if ir_oriented else trnh.real_start, l_seq)
                    return dict(strategy=strategy, lsc_offset=offset)

        assert False, seq_rec.name


def orient_by_trnH_GUG_by_data(seq_rec, offset, reverse, orientation, partition=None, starts=None, keep_offset=None):
    lsc_has_offset = (abs(starts[0]) > (keep_offset or 0))
    if orientation or lsc_has_offset:
        if not partition:
            partition = create_chloroplast_partition_all(len(seq_rec.seq), starts)
        if partition:
            seq_rec = orient_chloroplast_parts_by_data(seq_rec, orientation, partition=partition)
        else:
            print(f"Warning: sequence {seq_rec.name} has LSC offset but doesn't have partitions!")
    return (rotate_by_offset(seq_rec, offset, keep_offset=keep_offset, reverse=reverse) or seq_rec)


def orient_by_trnH_GUG(seq_rec, keep_offset=None):
    # Returns None if there is no need for rotation or new SeqRecord
    partition = find_chloroplast_partition(seq_rec)
    if ret := trnH_GUG_start(seq_rec, partition):
        reverse = ret.get('reverse', False)
        offset = ret['lsc_offset'] if 'lsc_offset' in ret else ret['zero_offset']
        if partition:
            orientation = chloroplast_parts_orientation(seq_rec, partition)
            orientation = [p for p, is_in in orientation.items() if is_in]
        else:
            orientation = []
        return orient_by_trnH_GUG_by_data(
            seq_rec, offset, reverse, orientation, partition=partition, keep_offset=keep_offset)


#
def chloroplast_alignment(step_data, annotations_step, sequences, to_align, run, alignment_program, keep_offset):
    from ..alignments.common_methods import AlignmentsStep, add_sequences, run_alignment_program

    # Check sequences
    if len(sequences) == 1 and sequences[0] == 'all':
        sequences = list(annotations_step.all_sequences())
    if len(sequences) < 2:
        print('At least 2 sequences should be specified!!!')
        return

    steps = AlignmentsStep(annotations_step.project, step_data, remove_data=True)
    annotations_step.propagate_step_name_prefix(steps)

    records = [annotations_step.get_sequence_record(seq_ident) for seq_ident in sequences]

    seq_files = []
    # Whole
    if 'w' in to_align:
        seq_files.append(add_sequences(
            steps.create_substep('whole'), 'whole',
            [(seq_ident, rec.seq, None) for seq_ident, rec in zip(sequences, records)]))

    # Parts and offset
    if 'p' in to_align or 'o' in to_align:
        with_parts = [(seq_ident, rec, parts) for seq_ident, rec in zip(sequences, records)
                      if (parts := find_chloroplast_partition(rec))]
        if len(with_parts) > 1:
            if 'p' in to_align:
                for part in ('lsc', 'ira', 'ssc'):
                    seq_files.append(add_sequences(
                        steps.create_substep(part), 'gene',
                        [(seq_ident, parts[part].extract(rec).seq, None) for seq_ident, rec, parts in with_parts]))
            if 'o' in to_align:
                seq_files.append(add_sequences(
                    steps.create_substep('offset'), 'whole',
                    [(seq_ident, (rotate_to_offset(rec, parts, keep_offset=keep_offset) or rec).seq, None)
                     for seq_ident, rec, parts in with_parts]))

    # # trnH-GUG
    # if 't' in to_align:
    #     seq_files.append(add_sequences(
    #         steps.create_substep('trnH-GUG'), 'whole',
    #         [(seq_ident, (orient_by_trnH_GUG(rec, keep_offset=keep_offset) or rec).seq, None)
    #          for seq_ident, rec in zip(sequences, records)]))

    #
    run_alignment_program(alignment_program, steps, seq_files, run)
    return steps
