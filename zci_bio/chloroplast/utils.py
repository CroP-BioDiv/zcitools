from common_utils.exceptions import ZCItoolsValueError
from ..utils.features import Feature, Partition


def find_chloroplast_irs(seq, check_size=True):
    # Finds the longest pair of inverted repeats
    _ir = ('inverted',)
    rep_regs = [f for f in seq.features
                if f.type == 'repeat_region' and
                f.qualifiers.get('rpt_type', _ir)[0] == 'inverted']
    if seq.name == 'NC_033344':
        print(seq.name, rep_regs)
    if len(rep_regs) >= 2:
        max_len = max(map(len, rep_regs)) - 3  # Some tolerance :-)
        max_regs = [f for f in rep_regs if len(f) >= max_len]
        if len(max_regs) != 2 and not check_size:
            # Backup
            max_regs = sorted(rep_regs, key=len)[-2:]
        if len(max_regs) == 2:
            ira, irb = max_regs
            diff_1 = (irb.location.parts[0].start - ira.location.parts[-1].end) % len(seq)
            diff_2 = (ira.location.parts[0].start - irb.location.parts[-1].end) % len(seq)
            return (ira, irb) if diff_1 < diff_2 else (irb, ira)


def irb_start(irb):
    return int(irb.location.parts[0].start)


def find_chloroplast_partition(seq):
    # Returns None or Partition object with parts named: lsc, ira, ssc, irb.
    irs = find_chloroplast_irs(seq)
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


# -----------------
def cycle_distance(a, b, cycle_len):
    d = abs(a - b)
    return d if d < (cycle_len / 2) else (cycle_len - d)


def rotate_by_offset(seq_rec, offset, keep_offset=None, reverse=False):
    # Returns None if there is no need for rotation or new SeqRecord
    if offset is None:
        return seq_rec.reverse_complement() if reverse else None

    offset %= len(seq_rec.seq)

    # Check if there is need to rotate
    if not offset or (keep_offset and cycle_distance(0, offset, len(seq_rec.seq)) <= keep_offset):
        return seq_rec.reverse_complement() if reverse else None

    new_seq = seq_rec[offset:] + seq_rec[:offset]
    return new_seq.reverse_complement() if reverse else new_seq


def trnH_GUG_start(seq_rec, partition):
    # Returns description how to offset sequence to get gene trnH-GUG on it's start, tuple (offset, bool to_reverse)
    # Returns None if there is no trnH-GUG gene in given record.
    trnhs = [f for f in seq_rec.features if f.type == 'gene' and f.qualifiers['gene'][0] == 'trnH-GUG']
    if not trnhs:
        print(f'Warning: no trnH-GUG found in sequence {seq_rec.name}!')
        return

    l_seq = len(seq_rec.seq)

    if partition and (lsc := partition['lsc']):
        # Take one that is the closest to some lsc end
        s, e = lsc.ends()
        mins = [min(min(cycle_distance(s, p, l_seq), cycle_distance(e, p, l_seq))
                    for p in (t.location.start, t.location.end)) for t in trnhs]
        gl_min = min(mins)
        idx = mins.index(gl_min)
        trnh = trnhs[idx]
        # trnH-GUG is located at LSC start. Check on which side it is now
        reverse = (cycle_distance(e, trnh.location.end, l_seq) == gl_min)
        offset = trnh.location.start if reverse else trnh.location.end
        return offset, reverse

    # Take one that is the closest to the origin
    trnh = min(trnhs, key=lambda f: min(cycle_distance(0, p, l_seq) for p in (f.location.start, f.location.end)))
    # If genome is good orineted than strand should be -1
    return (trnh.location.end, True) if trnh.location.strand > 0 else (trnh.location.start, False)


def rotate_to_trnH_GUG(seq_rec, keep_offset=None):
    # Returns None if there is no need for rotation or new SeqRecord
    if ret := trnH_GUG_start(seq_rec, find_chloroplast_partition(seq_rec)):
        offset, reverse = ret
        return rotate_by_offset(seq_rec, offset, keep_offset=keep_offset, reverse=reverse)


def rotate_to_offset(seq_rec, parts, keep_offset=None):
    return rotate_by_offset(seq_rec, parts['lsc'].real_start, keep_offset=keep_offset)


def chloroplast_alignment(step_data, annotations_step, sequences, to_align, run, alignment_program, keep_offset):
    from ..alignments.common_methods import AlignmentsStep, add_sequences, run_alignment_program

    # Check sequences
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

    # trnH-GUG
    if 't' in to_align:
        seq_files.append(add_sequences(
            steps.create_substep('trnH-GUG'), 'whole',
            [(seq_ident, (rotate_to_trnH_GUG(rec, keep_offset=keep_offset) or rec).seq, None)
             for seq_ident, rec in zip(sequences, records)]))

    #
    run_alignment_program(alignment_program, steps, seq_files, run)
    return steps
