from common_utils.exceptions import ZCItoolsValueError
from ..utils.features import Feature, Partition


def find_chloroplast_irs(seq):
    # Finds the longest pair of inverted repeats
    _ir = ('inverted',)
    rep_regs = [f for f in seq.features
                if f.type == 'repeat_region' and
                f.qualifiers.get('rpt_type', _ir)[0] == 'inverted']
    if rep_regs:
        max_len = max(map(len, rep_regs)) - 3  # Some tolerance :-)
        max_regs = [f for f in rep_regs if len(f) >= max_len]
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
def chloroplast_alignment(step_data, annotations_step, sequences, run, alignment_program):
    from ..alignments.common_methods import AlignmentsStep, add_sequences, run_alignment_program

    # Check sequences
    if len(sequences) < 2:
        print('At least 2 sequences should be specified!!!')
        return

    steps = AlignmentsStep(annotations_step.project, step_data, remove_data=True)
    annotations_step.propagate_step_name_prefix(steps)

    records = [annotations_step.get_sequence_record(seq_ident) for seq_ident in sequences]
    with_parts = [(seq_ident, rec, parts) for seq_ident, rec in zip(sequences, records)
                  if (parts := find_chloroplast_partition(rec))]
    print(with_parts)

    seq_files = []
    # Whole
    seq_files.append(add_sequences(
        steps.create_substep('whole'), 'whole',
        [(seq_ident, rec.seq, None) for seq_ident, rec in zip(sequences, records)]))

    # Parts
    if len(with_parts) > 1:
        for part in ('lsc', 'ira', 'ssc'):
            seq_files.append(add_sequences(
                steps.create_substep(part), 'gene',
                [(seq_ident, parts[part].extract(rec).seq, None) for seq_ident, rec, parts in with_parts]))

    #
    run_alignment_program(alignment_program, steps, seq_files, run)
    return steps
