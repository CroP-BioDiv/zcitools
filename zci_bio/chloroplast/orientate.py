from ..sequences.steps import SequencesStep
from .utils import find_regions
from ..utils.features import Feature

# Orientates chloroplast sequence in standard way.
# Uses Fast-Plast method. Check
#  - source file orientate_plastome_v.2.0.pl
#    (https://github.com/mrmckain/Fast-Plast/blob/master/bin/orientate_plastome_v.2.0.pl)
#  - explanation https://github.com/mrmckain/Fast-Plast/issues/22
# Consitent with Wikipedia image:
#  - https://en.wikipedia.org/wiki/Chloroplast_DNA#/media/File:Plastomap_of_Arabidopsis_thaliana.svg


def orientate_chloroplast(step_data, annotation_step):
    step = SequencesStep(annotation_step.project, step_data, remove_data=True)

    for seq_ident, seq in annotation_step._iterate_records():
        partition = find_regions(seq_ident, seq)
        if not partition:
            # ToDo: what to do?
            assert False, seq_ident
            continue

        l_seq = len(seq)
        in_parts = partition.put_features_in_parts(
            Feature(l_seq, feature=f) for f in seq.features if f.type == 'gene' and f.location)

        # Count gene orintation
        lsc_count = 0
        for f in in_parts.get('lsc', []):
            if any(x in f.name for x in ('rpl', 'rps')):
                lsc_count += f.feature.strand

        ssc_count = 0
        for f in in_parts.get('ssc', []):
            ssc_count += f.feature.strand

        ir_count = 0
        for f in in_parts.get('ira', []):
            if 'rrn' in f.name:
                ir_count += f.feature.strand

        # Revert if needed
        if lsc_count > 0:
            assert False, f'REVERT LSC {seq_ident}'
        if ssc_count > 0:
            print('REVERT SSC', seq_ident)
        if ir_count < 0:
            assert False, f'REVERT IRs {seq_ident}'

        # # Ima ili nema IR-ova
        # step.add_sequence_file(f)

    step.save()
    return step
