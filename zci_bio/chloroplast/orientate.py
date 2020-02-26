from common_utils.file_utils import copy_file, write_fasta
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
    step = GroupedSequencesStep(annotation_step.project, step_data, remove_data=True)
    good = step.create_substep('good')
    repaired = step.create_substep('repaired')
    bad = step.create_substep('bad')

    for seq_ident, seq in annotation_step._iterate_records():
        partition = find_regions(seq_ident, seq)
        if not partition:
            # ToDo: what to do? Check is orientation good with list of genes. If not raise 
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

        an_file = annotation_step.get_sequence_filename()
        if (lsc_count <= 0) and (ssc_count <= 0) and (ir_count >= 0):
            # No change
            step_f = good.step_file(os.path.basename(an_file))
            copy_file(an_file, step_f)
            good.add_sequence_file(step_f)
        else:
            parts = partition.extract(seq)  # dict name -> Seq object
            if ssc_count > 0:
                parts['ssc'] = parts['ssc'].reverse_complement()
            if lsc_count > 0:
                assert False, f'REVERT LSC {seq_ident}'
            if ir_count < 0:
                assert False, f'REVERT IRs {seq_ident}'

            new_seq = parts['lsc'] + parts['ira'] + parts['ssc'] + parts['irb']
            step_f = repaired.step_file(os.path.basename(an_file))
            write_fasta(step_f, [(seq_ident, new_seq.seq)])
            repaired.add_sequence_file(step_f)

    #
    good.save()
    repaired.save()
    bad.save()
    step.save()
    return step
