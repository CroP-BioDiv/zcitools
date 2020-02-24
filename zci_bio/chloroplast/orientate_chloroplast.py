from ..sequences.steps import SequencesStep

# Orientates chloroplast sequence in standard way.
# Uses Fast-Plast method. Check
#  - source file orientate_plastome_v.2.0.pl
#    (https://github.com/mrmckain/Fast-Plast/blob/master/bin/orientate_plastome_v.2.0.pl)
#  - explanation https://github.com/mrmckain/Fast-Plast/issues/22
# Consitent with Wikipedia image:
#  - https://en.wikipedia.org/wiki/Chloroplast_DNA#/media/File:Plastomap_of_Arabidopsis_thaliana.svg


def orientate_chloroplast(step_data, annotation_step):
    step = SequencesStep(step_data)

    for seq_ident, seq in annotation_step._iterate_records():
        # Ima ili nema IR-ova
        step.add_sequence_file(f)

    step.save()
    return step
