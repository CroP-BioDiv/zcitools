from ..steps.annotations import AnnotationsStep
from ..utils.file_utils import copy_file  # link_file

_instructions = """
Open web page: https://chlorobox.mpimp-golm.mpg.de/geseq.html

FASTA file(s) to annotate
 * Upload file: sequences.fa
 * (check) Circular
 * (check) Sequence source: Plastid

Options
 * (check) Generate multi-GenBank
 * uncheck other
 * ?(check) Generate multi-GFF3

Annotation
 BLAT search
   * (check) Annotate plastid IR
   * (check) Annotate plastid trans-spliced rps12
 3rd Party tRNA annotators
   * (check) ARAGORN (default settings)
   * (check) tRNAscan-SE (default settings)

BLAT Reference Sequences
   * (check) MPI-MP chloroplast references

Actions
 * Submit

When job is finished:
 - download Global multi-GenBank file into job directory ({abspath})
 - run zcit command: zcit finish {step_name}

Documentation:
https://chlorobox.mpimp-golm.mpg.de/gs_documentation.html
"""


def create_ge_seq_data(step_data, sequences_step):
    step = AnnotationsStep(step_data, remove_data=True)

    # Store sequence
    all_seqs_fa = sequences_step.get_all_seqs_fa()
    # Note: browser uploading sometimes do not work with links :-/
    copy_file(all_seqs_fa, step.step_file('sequences.fa'))

    # Store instructions
    with open(step.step_file('INSTRUCTIONS.txt'), 'w') as out:
        out.write(_instructions.format(abspath=step.absolute_path(), step_name=step_data['step_name']))

    #
    step.set_sequences(sequences_step.all_sequences())
    step.save(needs_editing=True)
    return step


def finish_ge_seq_data(step_obj):
    # Check file named: GeSeqJob-<num>-<num>_GLOBAL_multi-GenBank.gbff
    for f in step_obj.step_files():
        if f.startswith('GeSeqJob') and f.endswith('_GLOBAL_multi-GenBank.gbff'):
            filename = f
            break
    else:
        print("Warning: can't find GeSeq output file!")
        return

    # Leave original file
    copy_file(step_obj.step_file(filename), step_obj.get_all_annotation_filename())
    step_obj.save()
