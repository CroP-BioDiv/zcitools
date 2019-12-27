from step_project.steps.alignments import mVISTAStep
from step_project.utils.import_methods import import_bcbio_gff
from common_utils.file_utils import write_fasta

_instructions_no_run = """
Steps:
 - open web page http://genome.lbl.gov/vista/mvista/submit.shtml
 - in field 'Total number of sequences' set {num_sequences}
 - press Submit button
 - On next page
 - in field 'Your email address' set email address to get informed when calculation finish
 - change directory to {step_name}
 - in 'Sequence #<num>' upload fasta files one by one. Take a care about file order!
 - in 'Additional options' set
   - for 'Alignment program' set Shuffle-LAGAN (?)
   - in Sequence #<num> parts for each sequence
     - set it's filename (without .fa)
     - for Annotation uploaf corresponding gff3 file
 - press Submit button

"""


def create_mvista_data(step_data, annotations_step, cache, run, email):
    step = mVISTAStep(annotations_step.zcit, step_data, remove_data=True)
    sequences = sorted(annotations_step.all_sequences())

    # Store fasta and gff3 files
    gff = import_bcbio_gff()
    for seq_ident, seq_record in annotations_step._iterate_records():
        write_fasta(step.step_file(f'{seq_ident}.fa'), [(seq_ident, seq_record.seq)])
        # Create new SeqRecord with features of interest
        sr = seq_record.__class__(seq_record.seq, seq_ident)
        sr.features = [f for f in seq_record.features if f.type == 'gene']
        with open(step.step_file(f'{seq_ident}.gff3'), 'w') as output:
            gff.write([sr], output)

    # Jel run ili ne
    return step
