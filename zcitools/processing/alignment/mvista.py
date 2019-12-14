from zcitools.steps.alignments import mVISTAStep
from zcitools.utils.file_utils import write_fasta
# unzip_file, list_zip_files, write_yaml, read_yaml, \
from zcitools.utils.import_methods import import_bcbio_gff


def create_mvista_data(step_data, annotations_step, cache, run, email):
    step = mVISTAStep(step_data, remove_data=True)
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

    return step
