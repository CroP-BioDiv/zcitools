import os.path
from zcitools.utils.import_methods import import_bio_seq_io

_instructions = """
Open web page: https://chlorobox.mpimp-golm.mpg.de/OGDraw.html

For each GenBank file {calc_dir}/*.gb do:

FASTA file(s) to annotate
 * Upload file
 * (check) Circular
 * (check) Plastid
 * (check) Tidy up annotation

Inverted Repeat
 * (check) Auto

Output Options
 * check one output format type (PS prefered?)

Actions
 * Submit

When job is finished:
 - Download all results as zip (small disk icon in Results header) into {step_name}

When all files are processed:
 - run zcit command: zcit.py ogdraw {step_name}
"""


def calculate_ogdraw(step, calc_directory):
    # Check are instructions set
    i_f = step.step_calc_file(calc_directory, 'INSTRUCTIONS.txt')
    if not os.path.isfile(i_f):
        # Create gb files
        SeqIO = import_bio_seq_io()
        for seq_ident, seq_record in step._iterate_records():
            with open(step.step_calc_file(calc_directory, seq_ident + '.gb'), 'w') as out:
                SeqIO.write([seq_record], out, 'genbank')
        # Write instructions
        with open(i_f, 'w') as out:
            out.write(_instructions.format(calc_dir=step.step_file(calc_directory),
                                           step_name=step._step_name))
        return False

    # ToDo: unzip or not
    return True
