import os.path
from collections import OrderedDict
from common_utils.file_utils import write_fasta, run_module_script, set_run_instructions
from step_project.common.table.steps import TableStep
from zci_bio.annotations.steps import AnnotationsStep
from . import run_inverted_repeats

_instructions = """
Steps:
 - copy file calculate.zip onto server
 - unzip it
 - change directory to {step_name}
 - run script: python3 {script_name}
    - to specify number of threads to use run: python3 {script_name} <num_threads>
      default is number of cores.
 - copy file output.zip back into project's step directory {step_name}
 - run zcit command: zcit.py finish {step_name}

Notes:
 - Mummers repeat match executable (repeat-match) should be on the PATH or
   environment variable MUMMER_REPEAT_MATCH_EXE should point to it.
 - It is good to use command screen for running the script.
   screen -dm "python3 {script_name}"
"""


def create_irs_data(step_data, input_step, run, table_output=False):
    # List of dicts with attrs: filename, short, partitions (filename or None)
    # This data is used to optimize calculation
    if input_step.step_data_type == 'sequences':
        table_output = True

    if table_output:
        step = TableStep(input_step.project, step_data, remove_data=True)
    else:
        step = AnnotationsStep(input_step.project, step_data, remove_data=True)
        # Set sequences
        self.set_sequences(input_step.all_sequences())
        for seq_ident in input_step.all_sequences():
            in_gb = input_step.get_sequence_filename(seq_ident)
            copy_file(in_gb, step.step_file(os.path.basename(in_gb)))

    files_to_zip = []

    for seq_ident, seq_rec in input_step._iterate_records():
        files_to_zip.append(step.step_file(f'{seq_ident}.fa'))
        write_fasta(files_to_zip[-1], [(seq_ident, str(seq_rec.seq))])

    # Stores description.yml
    step.save(completed=run)

    if run:
        run_module_script(run_inverted_repeats, step)
        finish_irs_data(step)
    else:
        set_run_instructions(run_inverted_repeats, step, files_to_zip, _instructions)

    #
    return step


def _read_irs_files():
    for f in os.listdir('.'):
        if f.endswith('.irs'):
            with open(f, 'r') as r:
                yield (f[:-4],) + tuple(tuple(map(int, line.strip().split())) for line in r.readlines())


def finish_irs_data(step_obj):
    if step_obj.step_data_type == 'table':
        rows = [[seq_ident, *ira, *irb] for seq_ident, ira, irb in _read_irs_files()]
        step_obj.set_table_data(rows, [('seq_ident', 'seq_ident'),
                                       ('ira_start', 'int'), ('ira_end', 'int'),
                                       ('irb_start', 'int'), ('irb_end', 'int')])
    else:
        assert step_obj.step_data_type == 'annotations', step_obj.step_data_type
        from ..utils.import_methods import import_bio_seq_io
        SeqIO = import_bio_seq_io()
        from Bio.SeqFeature import SeqFeature, FeatureLocation

        for seq_ident, ira, irb in _read_irs_files():
            seq_rec = step_obj.get_sequence_record(seq_ident)
            qs = OrderedDict([('rpt_type', 'inverted')])
            seq_rec.features.append(SeqFeature(FeatureLocation(*ira), type='repeat_region'), qualifiers=qs)
            seq_rec.features.append(SeqFeature(FeatureLocation(*irb), type='repeat_region'), qualifiers=qs)
            #
            SeqIO.write([seq_rec], step_obj.step_file(f'{seq_ident}.gb'), 'genbank')

    step_obj.save()
