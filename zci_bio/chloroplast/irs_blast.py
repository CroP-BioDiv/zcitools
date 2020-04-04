import os.path
from common_utils.exceptions import ZCItoolsValueError
from common_utils.file_utils import ensure_directory, write_yaml, write_fasta, run_module_script, set_run_instructions
from common_utils.value_data_types import rows_2_excel
from ..utils.import_methods import import_bio_seq_io
from .utils import find_chloroplast_irs, find_referent_genome, irb_start
from .steps import ChloroplastSSCBlast
from . import run_irs_blast

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
 - BLAST+ toolkit should be on the PATH
 - It is good to use command screen for running the script.
   screen -dm "python3 {script_name}"
"""


def create_irs_data(step_data, annotation_step, params):
    SeqIO = import_bio_seq_io()

    seq_idents = annotation_step.all_sequences()  # set
    ref_ident = find_referent_genome(seq_idents, params.referent_genome)

    step = annotation_step.project.new_step(ChloroplastSSCBlast, step_data)
    ref_seq_rec = annotation_step.get_sequence_record(ref_ident)
    ssc_location = step.get_type_description_elem('ssc_location', default=dict())
    ensure_directory(step.step_file('run_dir'))
    files_to_zip = []
    calc_seq_idents = []

    # All sequences, to create database from
    for seq_ident in sorted(seq_idents):
        if not os.path.isfile(step.step_file('run_dir', f'{seq_ident}.xml')):
            fa_file = step.step_file('run_dir', f'{seq_ident}.fa')
            files_to_zip.append(fa_file)
            calc_seq_idents.append(seq_ident)
            if not os.path.isfile(fa_file):
                seq_rec = annotation_step.get_sequence_record(seq_ident)
                SeqIO.write([seq_rec], fa_file, 'fasta')
                # Store SSC position
                irs = find_chloroplast_irs(seq_rec)
                if irs:
                    ssc_location[seq_ident] = [len(seq_rec), int(irs[0].location.end), irb_start(irs[1])]
                else:
                    ssc_location[seq_ident] = [len(seq_rec), -1, -1]

    # Store query data
    query_file = step.step_file('run_dir', 'query.fa')
    if not os.path.isfile(query_file):
        irs = find_chloroplast_irs(ref_seq_rec)
        if not irs:
            raise ZCItoolsValueError(f"Referent genome ({ref_ident}) doesn't have IRS!")
        ira, irb = irs

        seq_str = str(ref_seq_rec.seq)
        bl = params.blast_length
        ssc_s = int(ira.location.end)
        ssc_e = irb_start(irb)
        write_fasta(query_file,
                    [('ira', seq_str[ssc_s - bl:ssc_s + bl]),
                     ('irb', seq_str[ssc_e - bl:ssc_e + bl])])

    if calc_seq_idents:
        # Store finish.yml
        finish_f = step.step_file('finish.yml')
        write_yaml(dict(calc_seq_idents=calc_seq_idents), finish_f)

        run = True  # ToDo: ...
        step.save(dict(ssc_location=ssc_location, blast_length=params.blast_length), completed=False)
        if run:
            run_module_script(run_irs_blast, step)
            finish_irs_data(step)
        else:
            files_to_zip.append(finish_f)
            set_run_instructions(run_irs_blast, step, files_to_zip, _instructions)
    #
    elif params.force_blast_parse:
        finish_irs_data(step)

    return step


def _middle(res):
    if not res.alignments:
        return -1
    hsp = max(res.alignments[0].hsps, key=lambda h: h.score)
    return (hsp.sbjct_start + hsp.sbjct_end) // 2


def finish_irs_data(step_object):
    from Bio.Blast import NCBIXML  # ToDo: in import methods
    ssc_location = step_object.get_type_description_elem('ssc_location')
    if not ssc_location:
        print('No SSC initial data!!!')
        return

    blast_length = step_object.get_type_description_elem('blast_length')
    check_ssc_len = 10000  # min(blast_length, 10000)
    rows = []
    for seq_ident, (seq_length, orig_start, orig_end) in sorted(ssc_location.items()):
        with open(step_object.step_file('run_dir', f'{seq_ident}.xml'), 'r') as result:
            # print(seq_ident, len(list(NCBIXML.parse(result))))
            ira_res, irb_res = list(NCBIXML.parse(result))
            ms = sorted([_middle(ira_res), _middle(irb_res)])

            rows.append([seq_ident, seq_length, orig_start, orig_end] + ms)
            rows[-1].append('NO' if ms[1] - ms[0] < check_ssc_len else '+')

    # Create table
    output_file = f'chloroplast_blast_{blast_length}.xlsx'
    columns = ['Seq', 'Length', 'Orig start', 'Orig end', 'Calc start', 'Calc end', 'OK']
    rows_2_excel(output_file, columns, rows)

    # #
    # step_object.save(None, completed=True)
