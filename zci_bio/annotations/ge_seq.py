import re
from zipfile import ZipFile
from .steps import AnnotationsStep
from ..utils.helpers import split_sequences
from common_utils.misc import split_list
from common_utils.file_utils import extract_from_zip, write_str_in_file

_re_zip_genbank = re.compile('GeSeqJob-[0-9]*-[0-9]*_(.*)_GenBank.gb')

_instructions = """
Open web page: https://chlorobox.mpimp-golm.mpg.de/geseq.html

For each FASTA file sequences_<num>.fa do:

FASTA file(s) to annotate
 * Upload file: {step_name}/sequences_<num>.fa
 * (check) Circular
 * (check) Sequence source: Plastid

Options
 * (check) Generate multi-GenBank
 * uncheck other

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
 - download all results as zip (small disk icon in Results header) into job directory ({step_name})
   - alternatively if OGDraw images not needed, download Global multi-GenBank file into job directory ({step_name})

When all files are processed:
 - run zcit command: zcit.py finish {step_name}

-----
Documentation: https://chlorobox.mpimp-golm.mpg.de/gs_documentation.html
"""


def create_ge_seq_data(step_data, sequences_step, common_db, num_sequences_in_file):
    step = AnnotationsStep(sequences_step.project, step_data, remove_data=True)
    sequences_step.propagate_step_name_prefix(step)
    all_sequences = list(sequences_step.all_sequences())

    # Fetch common DB sequences
    to_fetch = step.get_common_db_records(common_db, all_sequences, info=True)

    # Store sequence
    if to_fetch:
        for i, d in enumerate(split_list(to_fetch, num_sequences_in_file)):
            sequences_step.concatenate_seqs_fa(step.step_file(f'sequences_{i + 1}.fa'), d)

        # Store instructions
        write_str_in_file(step.step_file('INSTRUCTIONS.txt'), _instructions.format(step_name=step_data['step_name']))

    #
    step.set_sequences(all_sequences)
    step.save(completed=not to_fetch)
    return step


def finish_ge_seq_data(step_obj, common_db):
    # Note: original files are left in directory
    # ToDo: inverted_region 126081..1 !!! To_ind > from_ind!!!
    # First check job-results-<num>.zip file
    job_files = step_obj.step_files(matches='^job-results-[0-9]*.zip')
    added_seqs = []
    if job_files:
        # Extract GenBank files named job-results-<num>/GeSeqJob-<num>-<num>_<seq_ident>_GenBank.gb
        for filename in job_files:
            with ZipFile(step_obj.step_file(filename), 'r') as zip_f:
                for z_i in zip_f.infolist():
                    m = _re_zip_genbank.search(z_i.filename)
                    if m:
                        seq_ident = m.group(1)
                        added_seqs.append(seq_ident)
                        extract_from_zip(zip_f, z_i.filename, step_obj.step_file(seq_ident + '.gb'))

    else:
        # Check file named: GeSeqJob-<num>-<num>_GLOBAL_multi-GenBank.gbff
        job_files = step_obj.step_files(matches='^GeSeqJob-.*_GLOBAL_multi-GenBank.gbff')
        if not job_files:
            print("Warning: can't find any GeSeq output file! Nor job-results-*.zip, nor GeSeqJob-.*.gbff file(s).")
            return

        for filename in job_files:
            added_seqs.extend(split_sequences(step_obj.step_file(filename), '.gb'))

    # Set into the common DB
    # Note: step's check data can raise an exeption. It is better first to store into a cache
    if common_db:
        for seq_ident in added_seqs:
            common_db.set_record(seq_ident, step_obj.step_file(seq_ident + '.gb'))

    step_obj._check_data()
    step_obj.save(create=False)
