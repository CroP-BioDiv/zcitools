import re
from zipfile import ZipFile
from zcitools.steps.annotations import AnnotationsStep
from zcitools.utils.file_utils import copy_file  # link_file
from zcitools.utils.helpers import split_sequences, split_list

_re_zip_genbank = re.compile('GeSeqJob-[0-9]*-[0-9]*_(.*)_GenBank.gb')

_instructions = """
Open web page: https://chlorobox.mpimp-golm.mpg.de/geseq.html

For each FASTA file {calc_dir}/sequences_<num>.fa do:

FASTA file(s) to annotate
 * Upload file: {step_name}/sequences.fa
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


def create_ge_seq_data(step_data, sequences_step):
    step = AnnotationsStep(step_data, remove_data=True)
    cache = step.get_cache_object()
    all_sequences = list(sequences_step.all_sequences())

    # Fetch cached sequences
    to_fetch = step.get_cached_records(cache, all_sequences, info=True)

    # Store sequence
    if to_fetch:
        for i, d in enumerate(split_list(to_fetch, 30)):
            sequences_step.concatenate_seqs_fa(step.step_file(f'sequences_{i + 1}.fa'), d)

        # Store instructions
        with open(step.step_file('INSTRUCTIONS.txt'), 'w') as out:
            out.write(_instructions.format(step_name=step_data['step_name']))

    #
    step.set_sequences(all_sequences)
    step.save(needs_editing=True)
    return step


def finish_ge_seq_data(step_obj):
    # Note: original files are left in directory
    # ToDo: inverted_region 126081..1 !!! To_ind > from_ind!!!
    # First check job-results-<num>.zip file
    job_files = step_obj.step_files(matches='^job-results-[0-9]*.zip')
    if job_files:
        # Extract GenBank files named job-results-<num>/GeSeqJob-<num>-<num>_<seq_ident>_GenBank.gb
        for filename in job_files:
            with ZipFile(step_obj.step_file(filename), 'r') as zip_f:
                for z_i in zip_f.infolist():
                    m = _re_zip_genbank.search(z_i.filename)
                    if m:
                        seq_ident = m.group(1)
                        with open(step_obj.step_file(seq_ident + '.gb'), 'wb') as save_f:
                            save_f.write(zip_f.read(z_i.filename))

    else:
        # Check file named: GeSeqJob-<num>-<num>_GLOBAL_multi-GenBank.gbff
        job_files = step_files(matches='^GeSeqJob-.*_GLOBAL_multi-GenBank.gbff')
        if not job_files:
            print("Warning: can't find any GeSeq output file! Nor job-results-*.zip, nor GeSeqJob-.*.gbff file(s).")
            return

        for filename in job_files:
            split_sequences(filename, '.gb')

    step_obj._check_data()
    step_obj.save()
