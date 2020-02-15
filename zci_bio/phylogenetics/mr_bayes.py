import os
from . import run_mr_bayes
from zci_bio.phylogenetics.steps import MrBayesStep, MrBayesSteps
from common_utils.file_utils import copy_file, unzip_file, list_zip_files, write_yaml, read_yaml, \
    run_module_script, set_run_instructions
from common_utils.exceptions import ZCItoolsValueError

_NEXUS_DATA = """
begin mrbayes;
    set autoclose=yes nowarn=yes autoreplace=no;
    lset Nst=6 Rates=gamma;
    mcmcp ngen=10000000 printfreq=1000 samplefreq=1000 nchains=4
    savebrlens=yes filename=alignment;
    mcmc;
    sumt filename=alignment burnin=2500 contype=halfcompat;
end;

"""

_instructions = """
Steps:
 - copy file calculate.zip onto server
 - unzip it
 - change directory to {step_name}
 - run script: {script_name}
    - to specify number of threads to use run: {script_name} <num_threads>
      default is number of cores.
 - copy file output.zip back into project's step directory {step_name}
 - run zcit command: zcit.py finish {step_name}

Notes:
 - MrBayes executable (mr_bayes) should be on the PATH or
   environment variable MR_BAYES_EXE should point to it.
 - It is good to use command screen for running the script.
   screen -dm "python3 {script_name}"
"""


def _copy_alignment_file(align_step, in_step, files_to_proc):
    # ToDo: Nexus i dodati nes na kraj!
    a_f = in_step.step_file('alignment.nex')
    copy_file(align_step.get_nexus_file(), a_f)
    # Add MrBayes data
    with open(a_f, 'a') as output:
        output.write(_NEXUS_DATA)
    #
    files_to_proc.append(dict(filename=a_f, short=align_step.is_short()))


def create_mr_bayes_data(step_data, alignment_step, cache, run):
    # List of dicts with attrs: filename, short
    # This data is used to optimize calculation
    # ToDo: almost the same as raxml.py. Differs in class types, _copy_alignment_file() and file formats
    files_to_proc = []

    if alignment_step._IS_COLLECTION:
        step = MrBayesSteps(alignment_step.project, step_data, remove_data=True)
        for align_step in alignment_step.step_objects():
            substep = step.create_substep(align_step.get_local_name())
            substep.set_sequences(align_step.all_sequences())
            substep.seq_sequence_type(align_step.get_sequence_type())
            _copy_alignment_file(align_step, substep, files_to_proc)
            #
            substep.save(completed=False)
    else:
        step = MrBayesStep(alignment_step.project, step_data, remove_data=True)
        step.set_sequences(alignment_step.all_sequences())
        step.seq_sequence_type(alignment_step.get_sequence_type())
        _copy_alignment_file(alignment_step, step, files_to_proc)

    # Store files desc
    files_to_zip = [d['filename'] for d in files_to_proc]  # files to zip
    # Remove step directory from files since run script is called from step directory
    for d in files_to_proc:
        d['filename'] = step.strip_step_dir(d['filename'])
    finish_f = step.step_file('finish.yml')
    write_yaml(files_to_proc, finish_f)

    # Stores description.yml
    step.save(completed=run)

    if run:
        run_module_script(run_mr_bayes, step)
    else:
        files_to_zip.append(finish_f)
        set_run_instructions(run_mr_bayes, step, files_to_zip, _instructions)
    #
    return step


def finish_mr_bayes_data(step_obj, cache):
    output_f = step_obj.step_file('output.zip')
    if not os.path.isfile(output_f):
        raise ZCItoolsValueError('No calculation output file output.zip!')

    allowed_files = ('alignment.ckp', 'alignment.con.tre', 'alignment.parts', 'alignment.tstat', 'alignment.vstat')

    # Check are all file MrBayes outputs
    dirs = set(os.path.dirname(d['filename']) for d in read_yaml(step_obj.step_file('finish.yml')))
    for z_file in list_zip_files(output_f):
        parts = z_file.split('/')  # ZipFile uses '/' as separator
        _dir = '' if len(parts) == 1 else os.sep.join(parts[:-1])
        if _dir not in dirs:
            raise ZCItoolsValueError(f'Output contains file(s) in not step directory ({_dir})!')

        if parts[-1] not in allowed_files:
            raise ZCItoolsValueError(f'Not RAxML output file(s)found in the output ({parts[-1]})!')

    # Unzip data
    unzip_file(output_f, step_obj.directory)

    step_obj._check_data()
    step_obj.save(create=False)
