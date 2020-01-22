import os.path
import random
import itertools
from common_utils.file_utils import ensure_directory, copy_file, write_yaml, run_module_script, set_run_instructions
from .steps import NewHybridsStep
from . import run_new_hybrids

_MAX_SMALL_NUMBER = 9

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
 - NewHybrids executable (newhybrids-no-gui-linux.exe) should be on the PATH or
   environment variable NEW_HYBRIDS_EXE should point to it.
 - It is good to use command screen for running the script.
   screen -dm "python3 {script_name}"
"""


def create_new_hybrids_data(project, step_data, params):
    # Check input files
    if not os.path.isfile(params.data_file):
        raise ZCItoolsValueError(f"Input data file {params.data_file} doesn't exist!")
    if not os.path.isfile(params.gtyp_cat_file):
        raise ZCItoolsValueError(f"Input genotype category probabilities {params.gtyp_cat_file} doesn't exist!")
    data_file = os.path.basename(params.data_file)
    gtyp_cat_file = os.path.basename(params.gtyp_cat_file)

    step = NewHybridsStep(project, step_data, remove_data=True)
    step.set_data(data_file, gtyp_cat_file, params.theta_prior, params.pi_prior, params.burn_in,
                  params.num_sweeps)

    # Copy input files
    files_to_zip = [step.step_file(data_file), step.step_file(gtyp_cat_file)]
    copy_file(params.data_file, files_to_zip[0])
    copy_file(params.gtyp_cat_file, files_to_zip[1])

    # Create run directories
    seeds = random.sample(list(itertools.product(range(1, _MAX_SMALL_NUMBER + 1), repeat=2)), params.num_runs)

    for seed in seeds:
        files_to_zip.append(step.step_file(step.seed_dir(seed)))
        ensure_directory(files_to_zip[-1])

    files_to_zip.append(step.step_file('finish.yml'))
    write_yaml(dict(
        data_file=data_file,
        gtyp_cat_file=gtyp_cat_file,
        theta_prior=params.theta_prior,
        pi_prior=params.pi_prior,
        burn_in=params.burn_in,
        num_sweeps=params.num_sweeps), files_to_zip[-1])

    # Stores description.yml
    step.save(completed=params.run)

    # Run or set instructions
    if params.run:
        run_module_script(run_new_hybrids, step)
    else:
        set_run_instructions(run_new_hybrids, step, files_to_zip, _instructions)
    #
    return step


def finish_new_hybrids(step_obj):
    pass
