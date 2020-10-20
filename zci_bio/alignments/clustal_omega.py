import os.path
from . import run_clustal_omega
from .common_methods import create_alignment_data, finish_alignment_data

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
 - Clustal Omega executable (clustalo) should be on the PATH or
   environment variable CLUSTAL_OMEGA_EXE should point to it.
 - It is good to use command screen for running the script.
   screen -dm "python3 {script_name}"
"""


def create_clustal_data(step_data, annotations_step, alignments, whole_partition, run):
    return create_alignment_data(
        step_data, annotations_step, alignments, whole_partition, run, run_clustal_omega, _instructions)


def finish_clustal_data(step_obj):
    # Check are needed files in zip, not something strange
    files = set(d['filename'].replace('sequences.fa', 'alignment.phy')
                for d in read_yaml(step_obj.step_file('finish.yml')))
    finish_alignment_data(step_obj, files)
