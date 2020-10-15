import os.path
from . import run_muscle
from .common_methods import create_alignment_data

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
 - MUSCLE executable (muscle) should be on the PATH or
   environment variable MUSCLE_EXE should point to it.
 - It is good to use command screen for running the script.
   screen -dm "python3 {script_name}"
"""


def create_muscle_data(step_data, annotations_step, alignments, whole_partition, run):
    return create_alignment_data(
        step_data, annotations_step, alignments, whole_partition, run, run_muscle, _instructions)
