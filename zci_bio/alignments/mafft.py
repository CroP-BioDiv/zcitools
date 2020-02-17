import os.path
from . import run_mafft
from common_utils.misc import sets_equal
from common_utils.file_utils import unzip_file, list_zip_files, read_yaml
from common_utils.exceptions import ZCItoolsValueError

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
 - MAFFT executable (mafft) should be on the PATH or
   environment variable MAFFT_EXE should point to it.
 - It is good to use command screen for running the script.
   screen -dm "python3 {script_name}"
"""


def create_clustal_data(step_data, annotations_step, alignments, run):
    return create_alignment_data(step_data, annotations_step, alignments, run, run_mafft, _instructions)
