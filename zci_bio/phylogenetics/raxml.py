import re
import os
from . import run_raxml
from zci_bio.phylogenetics.steps import RAxMLStep, RAxMLSteps
from common_utils.file_utils import copy_file, unzip_file, list_zip_files, write_yaml, read_yaml, \
    run_module_script, set_run_instructions
from common_utils.exceptions import ZCItoolsValueError
from ..alignments.alignment_map_indices import AlignmentMapIndices
from ..utils.helpers import read_alignment


_re_raxml_output = re.compile(r'^RAxML_.*\.raxml_output')

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
 - RAxML executable (raxml_threads) should be on the PATH or
   environment variable RAXML_EXE should point to it.
   Real executable is probably raxmlHPC-PTHREADS-AVX2, and raxml_threads links to it.
 - It is good to use command screen for running the script.
   screen -dm "python3 {script_name}"
"""


def _copy_alignment_file(align_step, in_step, files_to_proc, set_partitions):
    a_f = in_step.step_file('alignment.phy')
    copy_file(align_step.get_phylip_file(), a_f)
    alignment = read_alignment(a_f)
    partitions = None
    #
    if set_partitions and align_step.is_composite():
        partitions = in_step.step_file('partitions.ind')
        ami = AlignmentMapIndices(alignment=alignment)
        ami = align_step.get_alignment_map_indices()
        ami.create_raxml_partitions(align_step.iterate_partitions(), partitions)
    #
    files_to_proc.append(dict(
        filename=a_f, short=align_step.is_short(), length=len(alignment[0]), partitions=partitions))


def create_raxml_data(step_data, alignment_step, cache, set_partitions, run):
    # List of dicts with attrs: filename, short, partitions (filename or None)
    # This data is used to optimize calculation
    files_to_proc = []

    if alignment_step._IS_COLLECTION:
        step = RAxMLSteps(alignment_step.project, step_data, remove_data=True)
        for align_step in alignment_step.step_objects():
            substep = step.create_substep(align_step.get_local_name())
            substep.set_sequences(align_step.all_sequences())
            substep.seq_sequence_type(align_step.get_sequence_type())
            _copy_alignment_file(align_step, substep, files_to_proc, set_partitions)
            #
            substep.save(completed=False)
    else:
        step = RAxMLStep(alignment_step.project, step_data, remove_data=True)
        step.set_sequences(alignment_step.all_sequences())
        step.seq_sequence_type(alignment_step.get_sequence_type())
        _copy_alignment_file(alignment_step, step, files_to_proc, set_partitions)

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
        run_module_script(run_raxml, step)
    else:
        files_to_zip.append(finish_f)
        set_run_instructions(run_raxml, step, files_to_zip, _instructions)
    #
    return step


def finish_raxml_data(step_obj, cache):
    output_f = step_obj.step_file('output.zip')
    if not os.path.isfile(output_f):
        raise ZCItoolsValueError('No calculation output file output.zip!')

    # Check are all file RAxML outputs, in same directories as files to process and
    # filenames matches RAxML_.*\.raxml_output
    dirs = set(os.path.dirname(d['filename']) for d in read_yaml(step_obj.step_file('finish.yml')))
    for z_file in list_zip_files(output_f):
        parts = z_file.split('/')  # ZipFile uses '/' as separator
        _dir = '' if len(parts) == 1 else os.sep.join(parts[:-1])
        if _dir not in dirs:
            raise ZCItoolsValueError(f'Output contains file(s) in not step directory ({_dir})!')

        if not _re_raxml_output.search(parts[-1]):
            raise ZCItoolsValueError(f'Not RAxML output file(s)found in the output ({parts[-1]})!')

    # Unzip data
    unzip_file(output_f, step_obj.directory)

    step_obj._check_data()
    step_obj.save(create=False)
