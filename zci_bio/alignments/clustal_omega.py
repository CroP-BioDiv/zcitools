import os.path
from . import run_clustal_omega
from .steps import AlignmentStep, AlignmentsStep
from common_utils.misc import sets_equal
from common_utils.file_utils import unzip_file, list_zip_files, write_yaml, read_yaml, \
    write_fasta, run_module_script, set_run_instructions
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
 - Clustal Omega executable (clustalo) should be on the PATH or
   environment variable CLUSTAL_OMEGA_EXE should point to it.
 - It is good to use command screen for running the script.
   screen -dm "python3 {script_name}"
"""

# Note: global cache is not used!


def _add_sequences(al_step, seq_type, sequences, sequence_data):
    seq_file = al_step.step_file('sequences.fa')
    write_fasta(seq_file, ((i, s) for i, s, _ in sequence_data))
    al_step.set_sequences(sequences)
    al_step.seq_sequence_type(seq_type)
    al_step.save(completed=False)
    # Annotations
    for seq_ident, seq, partitions in sequence_data:
        if partitions:
            al_step.store_partitions(seq_ident, partitions)
    #
    return dict(filename=seq_file, short=al_step.is_short(), max_seq_length=max(len(s) for _, s, _ in sequence_data))


def _lengths_2_partitions(data):  # data is list of strings (ident, string)
    ret = []
    end = 1
    for ident, d in data:
        end += len(d)
        ret.append((end, ident))
    return ret


def _create_al_step(created_steps, steps, name, annotations_step, step_data):
    if steps:
        al_step = steps.create_substep(name)
    else:
        al_step = AlignmentStep(annotations_step.project, step_data, remove_data=True)
    created_steps.append(al_step)
    return al_step


def _feature_sequences(
        created_steps, step_data, steps, annotations_step, sequences, feature_type, single, concatenated, seq_files):
    if not single and not concatenated:
        return

    same_features, parts = annotations_step.extract_shared_features(feature_type)
    # parts is dict ((seq_ident, gene) -> seq part))

    if single:
        for gene in same_features:
            al_step = _create_al_step(created_steps, steps, f'{feature_type}_{gene}', annotations_step, step_data)
            seq_files.append(_add_sequences(
                al_step, 'gene', sequences, [(seq_ident, parts[(seq_ident, gene, None)]) for seq_ident in sequences]))

    if concatenated:
        al_step = _create_al_step(created_steps, steps, f'{feature_type}_concatenated', annotations_step, step_data)
        same_features = sorted(same_features)  # ToDo: order?!
        seq_files.append(_add_sequences(
            al_step, 'genes', sequences,
            [(seq_ident,
              ''.join(parts[(seq_ident, gene)] for gene in same_features),
              _lengths_2_partitions(((gene, parts[(seq_ident, gene)]) for gene in same_features)))
             for seq_ident in sequences]))


def create_clustal_data(step_data, annotations_step, cache, alignments, run):
    alignments = set(alignments)
    assert all(a in ('w', 'gs', 'gc', 'cs', 'cc') for a in alignments)

    # List of dicts with attrs: filename, short, max_seq_length
    # This data is used to optimize calculation
    seq_files = []
    only_one_seq = sum(int(a in alignments) for a in ('w', 'gc', 'cc')) == 1 and \
        'gs' not in alignments and 'cs' not in alignments
    steps = None if only_one_seq else AlignmentsStep(annotations_step.project, step_data, remove_data=True)
    created_steps = []

    sequences = sorted(annotations_step.all_sequences())
    _feature_sequences(
        created_steps, step_data, steps, annotations_step, sequences,
        'gene', 'gs' in alignments, 'gc' in alignments, seq_files)
    _feature_sequences(
        created_steps, step_data, steps, annotations_step, sequences,
        'CDS', 'cs' in alignments, 'cc' in alignments, seq_files)
    if 'w' in alignments:
        al_step = _create_al_step(created_steps, steps, 'whole', annotations_step, step_data)
        # ToDo: kako anotirati cijelu sekvencu? Prema gene, CDS. Treba paziti da su sortirani tako da znamo rupe
        seq_files.append(_add_sequences(
            al_step, 'whole', sequences,
            [(seq_ident, annotations_step.get_sequence(seq_ident), None) for seq_ident in sequences]))

    assert created_steps
    base_step = steps or created_steps[0]

    # Store files desc
    files_to_zip = [d['filename'] for d in seq_files]  # files to zip
    # Remove step directory from files since run script is called from step directory
    for d in seq_files:
        d['filename'] = base_step.strip_step_dir(d['filename'])
    files_to_zip.append(base_step.step_file('finish.yml'))
    write_yaml(seq_files, files_to_zip[-1])

    # Stores description.yml
    base_step.save(completed=run)

    if run:
        run_module_script(run_clustal_omega, base_step)
    else:
        set_run_instructions(run_clustal_omega, base_step, files_to_zip, _instructions)
    #
    return base_step


def finish_clustal_data(step_obj, cache):
    output_f = step_obj.step_file('output.zip')
    if not os.path.isfile(output_f):
        raise ZCItoolsValueError('No calculation output file output.zip!')

    # Check are needed files in zip, not something strange
    files = set(d['filename'].replace('sequences.fa', 'alignment.phy')
                for d in read_yaml(step_obj.step_file('finish.yml')))
    # ToDo: possible problems with file separator
    sets_equal(files, set(list_zip_files(output_f)), 'file')  # raise exception if data is not good

    # Unzip data
    unzip_file(output_f, step_obj.directory)

    step_obj._check_data()
    step_obj.save(create=False)
