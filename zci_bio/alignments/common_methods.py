import os.path
import importlib
from .steps import AlignmentStep, AlignmentsStep
from ..utils.import_methods import import_bio_seq_io
from common_utils.file_utils import write_fasta, write_yaml, run_module_script, set_run_instructions, \
    unzip_file, list_zip_files, read_yaml
from common_utils.misc import sets_equal
from common_utils.exceptions import ZCItoolsValueError

"""
Alignment step contains one (type AlignmentStep) or more alignments (type AlignmentsStep).

For each alignment (step/substep), sequences to align are stored in file 'sequences.fa'.

For main step, file 'finish.yml' contains characteristics of all sets to align,
so that alignment script can optimize hot to run alignment program.
For each set to align, data stored are:
 - filename       : filename of data to align (<dir>/sequences.fa)
 - short          : should set be threated as 'short' task. It is only a hint, not request to the script.
 - max_seq_length : maximal length of sequence to align.
 - namelength     : max length of sequence names. Needed for phyllip format.
"""

_align_programs = dict(
    clustal_omega=dict(run_module='run_clustal_omega',
                       instructions="""
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
"""),
    mafft=dict(run_module='run_mafft',
               instructions="""
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
"""),
    muscle=dict(run_module='run_muscle',
                instructions="""
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
"""),
)
ALIGN_PROGRAMS = sorted(_align_programs.keys())


def run_alignment_program(alignment_program, base_step, seq_files, run):
    # Check alignment program
    if not (data := _align_programs.get(alignment_program)):
        raise ZCItoolsValueError(f'Alignment program {alignment_program} is not recognized!!!')

    # Store files desc
    files_to_zip = [d['filename'] for d in seq_files]  # files to zip
    # Remove step directory from files since run script is called from step directory
    for d in seq_files:
        d['filename'] = base_step.strip_step_dir(d['filename'])
    files_to_zip.append(base_step.step_file('finish.yml'))
    write_yaml(seq_files, files_to_zip[-1])

    # Stores description.yml
    base_step.save(completed=run)

    # Run or store script
    run_module = importlib.import_module('.' + data['run_module'], package='zci_bio.alignments')
    if run:
        run_module_script(run_module, base_step)
    else:
        set_run_instructions(run_module, base_step, files_to_zip, data['instructions'])


def add_sequences(al_step, seq_type, sequence_data):
    seq_file = al_step.step_file('sequences.fa')
    write_fasta(seq_file, ((i, s) for i, s, _ in sequence_data))
    al_step.set_sequences(i for i, _, _ in sequence_data)
    al_step.seq_sequence_type(seq_type)
    al_step.save(completed=False)
    # Annotations
    for seq_ident, seq, partition in sequence_data:
        if partition:
            al_step.store_partition(seq_ident, partition)
    #
    return dict(filename=seq_file,
                short=al_step.is_short(),
                namelength=max(len(i) for i, _, _ in sequence_data),
                max_seq_length=max(len(s) for _, s, _ in sequence_data))


def _lengths_2_features(data):
    # data is list of strings (ident, string)
    ret = []
    end = 0
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
            seq_files.append(add_sequences(
                al_step, 'gene', [(seq_ident, parts[(seq_ident, gene, None)]) for seq_ident in sequences]))

    if concatenated:
        al_step = _create_al_step(created_steps, steps, f'{feature_type}_concatenated', annotations_step, step_data)
        same_features = sorted(same_features)  # ToDo: order?!
        seq_files.append(add_sequences(
            al_step, 'genes',
            [(seq_ident,
              ''.join(parts[(seq_ident, gene)] for gene in same_features),
              _lengths_2_features(((gene, parts[(seq_ident, gene)]) for gene in same_features)))
             for seq_ident in sequences]))


def create_alignment_data(step_data, annotations_step, alignments, whole_partition, run, alignment_program):
    # Check alignments
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
        seq_files.append(add_sequences(
            al_step, 'whole',
            [(seq_ident, annotations_step.get_sequence(seq_ident), None) for seq_ident in sequences]))

        # Store filtered features
        if whole_partition:
            SeqIO = import_bio_seq_io()
            for seq_ident in sequences:
                seq_rec = annotations_step.get_sequence_record(seq_ident)
                # Copy SeqRecord and filter features
                sr_copy = seq_rec[:]
                sr_copy.features = [f for f in sr_copy.features if f.type == whole_partition]
                SeqIO.write([sr_copy], al_step.step_file(f'{seq_ident}.gb'), 'genbank')

    assert created_steps
    base_step = steps or created_steps[0]
    annotations_step.propagate_step_name_prefix(base_step)
    run_alignment_program(alignment_program, base_step, seq_files, run)
    return base_step


def finish_alignment_data(step_obj, files):
    output_f = step_obj.step_file('output.zip')
    if not os.path.isfile(output_f):
        raise ZCItoolsValueError('No calculation output file output.zip!')

    # ToDo: possible problems with file separator
    if files:
        sets_equal(files, set(list_zip_files(output_f)), 'file')  # raise exception if data is not good

    # Unzip data
    unzip_file(output_f, step_obj.directory)

    step_obj._check_data()
    step_obj.save(create=False)
