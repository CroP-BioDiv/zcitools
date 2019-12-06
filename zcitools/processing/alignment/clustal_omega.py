import os.path
from . import run_clustal_omega
from zcitools.steps.alignments import AlignmentStep, AlignmentsStep
from zcitools.utils.file_utils import zip_files, unzip_file, list_zip_files, write_yaml, read_yaml, \
    copy_file, write_str_in_file, write_fasta
from zcitools.utils.helpers import sets_equal

_instructions = """
Steps:
 - copy file calculate.zip onto server
 - unzip it
 - change directory to {step_name}
 - run script: {script_name}
 - copy file output.zip back into project's step directory {step_name}
 - run zcit command: zcit.py finish {step_name}

Notes:
 - Clustal Omega executable (clustalo) should be on the PATH or
   environment variable CLUSTAL_OMEGA_EXE should point to it.
"""

# Note: global cache is not used!


def _add_sequences(step, short, sequences, sequence_data):
    seq_file = step.step_file('sequences.fa')
    write_fasta(seq_file, sequence_data)
    step.set_sequences(sequences)
    step.save(needs_editing=True)
    return dict(filename=seq_file, short=short, max_seq_length=max(len(s) for _, s in sequence_data))


def _feature_sequences(step, annotations_step, sequences, feature_type, single, concatenated, seq_files):
    if not single and not concatenated:
        return

    same_features, parts = annotations_step.extract_shared_features(feature_type)
    # parts is dict ((seq_ident, gene) -> seq part))

    if single:
        for gene in same_features:
            substep = step.create_substep(AlignmentStep, f'{feature_type}_{gene}')
            seq_files.append(_add_sequences(
                substep, True, sequences, [(seq_ident, parts[(seq_ident, gene)]) for seq_ident in sequences]))

    if concatenated:
        substep = step.create_substep(AlignmentStep, f'{feature_type}_concatenated')
        same_features = sorted(same_features)  # ToDo: order?!
        seq_files.append(_add_sequences(
            substep, False, sequences,
            [(seq_ident, ''.join(parts[(seq_ident, gene)] for gene in same_features)) for seq_ident in sequences]))


def create_clustal_data(step_data, annotations_step, cache, alignments, run):
    alignments = set(alignments)
    assert all(a in ('w', 'gs', 'gc', 'cs', 'cc') for a in alignments)

    sequences = sorted(annotations_step.all_sequences())

    # List of dicts with attrs: filename, short, max_seq_length
    # This data is used to optimize calculation
    seq_files = []

    if sum(int(a in alignments) for a in ('w', 'gc', 'cc')) == 1 and \
            all(a not in alignments for a in ('w', 'gc', 'cc')):
        # Only one sequence is aligned
        step = AlignmentStep(step_data, remove_data=True)
        assert False, 'Not implemented!!!'

        if 'w' in alignments:
            sequence_data = [(seq_ident, annotations_step.get_sequence(seq_ident)) for seq_ident in sequences]
        else:
            feature_type = 'gene' if 'gc' in alignments else 'CDS'
            same_features, parts = annotations_step.extract_shared_features(feature_type)
            same_features = sorted(same_features)  # ToDo: order?!
            sequence_data = [(seq_ident, ''.join(parts[(seq_ident, gene)] for gene in same_features))
                             for seq_ident in sequences]

        seq_file = step.step_file('sequences.fa')
        seq_files.append(dict(filename=seq_file, short=False, max_seq_length=max(len(s) for _, s in sequence_data)))
        write_fasta(seq_file, sequence_data)
        #
        step.set_sequences(sequences)

    else:
        # Lot of sequence are aligned
        step = AlignmentsStep(step_data, remove_data=True)
        _feature_sequences(step, annotations_step, sequences, 'gene', 'gs' in alignments, 'gc' in alignments, seq_files)
        _feature_sequences(step, annotations_step, sequences, 'CDS', 'cs' in alignments, 'cc' in alignments, seq_files)
        if 'w' in alignments:
            substep = step.create_substep(AlignmentStep, 'whole')
            seq_files.append(_add_sequences(
                substep, False, sequences, [annotations_step.get_sequence(seq_ident) for seq_ident in sequences]))

    # Store files desc
    files_to_zip = [d['filename'] for d in seq_files]  # files to zip
    # Remove step directory from files since run script is called from step directory
    for d in seq_files:
        d['filename'] = step.strip_step_dir(d['filename'])
    finish_f = step.step_file('finish.yml')
    write_yaml(seq_files, finish_f)

    # Stores description.yml
    step.save(needs_editing=not run)

    if run:
        project_dir = os.getcwd()
        os.chdir(step.directory)
        run_clustal_omega.run(locale=True)
        os.chdir(project_dir)
    else:
        # Copy run script
        script_name = os.path.basename(run_clustal_omega.__file__)
        run_f = step.step_file(script_name)
        copy_file(run_clustal_omega.__file__, run_f)

        # Instructions
        inst_f = step.step_file('INSTRUCTIONS.txt')
        write_str_in_file(inst_f, _instructions.format(step_name=step_data['step_name'], script_name=script_name))

        # Zip needed files
        zip_files(step.step_file('calculate.zip'), files_to_zip + [inst_f, run_f, finish_f])

    #
    return step


def finish_clustal_data(step_obj, cache):
    output_f = step_obj.step_file('output.zip')
    # Check are needed files in zip, not something strange
    files = set(d['filename'].replace('sequences.fa', 'alignment.phy')
                for d in read_yaml(step_obj.step_file('finish.yml')))
    # ToDo: possible problems with file separator
    sets_equal(files, set(list_zip_files(output_f)), 'file')  # raise exception if data is not good

    # Unzip data
    unzip_file(step_obj.step_file('output.zip'), step_obj.directory)

    step_obj._check_data()
    step_obj.save(create=False)
