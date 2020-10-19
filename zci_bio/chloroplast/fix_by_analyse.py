import shutil
import os.path
from common_utils.file_utils import write_fasta
from .utils import find_chloroplast_partition, create_chloroplast_partition
from zci_bio.sequences.steps import SequencesStep


def _copy_from_origin(step, annotation_step, seq_ident):
    an_filename = annotation_step.get_sequence_filename(seq_ident)
    shutil.copy(an_filename, step.directory)
    step.add_sequence_file(os.path.basename(an_filename))


def _store_fasta(step, new_seq_ident, new_seq, common_db):
    fa_filename = step.step_file(f'{new_seq_ident}.fa')
    write_fasta(fa_filename, [(new_seq_ident, str(new_seq.seq))])
    step.add_sequence_file(os.path.basename(fa_filename))
    # ToDo: force stavljanja u common_db. Brisati u common_db_annot
    if common_db:
        common_db.set_record(new_seq_ident, fa_filename)


def fix_by_parts(step_data, analyse_step, common_db, omit_offset=10):
    project = analyse_step.project
    step = SequencesStep(project, step_data, remove_data=True)
    analyse_step.propagate_step_name_prefix(step)
    annotation_step = project.find_previous_step_of_type(analyse_step, 'annotations')

    #
    for row in analyse_step.rows_as_dicts():
        seq_ident = row['AccesionNumber']
        l_seq = row['Length']
        offset = row['Offset']
        orientation = row['Orientation']

        if (not_offset := (offset <= omit_offset or (l_seq - offset) <= omit_offset)) and \
           (not orientation):
            _copy_from_origin(step, annotation_step, seq_ident)
            continue

        seq = new_seq = annotation_step.get_sequence_record(seq_ident)
        new_seq_ident = step.seq_ident_of_our_change(seq_ident, 'p')

        if common_db and (f := common_db.get_record(new_seq_ident, step.directory, info=True)):
            step.add_sequence_file(os.path.basename(f))
            continue

        if orientation:  # Orientate parts
            if row['IRS took']:
                ir_l = row['IR']
                ira_end, irb_start = tuple(map(int, row['SSC ends'].split('-')))
                partition = create_chloroplast_partition(
                    l_seq, (ira_end - ir_l, ira_end), (irb_start, irb_start + ir_l), in_interval=True)
            else:
                partition = find_chloroplast_partition(seq)

            parts = partition.extract(seq)  # dict name -> Seq object
            if 'lsc' in orientation:  # LSC
                parts['lsc'] = parts['lsc'].reverse_complement()
            if 'ira' in orientation:  # IRs
                parts['ssc'] = parts['ssc'].reverse_complement()
            if 'ssc' in orientation:  # SSC
                print(f'  REVERT IRs: WHAT TO DO {seq_ident}?')
                continue

            new_seq = parts['lsc'] + parts['ira'] + parts['ssc'] + parts['irb']
            assert len(seq.seq) == len(new_seq.seq), \
                (seq_ident, len(seq.seq), len(new_seq.seq),
                    [(n, len(p)) for n, p in parts.items()],
                    [(n, len(p)) for n, p in partition.extract(seq).items()])

        if not not_offset:  # Offset sequence
            new_seq = new_seq[offset:] + new_seq[0:offset]

        # Store file
        _store_fasta(step, new_seq_ident, new_seq, common_db)

    #
    step.save()
    return step


def fix_by_trnh_gug(step_data, analyse_step, common_db, omit_offset=10):
    project = analyse_step.project
    step = SequencesStep(project, step_data, remove_data=True)
    analyse_step.propagate_step_name_prefix(step)
    annotation_step = project.find_previous_step_of_type(analyse_step, 'annotations')

    #
    for row in analyse_step.rows_as_dicts():
        seq_ident = row['AccesionNumber']
        l_seq = row['Length']
        offset = row['trnH-GUG']

        if not offset or offset <= omit_offset or (l_seq - offset) <= omit_offset:
            _copy_from_origin(step, annotation_step, seq_ident)
            continue

        seq = new_seq = annotation_step.get_sequence_record(seq_ident)
        new_seq_ident = step.seq_ident_of_our_change(seq_ident, 't')

        if common_db and (f := common_db.get_record(new_seq_ident, step.directory, info=True)):
            step.add_sequence_file(os.path.basename(f))
            continue

        new_seq = new_seq[offset:] + new_seq[0:offset]
        _store_fasta(step, new_seq_ident, new_seq, common_db)

    #
    step.save()
    return step
