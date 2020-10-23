import shutil
import os.path
from common_utils.file_utils import write_fasta
from .utils import create_chloroplast_partition_all
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
        starts = row['Took part starts'] or row['Part starts']
        if not starts:
            print(f"Warning: sequence {seq_ident} doesn't have starts!")
            _copy_from_origin(step, annotation_step, seq_ident)
            continue

        starts = [int(f.strip()) for f in starts.split(',')]
        offset = starts[0]
        l_seq = row['Length']
        orientation = row['Orientation']

        if (abs(offset) <= omit_offset) and not orientation:
            _copy_from_origin(step, annotation_step, seq_ident)
            continue

        seq = new_seq = annotation_step.get_sequence_record(seq_ident)
        new_seq_ident = step.seq_ident_of_our_change(seq_ident, 'p')

        if common_db and (f := common_db.get_record(new_seq_ident, step.directory, info=True)):
            step.add_sequence_file(os.path.basename(f))
            continue

        if orientation:  # Orientate parts
            partition = create_chloroplast_partition_all(l_seq, starts)
            parts = partition.extract(seq)  # dict name -> Seq object
            if 'lsc' in orientation:  # LSC
                parts['lsc'] = parts['lsc'].reverse_complement()
            if 'ssc' in orientation:  # SSC
                parts['ssc'] = parts['ssc'].reverse_complement()
            if 'ira' in orientation:  # IRs
                parts['ira'] = parts['ira'].reverse_complement()
                parts['irb'] = parts['irb'].reverse_complement()

            new_seq = parts['lsc'] + parts['ira'] + parts['ssc'] + parts['irb']
            assert len(seq.seq) == len(new_seq.seq), \
                (seq_ident, len(seq.seq), len(new_seq.seq), starts,
                    [(n, len(p)) for n, p in parts.items()],
                    [(n, len(p)) for n, p in partition.extract(seq).items()])

        # Offset sequence
        # Note: it is not needed to make offset if orientation was changed,
        # since parts concatenation orients the sequence
        elif offset:
            new_seq = new_seq[offset:] + new_seq[0:offset]

        # Store file
        _store_fasta(step, new_seq_ident, new_seq, common_db)

    #
    step.save()
    return step
