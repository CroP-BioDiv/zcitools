import shutil
import os.path
from common_utils.file_utils import write_fasta
from .utils import orient_chloroplast_parts_by_data, orient_by_trnH_GUG_by_data
from ..utils.import_methods import import_bio_seq_io
from zci_bio.sequences.steps import SequencesStep


def _copy_from_origin(step, annotation_step, seq_ident):
    an_filename = annotation_step.get_sequence_filename(seq_ident)
    shutil.copy(an_filename, step.directory)
    step.add_sequence_file(os.path.basename(an_filename))


def fix_by_parts(step_data, analyse_step, keep_offset, sequences_db, annotations_db):
    step = SequencesStep(analyse_step.project, step_data, remove_data=True)
    analyse_step.propagate_step_name_prefix(step)
    annotation_step = analyse_step.project.find_previous_step_of_type(analyse_step, 'annotations')

    #
    for row in analyse_step.rows_as_dicts():
        seq_ident = row['AccesionNumber']
        if not (starts := row['Part starts']):
            print(f"Warning: sequence {seq_ident} doesn't have parts!")
            _copy_from_origin(step, annotation_step, seq_ident)
            continue

        starts = [int(f.strip()) for f in starts.split(',')]
        lsc_offset = starts[0]
        orientation = row['Orientation']

        # If there is no need to orient parts or to do offseting, than copy some previous record
        if not orientation and (abs(lsc_offset) <= keep_offset):
            _copy_from_origin(step, annotation_step, seq_ident)
            continue

        # Check is sequence already in the CommonDB
        new_seq_ident = step.seq_ident_of_our_change(seq_ident, 'p')
        if sequences_db and (f := sequences_db.get_record(new_seq_ident, step.directory, info=True)):
            step.add_sequence_file(os.path.basename(f))
            continue

        seq_rec = annotation_step.get_sequence_record(seq_ident)
        if orientation:  # Orientate parts
            new_seq = orient_chloroplast_parts_by_data(seq_rec, orientation, starts=starts)

            # Store fasta file, in step and sequences DB
            fa_filename = step.step_file(f'{new_seq_ident}.fa')
            write_fasta(fa_filename, [(new_seq_ident, str(new_seq.seq))])
            step.add_sequence_file(os.path.basename(fa_filename))
            if sequences_db:
                sequences_db.set_record(new_seq_ident, fa_filename)

        else:  # Offset sequence
            # Note: it is not needed to make offset if orientation was changed,
            # since parts concatenation orients the sequence
            assert abs(lsc_offset) > keep_offset, (lsc_offset, keep_offset)
            new_seq = seq_rec[lsc_offset:] + seq_rec[:lsc_offset]  # Note: this keeps features
            _store_genbank(step, new_seq_ident, new_seq, seq_rec, sequences_db, annotations_db)

    #
    step.save()
    return step


def fix_by_trnH_GUG(step_data, analyse_step, keep_offset, sequences_db, annotations_db):
    step = SequencesStep(analyse_step.project, step_data, remove_data=True)
    analyse_step.propagate_step_name_prefix(step)
    annotation_step = analyse_step.project.find_previous_step_of_type(analyse_step, 'annotations')

    for row in analyse_step.rows_as_dicts():
        seq_ident = row['AccesionNumber']
        if (trnh_gug := row['trnH-GUG']) is None:
            print(f"Warning: sequence {seq_ident} doesn't have trnH-GUG gene!")
            # _copy_from_origin(step, annotation_step, seq_ident)  # ???
            continue

        fields = trnh_gug.split(':')
        offset = int(fields[0])
        zero_reverse = (len(fields) > 1)
        orientation = row['Orientation']
        if starts := row['Part starts']:
            starts = [int(f.strip()) for f in starts.split(',')]
            # Nothing to do!
            if (starts[0] + offset) <= keep_offset and not orientation:
                _copy_from_origin(step, annotation_step, seq_ident)
                continue

        # Check is sequence already in the CommonDB
        new_seq_ident = step.seq_ident_of_our_change(seq_ident, 't')
        if sequences_db and (f := sequences_db.get_record(new_seq_ident, step.directory, info=True)):
            step.add_sequence_file(os.path.basename(f))
            continue

        seq_rec = annotation_step.get_sequence_record(seq_ident)
        new_seq = orient_by_trnH_GUG_by_data(
            seq_rec, offset, zero_reverse, orientation, start=starts, keep_offset=keep_offset)

        _store_genbank(step, new_seq_ident, new_seq, seq_rec, sequences_db, annotations_db)

    #
    step.save()
    return step


def _store_genbank(step, new_seq_ident, seq_rec, copy_from_rec, sequences_db, annotations_db):
    # Store GenBank file, in step and both sequences and annotations DBs
    _copy_sequence_annotations(copy_from_rec, seq_rec)
    # Set step file
    loc_filename = f'{new_seq_ident}.gb'
    gb_filename = step.step_file(loc_filename)
    import_bio_seq_io().write([seq_rec], gb_filename, 'genbank')
    step.add_sequence_file(loc_filename)
    # Set common DB files
    if sequences_db:
        sequences_db.set_record(new_seq_ident, gb_filename)
    if annotations_db:
        annotations_db.set_record(new_seq_ident, gb_filename)


def _copy_sequence_annotations(from_rec, to_seq):
    f_ann = from_rec.annotations
    t_ann = to_seq.annotations
    for k in ('molecule_type', 'topology', 'accessions', 'organism', 'taxonomy'):
        if v := f_ann.get(k):
            t_ann[k] = v
