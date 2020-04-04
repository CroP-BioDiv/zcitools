import os.path
import itertools
from common_utils.file_utils import ensure_directory, write_yaml, write_fasta, \
    run_module_script, set_run_instructions
from common_utils.exceptions import ZCItoolsValueError
from common_utils.value_data_types import rows_2_excel
from ..utils.features import Feature
from ..utils.import_methods import import_bio_align_io
from .steps import ChloroplastOrientateStep
from .utils import find_chloroplast_partition, find_referent_genome
from . import run_orientate

# -------------------------------------------------------------------------
# Check chloroplast seqeunce parts orientation
# -------------------------------------------------------------------------
_part_names = ('lsc', 'ira', 'ssc')
_plus_minus = ('plus', 'minus')
_instructions = """
"""


def orientate_chloroplast_start(step_data, annotation_step, params):
    # Find referent genome
    # For each sequence, different than referent, directory is created named <seq_ident>.
    # It contains files:
    #  - {lsc|ira|ss}_{plus|minus}.fa       : input alignment files, contain 2 sequences.
    #  - align_{lsc|ira|ss}_{plus|minus}.fa : result alignment files.
    seq_idents = annotation_step.all_sequences()  # set
    ref_ident = find_referent_genome(seq_idents, params.referent_genome)
    #
    length = params.length_to_check
    step = annotation_step.project.new_step(ChloroplastOrientateStep, step_data, remove_data=False)
    sequence_data = step.get_type_description_elem('sequence_data', default=dict())
    #
    seq_rec = annotation_step.get_sequence_record(ref_ident)
    partition = find_chloroplast_partition(seq_rec)
    ref_parts = [str(partition.get_part_by_name(n).extract(seq_rec).seq)[:length] for n in _part_names]
    files_to_zip = []
    align_files = []

    #
    all_versions = ('plus', 'minus', 'plus_c', 'minus_c') if params.complement else ('plus', 'minus')
    for seq_ident in sorted(seq_idents):
        seq_rec = None
        if seq_ident not in sequence_data:
            seq_rec = annotation_step.get_sequence_record(seq_ident)
            partition = find_chloroplast_partition(seq_rec)

            # Count gene orientation
            l_seq = len(seq_rec)
            in_parts = partition.put_features_in_parts(
                Feature(l_seq, feature=f) for f in seq_rec.features if f.type == 'gene')

            lsc_count = sum(f.feature.strand if any(x in f.name for x in ('rpl', 'rps')) else 0
                            for f in in_parts.get('lsc', []))
            ssc_count = sum(f.feature.strand for f in in_parts.get('ssc', []))
            ira_count = sum(f.feature.strand if 'rrn' in f.name else 0 for f in in_parts.get('ira', []))

            sequence_data[seq_ident] = dict(
                length=len(seq_rec),
                lsc=(lsc_count <= 0), lsc_count=lsc_count, lsc_length=len(partition.get_part_by_name('lsc')),
                ssc=(ssc_count <= 0), ssc_count=ssc_count, ssc_length=len(partition.get_part_by_name('ssc')),
                ira=(ira_count >= 0), ira_count=ira_count, ira_length=len(partition.get_part_by_name('ira')))

        if all(all(step.is_file(seq_ident, f'align_{n}_{v}.fa') for v in all_versions) for n in _part_names):
            continue
        #
        if seq_rec is None:
            seq_rec = annotation_step.get_sequence_record(seq_ident)
            partition = find_chloroplast_partition(seq_rec)
        for n, ref_p in zip(_part_names, ref_parts):
            # Find missing output files
            _num = len(align_files)
            for x in all_versions:
                if not step.is_file(seq_ident, f'align_{n}_{x}.fa'):
                    files_to_zip.append(step.step_file(seq_ident, f'{n}_{x}.fa'))
                    align_files.append((seq_ident, n, x))
            if _num == len(align_files):
                continue

            # Store input files
            if all(step.is_file(seq_ident, f'align_{n}_{v}.fa') for v in all_versions):
                continue
            ensure_directory(step.step_file(seq_ident))
            part_s = partition.get_part_by_name(n).extract(seq_rec)

            f_p = step.step_file(seq_ident, f'{n}_plus.fa')
            f_p_c = step.step_file(seq_ident, f'{n}_plus_c.fa')
            if not os.path.isfile(f_p):
                write_fasta(f_p, [(ref_ident, ref_p), (seq_ident, str(part_s.seq)[:length])])
            if not os.path.isfile(f_p_c):
                write_fasta(f_p_c, [(ref_ident, ref_p),
                                    (seq_ident, str(part_s.reverse_complement().seq)[:(-length-1):-1])])

            f_m = step.step_file(seq_ident, f'{n}_minus.fa')
            f_m_c = step.step_file(seq_ident, f'{n}_minus_c.fa')
            if not os.path.isfile(f_m):
                write_fasta(f_m, [(ref_ident, ref_p), (seq_ident, str(part_s.reverse_complement().seq)[:length])])
            if not os.path.isfile(f_m_c):
                write_fasta(f_m_c, [(ref_ident, ref_p), (seq_ident, str(part_s.seq)[:(-length-1):-1])])

    #
    output_file = f"{params.output_file_prefix}_{length}{'_c' if params.complement else ''}.xlsx"
    data = dict(sequence_data=sequence_data, check_length=length, output_file=output_file, complement=params.complement)
    if align_files:
        # Store finish.yml
        finish_f = step.step_file('finish.yml')
        write_yaml(dict(align_files=align_files), finish_f)

        run = True  # ToDo: ...
        step.save(data, completed=False)
        if run:
            run_module_script(run_orientate, step)
            orientate_chloroplast_finish(step)  # , common_db, calc_seq_idents=calc_seq_idents)
        else:
            files_to_zip.append(finish_f)
            set_run_instructions(run_orientate, step, files_to_zip, _instructions)
    #
    elif params.force_parse:
        step.save(data)
        orientate_chloroplast_finish(step)  # , common_db, calc_seq_idents=calc_seq_idents)
    #
    else:
        step.save(data, completed=False)

    return step


def _start_length(seq_str):
    for i, s in enumerate(seq_str):
        if s != '-':
            return i
    return len(seq_str)


def _start_diff(align):
    seq0_str = str(align[0].seq)
    if seq0_str[0] == '-':
        return -_start_length(seq0_str)
    return _start_length(str(align[1].seq))


def _add_in_row(row, plus_f, minus_f):
    AlignIO = import_bio_align_io()
    a_p = AlignIO.read(plus_f, 'fasta')
    a_m = AlignIO.read(minus_f, 'fasta')
    l_p = a_p.get_alignment_length()
    l_m = a_m.get_alignment_length()
    is_plus = (l_p < l_m)
    row.append(f"{'Plus' if is_plus else 'Minus'} {l_p} / {l_m} ({_start_diff(a_p if is_plus else a_m)})")
    return is_plus


def orientate_chloroplast_finish(step_obj):
    # ToDo: sto bi sve htjeli znati? Duljinu sekvence, duljine djelova?
    type_desciption = step_obj.get_type_description()
    sequence_data = type_desciption['sequence_data']
    complement = type_desciption['complement']
    #
    rows = []
    for seq_ident in sorted(step_obj.step_subdirectories()):
        seq_data = sequence_data[seq_ident]
        row = [seq_ident, seq_data['length']]
        rows.append(row)
        all_ok = True
        for n in _part_names:
            row.append(seq_data[n + '_length'])
            row.append(f"{'Plus' if seq_data[n] else 'Minus'} / {seq_data[n + '_count']}")
            #
            is_plus = _add_in_row(
                row,
                step_obj.step_file(seq_ident, f'align_{n}_plus.fa'),
                step_obj.step_file(seq_ident, f'align_{n}_minus.fa'))
            if complement:
                _add_in_row(
                    row,
                    step_obj.step_file(seq_ident, f'align_{n}_plus_c.fa'),
                    step_obj.step_file(seq_ident, f'align_{n}_minus_c.fa'))
            #
            if seq_data[n] != is_plus:
                all_ok = False
        #
        row.append('OK' if all_ok else 'No')

    # Create table
    columns = ['Seq', 'Length']  # + list(_part_names)
    for n in _part_names:
        columns.append(f"{n} length")
        columns.append(f"NP {n}")
        columns.append(f"Ref {n}")
        if complement:
            columns.append(f"Complement {n}")

    columns.append('OK')
    rows_2_excel(type_desciption['output_file'], columns, rows)


# -------------------------------------------------------------------------
# Orientate chloroplast sequences
# -------------------------------------------------------------------------
# Orientates chloroplast sequence in standard way.
# Uses Fast-Plast method. Check
#  - source file orientate_plastome_v.2.0.pl
#    (https://github.com/mrmckain/Fast-Plast/blob/master/bin/orientate_plastome_v.2.0.pl)
#  - explanation https://github.com/mrmckain/Fast-Plast/issues/22
# Consitent with Wikipedia image:
#  - https://en.wikipedia.org/wiki/Chloroplast_DNA#/media/File:Plastomap_of_Arabidopsis_thaliana.svg
"""
CHLOROPLAST_ORIENTED_DB_NAME = 'chloroplast_oriented'


def _copy_annotation(to_step, from_file):
    base_file = os.path.basename(from_file)
    step_f = to_step.step_file(base_file)
    copy_file(from_file, step_f)
    to_step.get_sequence_filename(step_f)


def _create_step(step_cls, description, command_obj, step_data, cmd_args):
    project = command_obj.project
    step_data = dict(step_data)  # Copy step data
    step_data['step_name'] = project.new_step_name(
        command_obj, SimpleNamespace(step_num=cmd_args.step_num, step_description=description))
    return step_cls(project, step_data, remove_data=True)


def find_chloroplast_partition_force_calc(step, seq_ident, seq_rec):
    # Check are IRs present
    partition = find_chloroplast_partition(seq_rec)
    if not partition:
        # Try to find IRs with MUMmer
        # Note: calculation data is set in annotations step
        if calculate_and_add_irs_to_seq_rec(step, seq_ident, seq_rec):
            partition = find_chloroplast_partition(seq_rec)
    return partition


def orientate_chloroplast(command_obj, cmd_args, step_data, annotation_step, common_db):
    # common_db object's directory is set to chloroplast_oriented/sequences
    # common_db_annot's directory is set to annotation's directory (e.g. chloroplast_oriented/GeSeq)
    common_db_annot = common_db.get_relative_db('..', annotation_step.common_db_identifier()[-1])
    bad = to_repair = None
    good_files = []

    for seq_ident, seq in annotation_step._iterate_records():
        an_file = annotation_step.get_sequence_filename(seq_ident)

        partition = find_chloroplast_partition_force_calc(annotation_step, seq_ident, seq)
        if not partition:
            print(f'ERROR: no IRs for sequence {seq_ident}!')
            if not bad:
                bad = _create_step(AnnotationsStep, 'bad', command_obj, step_data, cmd_args)
            _copy_annotation(bad, an_file)
            continue

        # Fix sequence features
        seq.features = [f for f in seq.features if f.location]

        # Count gene orintation
        l_seq = len(seq)
        in_parts = partition.put_features_in_parts(Feature(l_seq, feature=f) for f in seq.features if f.type == 'gene')

        lsc_count = sum(f.feature.strand if any(x in f.name for x in ('rpl', 'rps')) else 0
                        for f in in_parts.get('lsc', []))
        ssc_count = sum(f.feature.strand for f in in_parts.get('ssc', []))
        ir_count = sum(f.feature.strand if 'rrn' in f.name else 0 for f in in_parts.get('ira', []))

        # ToDo: ako je pozicija IR
        lsc = partition.get_part_by_name('lsc')
        irb = partition.get_part_by_name('irb')
        # irs_good_positioned = min(l_seq - lsc.real_start, lsc.real_start) < 30 and \
        #     min(l_seq - irb.real_end, irb.real_end) < 30

        if (lsc_count <= 0) and (ssc_count <= 0) and (ir_count >= 0):  # and irs_good_positioned:  # No change
            good_files.append(an_file)
            common_db.set_record(seq_ident, an_file)
            common_db_annot.set_record(seq_ident, an_file)
            continue

        parts = partition.extract(seq)  # dict name -> Seq object
        if ssc_count > 0:
            parts['ssc'] = parts['ssc'].reverse_complement()
        if lsc_count > 0:
            parts['lsc'] = parts['lsc'].reverse_complement()
        if ir_count < 0:
            print(f'  REVERT IRs: WHAT TO DO {seq_ident}?')
            continue

        if not to_repair:
            to_repair = _create_step(SequencesStep, 'to_repair', command_obj, step_data, cmd_args)

        new_seq = parts['lsc'] + parts['ira'] + parts['ssc'] + parts['irb']
        assert len(seq.seq) == len(new_seq.seq), \
            (seq_ident, len(seq.seq), len(new_seq.seq),
                [(n, len(p)) for n, p in parts.items()],
                [(n, len(p)) for n, p in partition.extract(seq).items()])
        base_f = os.path.basename(an_file)
        base_f = base_f.split('.')[0] + '.fa'
        fa_filename = to_repair.step_file(base_f)
        write_fasta(fa_filename, [(seq_ident, str(new_seq.seq))])
        to_repair.add_sequence_file(base_f)
        # ToDo: force stavljanja u common_db. Brisati u common_db_annot
        common_db.set_record(seq_ident, fa_filename)
        # common_db_annot.remove_record(seq_ident)

    #
    mess = []
    ret_steps = []
    if bad:
        ret_steps.append(bad)
        bad.save()
        mess.append(f'Error: sequences in step {bad.directory} is not possible to repair!')

    if not to_repair:
        if not bad:
            mess.append(f'All sequences are oriented good!')
    else:
        if good_files:
            good = _create_step(AnnotationsStep, 'good', command_obj, step_data, cmd_args)
            ret_steps.append(good)
            for an_file in good_files:
                _copy_annotation(good, an_file)
            good.save()
            mess.append(f'Good oriented files are set in step {good.directory}!')

        to_repair.save()
        ret_steps.append(to_repair)
        mess.append(f'Sequences to orient are stored in step {to_repair.directory}!')

    if mess:
        print()
        print('\n'.join(mess))
        print()

    if to_repair:
        # Run annotation with same command that was used for input annotation
        annotation_command = cmd_args.annotation_command or annotation_step.get_command()
        args = SimpleNamespace(step=to_repair.directory, sequence_db=CHLOROPLAST_ORIENTED_DB_NAME,
                               step_num=None, step_description=None)
        annotation_step.project._run_command(annotation_command, args)

    return ret_steps
"""
