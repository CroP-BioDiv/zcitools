import os.path
from types import SimpleNamespace
from common_utils.file_utils import copy_file, write_str_in_file, write_fasta
from ..sequences.steps import SequencesStep
from ..annotations.steps import AnnotationsStep
from .utils import find_regions
from ..utils.features import Feature

# Orientates chloroplast sequence in standard way.
# Uses Fast-Plast method. Check
#  - source file orientate_plastome_v.2.0.pl
#    (https://github.com/mrmckain/Fast-Plast/blob/master/bin/orientate_plastome_v.2.0.pl)
#  - explanation https://github.com/mrmckain/Fast-Plast/issues/22
# Consitent with Wikipedia image:
#  - https://en.wikipedia.org/wiki/Chloroplast_DNA#/media/File:Plastomap_of_Arabidopsis_thaliana.svg

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


def orientate_chloroplast(command_obj, cmd_args, step_data, annotation_step, common_db):
    bad = to_repair = None
    good_files = []
    # ToDo: is this good? What if common_db_identifier() is a list?
    common_db_annot = common_db.get_relative_db('..', annotation_step.common_db_identifier()[-1])

    for seq_ident, seq in annotation_step._iterate_records():
        an_file = annotation_step.get_sequence_filename(seq_ident)

        partition = find_regions(seq_ident, seq)
        if not partition:
            # ToDo: what to do? Check is orientation good with list of genes. If not raise 
            print('  WHAT TO DO {seq_ident}?')
            if not bad:
                bad = _create_step(AnnotationsStep, 'bad', command_obj, step_data, cmd_args)
            _copy_annotation(bad, an_file)
            continue

        # Count gene orintation
        l_seq = len(seq)
        in_parts = partition.put_features_in_parts(
            Feature(l_seq, feature=f) for f in seq.features if f.type == 'gene' and f.location)

        lsc_count = sum(f.feature.strand if any(x in f.name for x in ('rpl', 'rps')) else 0
                        for f in in_parts.get('lsc', []))
        ssc_count = sum(f.feature.strand for f in in_parts.get('ssc', []))
        ir_count = sum(f.feature.strand if 'rrn' in f.name else 0 for f in in_parts.get('ira', []))

        if (lsc_count <= 0) and (ssc_count <= 0) and (ir_count >= 0):  # No change
            good_files.append(an_file)
            common_db.set_record(seq_ident, an_file)
            common_db_annot.set_record(seq_ident, an_file)
            continue

        parts = partition.extract(seq)  # dict name -> Seq object
        if ssc_count > 0:
            parts['ssc'] = parts['ssc'].reverse_complement()
        if lsc_count > 0:
            assert False, f'REVERT LSC {seq_ident}'
        if ir_count < 0:
            assert False, f'REVERT IRs {seq_ident}'

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
        common_db.set_record(seq_ident, fa_filename)

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
