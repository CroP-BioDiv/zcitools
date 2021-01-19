import os.path
import re
from common_utils.file_utils import link_file
from common_utils.exceptions import ZCItoolsValueError


def fetch_common_db_data(step_data, table_step, step_type, common_db):
    step = table_step.project.new_step_by_type(step_type, step_data, remove_data=True)

    for seq_ident in table_step.get_column_values_by_type('seq_ident'):
        f = common_db.get_record(seq_ident, step.directory, info=True)
        if not f:
            raise ZCItoolsValueError(f"There is not CommonDB record for seq ident {seq_ident}!")
        step.add_sequence_file(os.path.basename(f))

    step.save()
    return step


def create_subset(step_data, input_step, args):
    assert not args.analyses_subset or args.analyses_subset in ('sum', 'ge_seq', 'ncbi'), args.analyses_subset
    if args.analyses_with_irs and not args.analyses_subset:
        raise ZCItoolsValueError('Error: analyses subset was not set!')

    project = input_step.project
    step = project.new_step_by_type(input_step.step_data_type, step_data, remove_data=True)

    wo = set(args.without_seq_idents or [])
    wo_re = [re.compile(r) for r in (args.without_seq_idents_re or [])]
    seqs = set(args.seq_idents or [])
    seqs_re = [re.compile(r) for r in (args.seq_idents_re or [])]

    if args.analyses_with_irs:
        analyses_step = project.read_step(args.analyses_with_irs, check_data_type='table', no_check=True)
        if args.analyses_subset == 'sum':
            seqs.update(s for s, g, n in analyses_step.select(
                ['AccesionNumber', 'GeSeq part starts', 'NCBI part starts']) if g or n)
        elif args.analyses_subset == 'ge_seq':
            seqs.update(s for s, g in analyses_step.select(['AccesionNumber', 'GeSeq part starts']) if g)
        elif args.analyses_subset == 'ncbi':
            seqs.update(s for s, g in analyses_step.select(['AccesionNumber', 'NCBI part starts']) if g)

    seqs.difference_update(wo)
    use_all = all(not x for x in (seqs, seqs_re, wo, wo_re))

    for s in input_step.all_sequences():
        if use_all or \
           ((s in seqs or any(r.match(s) for r in seqs_re)) and not (s in wo or any(r.match(s) for r in wo_re))):
            input_filename = input_step.get_sequence_filename(s)
            loc_f = os.path.basename(input_filename)
            link_file(input_filename, os.path.join(step.directory, loc_f))
            step.add_sequence_file(loc_f)

    step.save()
    return step
