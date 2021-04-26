import shutil
import os.path
from common_utils.file_utils import write_fasta
from .utils import orient_chloroplast_parts_by_data, orient_by_trnF_GAA_by_data  # , orient_by_trnH_GUG_by_data
from .constants import DEFAULT_KEEP_OFFSET
from ..utils.import_methods import import_bio_seq_io
from zci_bio.sequences.steps import SequencesStep


def _copy_from_origin(step, annotation_step, seq_ident):
    an_filename = annotation_step.get_sequence_filename(seq_ident)
    shutil.copy(an_filename, step.directory)
    step.add_sequence_file(os.path.basename(an_filename))


def fix_by_parts(step_data, analyses_step, subset, keep_offset, sequences_db):
    assert subset in ('all', 'sum', 'ge_seq', 'ncbi'), subset
    step = SequencesStep(analyse_step.project, step_data, remove_data=True)
    analyse_step.propagate_step_name_prefix(step)
    annotation_step = analyse_step.project.find_previous_step_of_type(analyse_step, 'annotations')
    fix_seq_prefix = 'n' if (subset == 'ncbi') else 'p'

    #
    for row in analyse_step.rows_as_dicts():
        seq_ident = row['AccesionNumber']
        if subset == 'ncbi':
            starts = row['NCBI part starts']
            orientation = row['NCBI orientation']
        else:
            starts = row['GeSeq part starts']
            orientation = row['GeSeq orientation']
            if not starts and subset != 'ge_seq':
                starts = row['NCBI part starts']
                orientation = row['NCBI orientation']

        if not starts:
            if subset == 'all':
                print(f"Warning: sequence {seq_ident} doesn't have parts!")
                _copy_from_origin(step, annotation_step, seq_ident)
            continue

        starts = [int(f.strip()) for f in starts.split(',')]
        lsc_offset = starts[0]

        # If there is no need to orient parts or to do offseting, than copy some previous record
        if not orientation and (abs(lsc_offset) <= keep_offset):
            _copy_from_origin(step, annotation_step, seq_ident)
            continue

        # Check is sequence already in the CommonDB
        new_seq_ident = step.seq_ident_of_our_change(seq_ident, fix_seq_prefix)
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
            new_seq = seq_rec[lsc_offset:] + seq_rec[:lsc_offset]
            _store_genbank(step, new_seq_ident, new_seq, seq_rec, sequences_db)

    #
    step.save()
    return step


def fix_by_trnF_GAA(step_data, analyse_step, subset, keep_offset, sequences_db):
    step = SequencesStep(analyse_step.project, step_data, remove_data=True)
    analyse_step.propagate_step_name_prefix(step)
    annotation_step = analyse_step.project.find_previous_step_of_type(analyse_step, 'annotations')

    for row in analyse_step.rows_as_dicts():
        seq_ident = row['AccesionNumber']
        starts = row['GeSeq part starts']
        offset = row['GeSeq trnF-GAA']
        orientation = row['GeSeq orientation']
        if not starts and subset != 'ge_seq':
            starts = row['NCBI part starts']
            offset = row['NCBI trnF-GAA']
            orientation = row['NCBI orientation']

        if offset in (None, ''):
            if subset == 'all':
                print(f"Warning: sequence {seq_ident} doesn't have trnF-GAA gene!")
                _copy_from_origin(step, annotation_step, seq_ident)  # ???
            continue

        if starts:
            starts = [int(f.strip()) for f in starts.split(',')]
            # Nothing to do!
            if abs(starts[0]) <= keep_offset and \
               abs(offset) <= keep_offset and \
               abs(starts[0] + offset) <= keep_offset \
               and not orientation:
                _copy_from_origin(step, annotation_step, seq_ident)
                continue

        # Check is sequence already in the CommonDB
        new_seq_ident = step.seq_ident_of_our_change(seq_ident, 'f')
        if sequences_db and (f := sequences_db.get_record(new_seq_ident, step.directory, info=True)):
            step.add_sequence_file(os.path.basename(f))
            continue

        seq_rec = annotation_step.get_sequence_record(seq_ident)
        new_seq = orient_by_trnF_GAA_by_data(seq_rec, offset, orientation, starts=starts, keep_offset=keep_offset)
        _store_genbank(step, new_seq_ident, new_seq, seq_rec, sequences_db)

    #
    step.save()
    return step


# def fix_by_trnH_GUG(step_data, analyse_step, keep_offset, sequences_db, annotations_db):
#     step = SequencesStep(analyse_step.project, step_data, remove_data=True)
#     analyse_step.propagate_step_name_prefix(step)
#     annotation_step = analyse_step.project.find_previous_step_of_type(analyse_step, 'annotations')

#     for row in analyse_step.rows_as_dicts():
#         seq_ident = row['AccesionNumber']
#         if (trnh_gug := row['trnH-GUG']) in (None, ''):
#             print(f"Warning: sequence {seq_ident} doesn't have trnH-GUG gene!")
#             # _copy_from_origin(step, annotation_step, seq_ident)  # ???
#             continue

#         fields = trnh_gug.split(':')
#         offset = int(fields[0])
#         zero_reverse = False
#         orientation = row['Orientation']
#         if starts := row['Part starts']:
#             starts = [int(f.strip()) for f in starts.split(',')]
#             # Nothing to do!
#             if abs(starts[0]) <= keep_offset and \
#                abs(offset) <= keep_offset and \
#                abs(starts[0] + offset) <= keep_offset \
#                and not orientation:
#                 _copy_from_origin(step, annotation_step, seq_ident)
#                 continue
#         else:
#             zero_reverse = (len(fields) > 1)

#         # Check is sequence already in the CommonDB
#         new_seq_ident = step.seq_ident_of_our_change(seq_ident, 'h')
#         if sequences_db and (f := sequences_db.get_record(new_seq_ident, step.directory, info=True)):
#             step.add_sequence_file(os.path.basename(f))
#             continue

#         seq_rec = annotation_step.get_sequence_record(seq_ident)
#         new_seq = orient_by_trnH_GUG_by_data(
#             seq_rec, offset, zero_reverse, orientation, starts=starts, keep_offset=keep_offset)

#         # ODKOMENTIRATI KAD PRORADI!!!
#         # _store_genbank(step, new_seq_ident, new_seq, seq_rec, sequences_db, annotations_db)
#         _store_genbank(step, new_seq_ident, new_seq, seq_rec, None, None)

#     #
#     step.save()
#     return step


def _store_genbank(step, new_seq_ident, seq_rec, copy_from_rec, sequences_db):
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


def _copy_sequence_annotations(from_rec, to_seq):
    if from_rec and to_seq:
        f_ann = from_rec.annotations
        t_ann = to_seq.annotations
        for k in ('molecule_type', 'topology', 'accessions', 'organism', 'taxonomy'):
            if v := f_ann.get(k):
                t_ann[k] = v


#
class _Stat:
    def __init__(self):
        self.taxid_2_stat = dict()  # Taxid -> dict (parent, nums)

    def add(self, add_to, ge_seq_parts, ge_seq_orient, ncbi_parts):
        # print(add_to, self.taxid_2_stat)
        for taxid, parent in zip(add_to, [None] + add_to[:-1]):
            if node := self.taxid_2_stat.get(taxid):
                assert node['parent'] == parent, (taxid, node['parent'], parent)
            else:
                self.taxid_2_stat[taxid] = node = dict(
                    taxid=taxid, parent=parent, num=0,
                    ge_seq_irs=0, ge_seq_offset=0, ge_seq_ssc=0, ge_seq_lsc=0, ge_seq_ir=0,
                    ncbi_irs=0)

            #
            node['num'] += 1
            if ge_seq_parts:
                node['ge_seq_irs'] += 1
                if abs(ge_seq_parts[0]) > DEFAULT_KEEP_OFFSET:
                    node['ge_seq_offset'] += 1
                if ge_seq_orient:
                    for p in ('ssc', 'lsc', 'ir'):
                        if p in ge_seq_orient:
                            node[f'ge_seq_{p}'] += 1
            if ncbi_parts:
                node['ncbi_irs'] += 1
            # print('  ', taxid, node)

    def get(self, parent):
        # ToDo: more efficient
        return [(s, self.get(s['taxid'])) for s in self.taxid_2_stat.values() if s['parent'] == parent]

    def print_all(self, taxid_translator, minimum_sequences):
        self._print_all(taxid_translator, self.get(None), 0, minimum_sequences)

    def _print_all(self, taxid_translator, stats, ident, minimum_sequences):
        names = taxid_translator([s['taxid'] for s, _ in stats])
        if minimum_sequences:
            rest = [s for s, _ in stats if s['num'] < minimum_sequences]
            if rest:
                stats = sorted(((s, sub_s) for s, sub_s in stats if s['num'] >= minimum_sequences),
                               key=lambda x: names[x[0]['taxid']])
                stats.append((self._sum(rest), []))
        else:
            stats = sorted(stats, key=lambda x: names[x[0]['taxid']])
        for s, sub_s in stats:
            n = '  ' * ident + (names[s['taxid']] if s['taxid'] else f"...({s['num_groups']})")
            parts = f"{s['ge_seq_ssc']}/{s['ge_seq_lsc']}/{s['ge_seq_ir']}"
            print(f"{n:18} {s['num']:3} : {s['ncbi_irs']:3}, {s['ge_seq_irs']:3} {s['ge_seq_offset']:3} {parts:>8}")
            if sub_s:
                self._print_all(taxid_translator, sub_s, ident + 1, minimum_sequences)

    @staticmethod
    def _sum(stats):
        assert stats
        ss = dict(taxid=None, parent=stats[0]['parent'], num=0, num_groups=len(stats),
                  ge_seq_irs=0, ge_seq_offset=0, ge_seq_ssc=0, ge_seq_lsc=0, ge_seq_ir=0,
                  ncbi_irs=0)
        for s in stats:
            for x in ('num', 'ge_seq_irs', 'ge_seq_offset', 'ge_seq_ssc', 'ge_seq_lsc', 'ge_seq_ir', 'ncbi_irs'):
                ss[x] += s[x]
        return ss


def statistics_by_taxa(project, analyse_step, taxa_ranks, minimum_sequences):
    from itertools import chain
    from ..utils.import_methods import import_ete3_NCBITaxa

    # Collect taxonomy data
    ncbi_taxonomy = import_ete3_NCBITaxa()()
    table_step = project.find_previous_step_of_type(analyse_step, 'table')
    nc_2_taxid = table_step.mapping_between_columns('ncbi_ident', 'tax_id')
    all_taxids = set(nc_2_taxid.values())
    taxid_2_lineage = dict((t, ncbi_taxonomy.get_lineage(t)) for t in all_taxids)
    all_taxids.update(chain.from_iterable(taxid_2_lineage.values()))
    taxid_2_rank = ncbi_taxonomy.get_rank(all_taxids)

    # Test
    # for lin in taxid_2_lineage.values():
    #     print([taxid_2_rank[t] for t in lin])

    # Extract only taxid of interest
    taxid_2_rank = dict((t, r) for t, r in taxid_2_rank.items() if r in taxa_ranks)
    print(taxid_2_rank)

    # Collect stat data
    _parts = lambda ps: list(map(int, ps.split(','))) if ps else None
    stat = _Stat()
    annotation_step = project.find_previous_step_of_type(analyse_step, 'annotations')
    for row in analyse_step.rows_as_dicts():
        seq_ident = row['AccesionNumber']
        taxid = nc_2_taxid[seq_ident]
        parents = [t for t in taxid_2_lineage[taxid] if t in taxid_2_rank]
        assert len(parents) == len(taxa_ranks), (len(parents), len(taxa_ranks), taxid, nc_2_taxid[seq_ident], parents)
        stat.add(parents, _parts(row['GeSeq part starts']), row['GeSeq orientation'], _parts(row['NCBI part starts']))

    #
    stat.print_all(ncbi_taxonomy.get_taxid_translator, minimum_sequences)
