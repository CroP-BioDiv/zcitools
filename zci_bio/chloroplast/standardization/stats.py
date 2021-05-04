from itertools import chain
from zci_bio.utils.ncbi_taxonomy import get_ncbi_taxonomy
from .constants import DEFAULT_KEEP_OFFSET


class StatByTaxonomy:
    _STAT_ATTRS = []
    _EXCEL_COLUMNS = []

    def __init__(self, taxids, ranks=None, names=None):
        self.ncbi_taxonomy = get_ncbi_taxonomy()
        self.taxid_2_stat = dict()  # Taxid -> dict (parent, nums)

        # Find all taxids, species and parent clades
        all_taxids = set(taxids)
        gl = self.ncbi_taxonomy.get_lineage
        self.taxid_2_lineage = dict((t, gl(t)) for t in all_taxids)
        all_taxids.update(chain.from_iterable(self.taxid_2_lineage.values()))
        self.taxid_2_rank = self.ncbi_taxonomy.get_rank(all_taxids)

        # Extract only taxid of interest
        self.group_taxids = set()
        if ranks:
            self.group_taxids.update(t for t, r in self.taxid_2_rank.items() if r in ranks)
        if names:
            for name in names:
                if n_taxid := self.ncbi_taxonomy.get_exact_name_translation(name):
                    self.group_taxids.add(n_taxid)

    def add(self, taxid, *args, **kwargs):
        parents = [t for t in self.taxid_2_lineage[taxid] if t in self.group_taxids]
        if not parents:
            print(f'Warning: stat for taxid {taxid} is not set, because there is no parent to add to!')
            return

        # Add stat to parent nodes
        for taxid, parent in zip(parents, [None] + parents[:-1]):
            if node := self.taxid_2_stat.get(taxid):
                assert node['parent'] == parent, (taxid, node['parent'], parent)
            else:
                node = dict(taxid=taxid, parent=parent, num=0)
                node.update((a, 0) for a in self._STAT_ATTRS)
                self.taxid_2_stat[taxid] = node
            #
            node['num'] += 1
            self._add(node, *args, **kwargs)

    @classmethod
    def _sum(cls, stats):
        # Note: contains num_groups attr!
        assert stats
        ss = dict(taxid=None, parent=stats[0]['parent'], num=0, num_groups=len(stats))
        ss.update((a, 0) for a in cls._STAT_ATTRS)
        for a in ['num'] + cls._STAT_ATTRS:
            for s in stats:
                ss[a] += s[a]
        return ss

    #
    def get(self, parent):
        # ToDo: more efficient
        return [(s, self.get(s['taxid'])) for s in self.taxid_2_stat.values() if s['parent'] == parent]

    def _depth(self, stats):
        return (max(self._depth(sub_s) for _, sub_s in stats) + 1) if stats else 0

    def to_table(self, minimum_sequences):
        stats = self.get(None)
        max_depth = self._depth(stats)
        return max_depth, self._to_table(stats, minimum_sequences, max_depth, 0)

    def _to_table(self, stats, minimum_sequences, max_depth, depth):
        names = self.ncbi_taxonomy.get_taxid_translator([s['taxid'] for s, _ in stats])
        # Collaps groups with sequences less than given minimum
        # Note: Single leaf is not collapsed. Single inner node is collapsed.
        if minimum_sequences and \
           (rest := [s for s, _ in stats if s['num'] < minimum_sequences]) \
           and (len(rest) > 1 or depth < max_depth - 1):
            stats = sorted(((s, sub_s) for s, sub_s in stats if s['num'] >= minimum_sequences),
                           key=lambda x: names[x[0]['taxid']])
            stats.append((self._sum(rest), []))
        else:
            stats = sorted(stats, key=lambda x: names[x[0]['taxid']])

        rows = []
        for s, sub_s in stats:
            rows.append(self._to_row(s, names, max_depth, depth))
            if sub_s:
                rows += self._to_table(sub_s, minimum_sequences, max_depth, depth + 1)

        # Add sum row
        if depth == 0 and len(stats) > 1:
            rows.append(self._to_sum_row(self._sum([s for s, _ in stats]), max_depth))

        #
        return rows

    def _to_row(self, s, names, max_depth, depth):
        row = [''] * max_depth + [s['num']] + [s[a] for a in self._STAT_ATTRS]
        row[depth] = (names[s['taxid']] if s['taxid'] else f"... ({s['num_groups']})")
        return row

    def _to_sum_row(self, s, max_depth):
        row = [''] * max_depth + [s['num']] + [s[a] for a in self._STAT_ATTRS]
        row[0] = 'all'
        return row

    def export_excel(self, filename, minimum_sequences):
        from common_utils.value_data_types import rows_2_excel
        md, rows = self.to_table(minimum_sequences)
        columns = [''] * md + ['Num sequences'] + self._EXCEL_COLUMNS
        rows_2_excel(filename, columns, rows, index=False)


#
class _Stat(StatByTaxonomy):
    _STAT_ATTRS = ['ncbi_irs', 'ge_seq_irs', 'ge_seq_offset', 'ge_seq_ssc', 'ge_seq_lsc', 'ge_seq_ir']
    _EXCEL_COLUMNS = ['NCBI IRs', 'GeSeq IRs', 'Offset', 'SSC', 'LSC', 'IRs']

    def _add(self, node, ge_seq_parts, ge_seq_orient, ncbi_parts):
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

    #
    def print_all(self, minimum_sequences):
        md, rows = self.to_table(minimum_sequences)
        for row in rows:
            for ident, x in enumerate(row):
                if x:
                    break
            n = '  ' * ident + row[ident]
            parts = f"{row[md + 4]}/{row[md + 5]}/{row[md + 6]}"
            print(f"{n:18} {row[md]:3} : {row[md + 1]:3}, {row[md + 2]:3} {row[md + 3]:3} {parts:>8}")


def statistics_by_taxa(project, analyse_step, taxa_ranks, taxa_names, minimum_sequences, output_excel):
    # Collect taxonomy data
    table_step = project.find_previous_step_of_type(analyse_step, 'table')
    nc_2_taxid = table_step.mapping_between_columns('ncbi_ident', 'tax_id')
    stat = _Stat(nc_2_taxid.values(), ranks=taxa_ranks, names=taxa_names)

    # Collect stat data
    _parts = lambda ps: list(map(int, ps.split(','))) if ps else None
    # annotation_step = project.find_previous_step_of_type(analyse_step, 'annotations')
    for row in analyse_step.rows_as_dicts():
        seq_ident = row['AccesionNumber']
        stat.add(nc_2_taxid[seq_ident], _parts(row['GeSeq part starts']), row['GeSeq orientation'], _parts(row['NCBI part starts']))

    #
    if output_excel:
        stat.export_excel(output_excel, minimum_sequences)
    else:
        stat.print_all(minimum_sequences)
