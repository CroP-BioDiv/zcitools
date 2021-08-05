from itertools import chain
from common_utils.cache import cache
from zci_bio.utils.ncbi_taxonomy import get_ncbi_taxonomy


class StatByTaxonomy:
    _STAT_ATTRS = []
    _EXCEL_COLUMNS = []

    def __init__(self, taxids, ranks=None, names=None):
        self.taxid_2_stat = dict()  # Taxid -> dict (parent, nums)

        # Find all taxids, species and parent clades
        ncbi_taxonomy = get_ncbi_taxonomy()
        all_taxids = set(taxids)
        gl = ncbi_taxonomy.get_lineage
        self.taxid_2_lineage = dict((t, gl(t)) for t in all_taxids)
        all_taxids.update(chain.from_iterable(self.taxid_2_lineage.values()))
        self.taxid_2_rank = ncbi_taxonomy.get_rank(all_taxids)
        self.taxid_2_name = ncbi_taxonomy.get_taxid_translator(all_taxids)

        # Extract only taxid of interest
        self.group_taxids = set()
        if ranks:
            self.group_taxids.update(t for t, r in self.taxid_2_rank.items() if r in ranks)
        if names:
            for name in names:
                if n_taxid := ncbi_taxonomy.get_exact_name_translation(name):
                    self.group_taxids.add(n_taxid)

    def stat_attrs(self):
        return self._STAT_ATTRS

    def excel_columns(self):
        return self._EXCEL_COLUMNS

    @cache
    def name_2_taxid(self):
        return dict((n, t) for t, n in self.taxid_2_name.items())

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
                self._init_new_node(node)
                self.taxid_2_stat[taxid] = node
            #
            node['num'] += 1
            self._add(node, *args, **kwargs)

    def _init_new_node(self, node):
        node.update((a, 0) for a in self.stat_attrs())

    def _sum(self, stats):
        # Note: contains num_groups attr!
        assert stats
        ss = dict(taxid=None, parent=stats[0]['parent'], num=0, num_groups=len(stats))
        self._init_new_node(ss)
        for a in ['num'] + self.stat_attrs():
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
        # Collaps groups with sequences less than given minimum
        # Note: Single leaf is not collapsed. Single inner node is collapsed.
        if minimum_sequences and \
           (rest := [s for s, _ in stats if s['num'] < minimum_sequences]) \
           and (len(rest) > 1 or depth < max_depth - 1):
            stats = sorted(((s, sub_s) for s, sub_s in stats if s['num'] >= minimum_sequences),
                           key=lambda x: self.taxid_2_name[x[0]['taxid']])
            stats.append((self._sum(rest), []))
        else:
            stats = sorted(stats, key=lambda x: self.taxid_2_name[x[0]['taxid']])

        rows = []
        for s, sub_s in stats:
            rows.append(self._to_row(s, max_depth, depth))
            if sub_s:
                rows += self._to_table(sub_s, minimum_sequences, max_depth, depth + 1)

        # Add sum row
        if depth == 0 and len(stats) > 1:
            rows.append(self._to_sum_row(self._sum([s for s, _ in stats]), max_depth))

        #
        return rows

    def _to_row_label(self, s, max_depth, depth, label):
        row = [''] * max_depth + [s['num']] + [s[a] for a in self.stat_attrs()]
        row[depth] = label
        return row

    def _to_row(self, s, max_depth, depth):
        return self._to_row_label(
            s, max_depth, depth, (self.taxid_2_name[s['taxid']] if s['taxid'] else f"... ({s['num_groups']})"))

    def _to_sum_row(self, s, max_depth):
        return self._to_row_label(s, max_depth, 0, 'all')

    def export_excel(self, filename, minimum_sequences):
        from common_utils.value_data_types import rows_2_excel
        md, rows = self.to_table(minimum_sequences)
        columns = [''] * md + ['Num sequences'] + self.excel_columns()
        rows_2_excel(filename, columns, rows, index=False)

    def export_sheet(self, sheet_name, minimum_sequences):
        md, rows = self.to_table(minimum_sequences)
        columns = [''] * md + ['Num sequences'] + self.excel_columns()
        return (sheet_name, columns, rows)

    #
    def get_node_from_name(self, taxa_name):
        if not (taxid := self.name_2_taxid().get(taxa_name)):
            print(f'No taxon with name {taxa_name}!')
            return
        if taxid not in self.group_taxids:
            print(f'Taxon {taxa_name} is not used for grouping!')
            return
        if not (node := self.taxid_2_stat.get(taxid)):
            print(f'Taxon {taxa_name} has no stat data!')
            return
        return node
