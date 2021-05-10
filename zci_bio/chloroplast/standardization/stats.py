from itertools import chain
from collections import namedtuple
import math
from common_utils.cache import cache
from common_utils.import_method import import_matplotlib_pyplot
from zci_bio.utils.ncbi_taxonomy import get_ncbi_taxonomy
from .constants import DEFAULT_KEEP_OFFSET


class _Value(namedtuple('_Value', 'label, value, color')):
    def __new__(cls, label, value, color=None):
        return super(_Value, cls).__new__(cls, label, value, color)


def _values_2_lists(values, normalize=True, revert=False):
    labels = [v.label for v in values if v.value]
    ratios = [v.value for v in values if v.value]
    colors = [v.color for v in values if v.value and v.color]
    if revert:
        labels = labels[::-1]
        ratios = ratios[::-1]
        colors = colors[::-1]
    #
    if len(colors) != len(ratios):
        if colors:
            print('Warning: not all colors are specified!')
        colors = None

    if normalize:
        s_ratios = sum(ratios)
        if abs(s_ratios - 1) > 0.01:
            ratios = [x / s_ratios for x in ratios]
    return ratios, labels, colors


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

    def _to_row(self, s, max_depth, depth):
        row = [''] * max_depth + [s['num']] + [s[a] for a in self._STAT_ATTRS]
        row[depth] = (self.taxid_2_name[s['taxid']] if s['taxid'] else f"... ({s['num_groups']})")
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


#
class _Stat(StatByTaxonomy):
    _STAT_ATTRS = ['with_irs', 'standard', 'wrong',
                   'wrong_more',
                   'wrong_offset', 'wrong_ssc', 'wrong_lsc', 'wrong_ir',
                   'only_offset', 'only_ssc', 'only_lsc', 'only_ir']
    _EXCEL_COLUMNS = ['Annotated IRs', 'Standardized', 'Wrong',
                      'More',
                      'Offset', 'SSC', 'LSC', 'IRs',
                      'Only offset', 'Only SSC', 'Only LSC', 'Only IRs']

    def _add(self, node, parts, orient):
        if parts:
            node['with_irs'] += 1
            offset = abs(parts[0]) > DEFAULT_KEEP_OFFSET
            ps = [p for p in ('ssc', 'lsc', 'ir') if p in orient] if orient else []
            wrong_num = int(offset) + len(ps)

            if not wrong_num:
                node['standard'] += 1
            elif wrong_num == 1:
                node['wrong'] += 1
                if offset:
                    node['wrong_offset'] += 1
                for p in ps:
                    node[f'wrong_{p}'] += 1
            else:
                node['wrong'] += 1
                node['wrong_more'] += 1
                if offset:
                    node['only_offset'] += 1
                for p in ps:
                    node[f'only_{p}'] += 1

    #
    def print_all(self, minimum_sequences):
        md, rows = self.to_table(minimum_sequences)
        for row in rows:
            for ident, x in enumerate(row):
                if x:
                    break
            n = '  ' * ident + row[ident]
            parts = f"{row[md + 5]}/{row[md + 6]}/{row[md + 7]}"
            print(f"{n:18} {row[md]:3} : {row[md + 1]:3}, {row[md + 2]:3} {row[md + 3]:3} {parts:>8} {row[md + 4]:3}")

    def _pie_chart_data(self, taxa_name):
        if node := self.get_node_from_name(taxa_name):
            bar_ratios = [node['wrong_more'], node['wrong_offset'],
                          node['wrong_lsc'], node['wrong_ssc'], node['wrong_ir']]
            bar_labels = ['More', 'Offset', 'LSC', 'SSC', 'IR']
            bar_labels = [f'{l} ({n})' for l, n in zip(bar_labels, bar_ratios)]

            return dict(
                pie=[_Value('Wrong', node['wrong']),
                     _Value('Standard', node['standard']),
                     _Value('No IRs', node['num'] - node['with_irs'])],
                bar=[_Value(l, v) for l, v in zip(bar_labels, bar_ratios)],
                title=f"Clade: {taxa_name} ({node['num']})",
            )

    def make_pie_chart(self, taxa_name):
        if d := self._pie_chart_data(taxa_name):
            pie_chart_bar_of_pie(d, explode=0.06, output_filename=f"pie_chart_{taxa_name}")

    def make_pie_charts(self, taxa_names):
        ds = [self._pie_chart_data(t) for t in taxa_names]
        if all(ds):
            pie_charts_bar_of_pie(ds, explode=0.06, output_filename='pie_charts')


def statistics_by_taxa(project, analyse_step, taxa_ranks, taxa_names, minimum_sequences,
                       output_excel, pie_chart_name, merge_pie_charts, print_output):
    # Collect taxonomy data
    table_step = project.find_previous_step_of_type(analyse_step, 'table')
    nc_2_taxid = table_step.mapping_between_columns('ncbi_ident', 'tax_id')
    stat = _Stat(nc_2_taxid.values(), ranks=taxa_ranks, names=taxa_names)

    # Collect stat data
    _parts = lambda ps: list(map(int, ps.split(','))) if ps else None
    # annotation_step = project.find_previous_step_of_type(analyse_step, 'annotations')
    for row in analyse_step.rows_as_dicts():
        seq_ident = row['AccesionNumber']
        stat.add(nc_2_taxid[seq_ident], _parts(row['GeSeq part starts']), row['GeSeq orientation'])
        #, _parts(row['NCBI part starts']))

    #
    if output_excel:
        stat.export_excel(output_excel, minimum_sequences)
    if print_output:
        stat.print_all(minimum_sequences)
    if pie_chart_name:
        if merge_pie_charts:
            stat.make_pie_charts(pie_chart_name)
        else:
            for pc in pie_chart_name:
                stat.make_pie_chart(pc)


def _normalize_ratios(ratios, labels):
    labels = [l for l, x in zip(labels, ratios) if x]
    ratios = [x for x in ratios if x]
    s_ratios = sum(ratios)
    if abs(s_ratios - 1) > 0.01:
        ratios = [x / s_ratios for x in ratios]
    return ratios, labels


# Pie charts
def _save_chart(pyplot, output_filename):
    if output_filename:
        pyplot.savefig(f'{output_filename}.svg')
        pyplot.savefig(f'{output_filename}.png', dpi=150)

    pyplot.show()


def pie_chart_bar_of_pie(pie_data, explode=None, output_filename=None):
    # https://matplotlib.org/stable/gallery/pie_and_polar_charts/bar_of_pie.html
    # pie_data is dict with attrs:
    #  - pie : list _Value objects
    #  - bar : list _Value objects
    #  - title : string
    #  - output_filename : string
    pyplot = import_matplotlib_pyplot()

    # make figure and assign axis objects
    fig, (ax1, ax2) = pyplot.subplots(1, 2, figsize=(9, 5))
    fig.subplots_adjust(wspace=0)

    _pie_chart_bar_of_pie(ax1, ax2, pie_data, explode=explode)
    _save_chart(pyplot, output_filename)


def pie_charts_bar_of_pie(pies_data, explode=None, output_filename=None):
    pyplot = import_matplotlib_pyplot()

    # make figure and assign axis objects
    fig, axis = pyplot.subplots(1, 3 * len(pies_data) - 1, figsize=(9, 3))
    fig.subplots_adjust(wspace=0)
    for i in range(1, len(pies_data)):
        axis[i * 3 - 1].axis('off')

    for i, pie_data in enumerate(pies_data):
        _pie_chart_bar_of_pie(axis[i * 3], axis[i * 3 + 1], pie_data, explode=explode)

    _save_chart(pyplot, output_filename)


def _pie_chart_bar_of_pie(ax1, ax2, pie_data, explode=None):
    from matplotlib.patches import ConnectionPatch

    # pie chart parameters
    pie_ratios, pie_labels, pie_colors = _values_2_lists(pie_data['pie'])
    num_parts = len(pie_ratios)

    # Explode
    if explode:
        _explode = [0.1 if explode is True else explode] + [0] * (num_parts - 1)
    else:
        _explode = None

    # rotate so that first wedge is split by the x-axis
    angle = -180 * pie_ratios[0]
    ax1.pie(pie_ratios, autopct='%1.1f%%', startangle=angle, labels=pie_labels, explode=_explode, colors=pie_colors)
    if title := pie_data.get('title'):
        ax1.set_title(title)

    # bar chart parameters
    xpos = 0
    bottom = 0
    bar_ratios, bar_labels, bar_colors = _values_2_lists(pie_data['bar'], revert=True)
    num_bars = len(bar_ratios)
    bar_width = pie_data.get('bar_width', .2)
    if not (bar_colors := pie_data.get('bar_colors')):
        d = (0.9 - 0.3) / (num_bars - 1)
        bar_colors = [[0.1, 0.3, 0.3 + d * x] for x in range(num_bars)]

    x_label = bar_width * 0.6
    for j, height in enumerate(bar_ratios):
        ax2.bar(xpos, height, bar_width, bottom=bottom, color=bar_colors[j])
        ypos = bottom + ax2.patches[j].get_height() / 2
        bottom += height
        ax2.text(x_label, ypos, bar_labels[j], ha='left', va='center')

    # ax2.legend(('50-65', 'Over 65', '35-49', 'Under 35'))
    ax2.axis('off')
    ax2.set_xlim(- 2.5 * bar_width, 2.5 * bar_width)

    # use ConnectionPatch to draw lines between the two plots
    # get the wedge data
    theta1, theta2 = ax1.patches[0].theta1, ax1.patches[0].theta2
    center, r = ax1.patches[0].center, ax1.patches[0].r
    bar_height = sum([item.get_height() for item in ax2.patches])

    # draw top connecting line
    x = r * math.cos(math.pi / 180 * theta2) + center[0]
    y = r * math.sin(math.pi / 180 * theta2) + center[1]
    con = ConnectionPatch(xyA=(-bar_width / 2, bar_height), coordsA=ax2.transData,
                          xyB=(x, y), coordsB=ax1.transData)
    con.set_color('black')
    con.set_linewidth(1)
    con.set_linestyle((0, (1, 4)))
    ax2.add_artist(con)

    # draw bottom connecting line
    x = r * math.cos(math.pi / 180 * theta1) + center[0]
    y = r * math.sin(math.pi / 180 * theta1) + center[1]
    con = ConnectionPatch(xyA=(-bar_width / 2, 0), coordsA=ax2.transData,
                          xyB=(x, y), coordsB=ax1.transData)
    con.set_color('black')
    ax2.add_artist(con)
    con.set_linewidth(1)
    con.set_linestyle((0, (1, 4)))
