from common_utils.pie_charts import pie_chart_bar_of_pie, pie_charts_bar_of_pie, Value
from .constants import DEFAULT_KEEP_OFFSET
from zci_bio.utils.stat_by_taxonomy import StatByTaxonomy


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
            bar_ratios = [node['wrong_offset'], node['wrong_lsc'], node['wrong_ssc'],
                          node['wrong_ir'], node['wrong_more']]
            bar_labels = ['Cyclic shift', 'Inverted LSC', 'Inverted SSC', 'Inverted IR', 'Combination']
            bar_labels = [f'{l} (n={n})' for l, n in zip(bar_labels, bar_ratios)]
            colors = ['#FF2222', '#DD0000', '#FF0000', '#BB0000', '#FF4444']
            # colors = ['darksalmon', 'xkcd:cherry red', 'salmon', 'xkcd:strawberry', 'lightsalmon']

            return dict(
                pie=[Value('Not in\nstandard\nform', node['wrong'], 'red'),
                     Value('In\nstandard\nform', node['standard'], 'green'),
                     Value('IRs not\nannotated', node['num'] - node['with_irs'], 'grey')],
                bar=[Value(l, v, c) for l, v, c in zip(bar_labels, bar_ratios, colors)],
                title=f"{taxa_name[0].upper()}{taxa_name[1:]} (n={node['num']})",
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
