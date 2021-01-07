from math import floor, log10
from itertools import product
from step_project.common.table.steps import TableStep
from common_utils.file_utils import get_settings, write_str_in_file
from common_utils.exceptions import ZCItoolsValueError
from common_utils.cache import cache_args
from common_utils.terminal_layout import StringColumns, fill_rows
from ..utils.phylogenetic_tree import PhylogeneticTree


def _most_sign(num, digits):
    return str(round(num, digits - int(floor(log10(abs(num)))) - 1))


class _TreeDiffs:
    def __init__(self, norm_result, seq_type):
        self.seq_type = s = seq_type
        get_tree = norm_result.trees.get

        # MrBayes - RAxML; Same alignments, same partition data, diff phylogenetics methods
        self.mr_bayes_2_raxml = dict(
            (on + wg, get_tree(f'{s}{on}_04_{wg}_MrBayes').distance_robinson_foulds(get_tree(f'{s}{on}_04_{wg}_RAxML')))
            for on, wg in product('on', 'WG'))

        # Whole - Genes; Same alignments, diff partition data, same phylogenetics methods
        on_MR = list(product('on', ('MrBayes', 'RAxML')))
        self.whole_2_genes_RF = dict(
            (on + m[0], get_tree(f'{s}{on}_04_W_{m}').distance_robinson_foulds(get_tree(f'{s}{on}_04_G_{m}')))
            for on, m in on_MR)
        self.whole_2_genes_BS = dict(
            (on + m[0], get_tree(f'{s}{on}_04_W_{m}').distance_branche_score(get_tree(f'{s}{on}_04_G_{m}')))
            for on, m in on_MR)
        self.whole_2_genes_KC = dict(
            (on + m[0], get_tree(f'{s}{on}_04_W_{m}').distance_kendall_colijn(get_tree(f'{s}{on}_04_G_{m}')))
            for on, m in on_MR)

        # original-normalized
        WG_MR = list(product('WG', ('MrBayes', 'RAxML')))
        self.orig_2_norm_RF = dict(
            (wg + m[0], get_tree(f'{s}o_04_{wg}_{m}').distance_robinson_foulds(get_tree(f'{s}n_04_{wg}_{m}')))
            for wg, m in WG_MR)
        self.orig_2_norm_BS = dict(
            (wg + m[0], get_tree(f'{s}o_04_{wg}_{m}').distance_branche_score(get_tree(f'{s}n_04_{wg}_{m}')))
            for wg, m in WG_MR)
        self.orig_2_norm_KC = dict(
            (wg + m[0], get_tree(f'{s}o_04_{wg}_{m}').distance_kendall_colijn(get_tree(f'{s}n_04_{wg}_{m}')))
            for wg, m in WG_MR)

    def get_rows(self):
        # Distance result format methods

        # Robinson-Foulds distance
        # Keys: rf, max_rf, ref_edges_in_source, source_edges_in_ref, effective_tree_size,
        #       norm_rf, treeko_dist, source_subtrees, common_edges, source_edges, ref_edges
        def _rf(d):
            return f"{int(d['rf'])}/{int(d['max_rf'])}"

        # Branch score distance
        def _bs(d):
            av_l = (d['stat_1']['average_length'] + d['stat_2']['average_length']) / 2
            return f"{_most_sign(d['bs'], 3)}/{_most_sign(av_l, 3)}"

        # Kendall-Colijn distance
        def _kc(d):
            return _most_sign(d, 3)

        rows = [
            [f'Seqs {self.seq_type}'],
            ['MrBayes-RAxML'],
            ['', 'Distance', 'Whole (no partition)', 'Genes (partition)'],
            ['original', 'RF', _rf(self.mr_bayes_2_raxml['oW']), _rf(self.mr_bayes_2_raxml['oG'])],
            ['normalized', 'RF', _rf(self.mr_bayes_2_raxml['nW']), _rf(self.mr_bayes_2_raxml['nG'])],
            #
            [],
            ['Whole-Genes'],
            ['', 'Distance', 'MrBayes', 'RAxML'],
            ['original', 'RF', _rf(self.whole_2_genes_RF['oM']), _rf(self.whole_2_genes_RF['oR'])],
            ['', 'BS', _bs(self.whole_2_genes_BS['oM']), _bs(self.whole_2_genes_BS['oR'])],
            ['', 'KC', _kc(self.whole_2_genes_KC['oM']), _kc(self.whole_2_genes_KC['oR'])],
            ['normalized', 'RF', _rf(self.whole_2_genes_RF['nM']), _rf(self.whole_2_genes_RF['nR'])],
            ['', 'BS', _bs(self.whole_2_genes_BS['nM']), _bs(self.whole_2_genes_BS['nR'])],
            ['', 'KC', _kc(self.whole_2_genes_KC['nM']), _kc(self.whole_2_genes_KC['nR'])],
            #
            [],
            ['original-normalized'],
            ['', 'Distance', 'MrBayes', 'RAxML'],
            ['Whole', 'RF', _rf(self.orig_2_norm_RF['WM']), _rf(self.orig_2_norm_RF['WR'])],
            ['', 'BS', _bs(self.orig_2_norm_BS['WM']), _bs(self.orig_2_norm_BS['WR'])],
            ['', 'KC', _kc(self.orig_2_norm_KC['WM']), _kc(self.orig_2_norm_KC['WR'])],
            ['Genes', 'RF', _rf(self.orig_2_norm_RF['GM']), _rf(self.orig_2_norm_RF['GR'])],
            ['', 'BS', _bs(self.orig_2_norm_BS['GM']), _bs(self.orig_2_norm_BS['GR'])],
            ['', 'KC', _kc(self.orig_2_norm_KC['GM']), _kc(self.orig_2_norm_KC['GR'])],
        ]
        fill_rows(rows)
        return rows


class NormalizationResult:
    def __init__(self, project):
        self.project = project
        self.outgroup = None
        self.analyses_step = None
        self.has_A = False
        self.tree_steps = dict()  # step_name -> step object
        self.trees = dict()       # step_name -> PhylogeneticTree object
        #
        self.S_tree_diffs = None  # Of type _TreeDiffs
        self.A_tree_diffs = None  # Of type _TreeDiffs

    def run(self, step_data):
        self._find_project_data()

        self.S_tree_diffs = _TreeDiffs(self, 'S')
        rows = self.S_tree_diffs.get_rows()
        text = str(StringColumns(rows))
        if self.has_A:
            self.A_tree_diffs = _TreeDiffs(self, 'A')
            rows_a = self.A_tree_diffs.get_rows()
            text += '\n\n' + str(StringColumns(rows_a))
            rows.extend(rows_a)

        # Create step and collect data
        step = TableStep(self.project, step_data, remove_data=True)
        self.analyses_step.propagate_step_name_prefix(step)
        step.set_table_data(rows, [(f'c_{i}', 'str') for i in range(1, 5)])
        step.save()  # completed=True
        #
        write_str_in_file(step.step_file('summary.txt'), text)
        print(text)

        return step

    def _find_project_data(self):
        # Outgroup
        if (settings := get_settings()) and (wf := settings.get('workflow_parameters')):
            self.outgroup = wf.get('outgroup')
        if not self.outgroup:
            print('Info: No outgroup specified?!', settings)

        # Analyses chloroplast step
        self.analyses_step = self.project.read_step_if_in(
                '04_AnalyseChloroplast', check_data_type='table', no_check=True)
        if not self.analyses_step:
            raise ZCItoolsValueError('No analyse chloroplast step (04_AnalyseChloroplast)!')
        self.has_A = wf.get('calc_all', 0) and not all(self.analyses_step.get_column_values('Part starts'))

        # Find all phylogenetic steps
        for sa in ('SA' if self.has_A else 'S'):
            for on in 'on':
                for wg in 'WG':
                    for phylo, dt in (('MrBayes', 'mr_bayes'), ('RAxML', 'raxml')):
                        step_name = f'{sa}{on}_04_{wg}_{phylo}'
                        step = self.project.read_step_if_in(step_name, check_data_type=dt, no_check=True)
                        self.tree_steps[step_name] = step
                        self.trees[step_name] = PhylogeneticTree(
                            step.get_consensus_file(), self.outgroup,
                            rename_nodes=lambda name: f'NC_{name[2:]}' if name.startswith('p_') else name)

        if (no_steps := [k for k, v in self.tree_steps.items() if not v]):
            raise ZCItoolsValueError(f"No tree step(s): {', '.join(sorted(no_steps))}!")
