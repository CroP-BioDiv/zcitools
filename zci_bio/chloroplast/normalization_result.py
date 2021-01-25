from math import floor, log10
from itertools import product, chain
from step_project.common.table.steps import TableStep
from common_utils.file_utils import get_settings, write_str_in_file
from common_utils.exceptions import ZCItoolsValueError
from common_utils.cache import cache_args
from common_utils.terminal_layout import StringColumns, fill_rows
from ..utils.phylogenetic_tree import PhylogeneticTree


def _most_sign(num, digits):
    if num >= 10 ** (digits - 1):
        return int(num)
    return str(round(num, digits - int(floor(log10(abs(num)))) - 1))


# Distance result format methods
# Robinson-Foulds distance
# Keys: rf, max_rf, ref_edges_in_source, source_edges_in_ref, effective_tree_size,
#       norm_rf, treeko_dist, source_subtrees, common_edges, source_edges, ref_edges
def _rf(d):
    return f"{int(d['rf'])}/{int(d['max_rf'])}"


# Branch score distance
def _bs(d):
    av_l = (d['stat_1']['average_length'] + d['stat_2']['average_length']) / 2
    num_edges = d['bs'] / av_l
    return f"{_most_sign(d['bs'], 2)}/{_most_sign(av_l, 2)} ({_most_sign(num_edges, 3)})"


# Kendall-Colijn distance
def _kc(d):
    return _most_sign(d, 3)


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

    def get_rows(self, *tree_diffs):
        assert all(isinstance(t, _TreeDiffs) for t in tree_diffs), tree_diffs
        tds = (self,) + tree_diffs
        nt = len(tds)
        cc = lambda x: list(chain(*x))
        return fill_rows([
            ['', ''] + cc([f'Seqs {t.seq_type}', ''] for t in tds),
            ['MrBayes-RAxML'],
            ['', 'Distance'] + ['Whole (no partition)', 'Genes (partition)'] * nt,
            ['original', 'RF'] + cc([_rf(t.mr_bayes_2_raxml['oW']), _rf(t.mr_bayes_2_raxml['oG'])] for t in tds),
            ['normalized', 'RF'] + cc([_rf(t.mr_bayes_2_raxml['nW']), _rf(t.mr_bayes_2_raxml['nG'])] for t in tds),
            #
            [],
            ['Whole-Genes'],
            ['', 'Distance'] + ['MrBayes', 'RAxML'] * nt,
            ['original', 'RF'] + cc([_rf(t.whole_2_genes_RF['oM']), _rf(t.whole_2_genes_RF['oR'])] for t in tds),
            ['', 'BS'] + cc([_bs(t.whole_2_genes_BS['oM']), _bs(t.whole_2_genes_BS['oR'])] for t in tds),
            ['', 'KC'] + cc([_kc(t.whole_2_genes_KC['oM']), _kc(t.whole_2_genes_KC['oR'])] for t in tds),
            ['normalized', 'RF'] + cc([_rf(t.whole_2_genes_RF['nM']), _rf(t.whole_2_genes_RF['nR'])] for t in tds),
            ['', 'BS'] + cc([_bs(t.whole_2_genes_BS['nM']), _bs(t.whole_2_genes_BS['nR'])] for t in tds),
            ['', 'KC'] + cc([_kc(t.whole_2_genes_KC['nM']), _kc(t.whole_2_genes_KC['nR'])] for t in tds),
            #
            [],
            ['original-normalized'],
            ['', 'Distance'] + ['MrBayes', 'RAxML'] * nt,
            ['Whole', 'RF'] + cc([_rf(t.orig_2_norm_RF['WM']), _rf(t.orig_2_norm_RF['WR'])] for t in tds),
            ['', 'BS'] + cc([_bs(t.orig_2_norm_BS['WM']), _bs(t.orig_2_norm_BS['WR'])] for t in tds),
            ['', 'KC'] + cc([_kc(t.orig_2_norm_KC['WM']), _kc(t.orig_2_norm_KC['WR'])] for t in tds),
            ['Genes', 'RF'] + cc([_rf(t.orig_2_norm_RF['GM']), _rf(t.orig_2_norm_RF['GR'])] for t in tds),
            ['', 'BS'] + cc([_bs(t.orig_2_norm_BS['GM']), _bs(t.orig_2_norm_BS['GR'])] for t in tds),
            ['', 'KC'] + cc([_kc(t.orig_2_norm_KC['GM']), _kc(t.orig_2_norm_KC['GR'])] for t in tds),
        ])


class NormalizationResult:
    def __init__(self, project):
        self.project = project
        self.outgroup = None
        self.analyses_step = None
        self.tree_steps = dict()  # step_name -> step object
        self.trees = dict()       # step_name -> PhylogeneticTree object

    def run(self, step_data):
        self._find_project_data()

        G_tree_diffs = _TreeDiffs(self, 'G')
        # rows_g = G_tree_diffs.get_rows()
        # text = str(StringColumns(rows_g))

        N_tree_diffs = _TreeDiffs(self, 'N')
        # rows_n = N_tree_diffs.get_rows()
        # text += '\n\n' + str(StringColumns(rows_n))

        rows = G_tree_diffs.get_rows(N_tree_diffs)
        text = str(StringColumns(rows))

        # Create step and collect data
        step = TableStep(self.project, step_data, remove_data=True)
        self.analyses_step.propagate_step_name_prefix(step)
        step.set_table_data(rows, [(f'c_{i}', 'str') for i in range(1, 5)])
        # rows_g + rows_n
        step.save()  # completed=True
        #
        write_str_in_file(step.step_file('summary.txt'), text)
        print(text)

        return step

    def _find_project_data(self):
        from ..workflows.chloroplast_normalization import workflow_branches

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

        analyses_branches = workflow_branches(wf, self.analyses_step)

        # Find all phylogenetic steps
        for sa in analyses_branches:
            for on in 'on':
                for wg in 'WG':
                    for phylo, dt in (('MrBayes', 'mr_bayes'), ('RAxML', 'raxml')):
                        step_name = f'{sa}{on}_04_{wg}_{phylo}'
                        step = self.project.read_step_if_in(step_name, check_data_type=dt, no_check=True)
                        self.tree_steps[step_name] = step
                        self.trees[step_name] = PhylogeneticTree(
                            step.get_consensus_file(), self.outgroup,
                            rename_nodes=_rename_seq_name)

        if (no_steps := [k for k, v in self.tree_steps.items() if not v]):
            raise ZCItoolsValueError(f"No tree step(s): {', '.join(sorted(no_steps))}!")


def _rename_seq_name(name):
    if name.startswith('p_') or name.startswith('n_'):
        return f'NC_{name[2:]}'
    return name
