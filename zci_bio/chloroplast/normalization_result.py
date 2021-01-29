from math import floor, ceil, log10
from itertools import product, chain
from step_project.common.table.steps import TableStep
from common_utils.file_utils import get_settings, write_str_in_file
from common_utils.exceptions import ZCItoolsValueError
from common_utils.cache import cache_args
from common_utils.terminal_layout import StringColumns, fill_rows
from common_utils.import_method import import_matplotlib_pylot
from ..utils.phylogenetic_tree import PhylogeneticTree

# In inches. A4 is 8-1/4 * 11-3/4.
# Note: figsize is value used to calibrate all other values
figsize_x = 7
figsize_y = 3
bottom_space = 0.2
figsize_y_no_x = figsize_y * (1 - bottom_space)
right_axis_space = 0.05
x_label_y = [-0.03, -0.08, -0.14]


def _rename_seq_name(name):
    if name.startswith('p_') or name.startswith('n_'):
        return f'NC_{name[2:]}'
    return name


def _most_sign(num, digits):
    if num == 0:
        return 0
    if num >= 10 ** (digits - 1):
        return int(num)
    return str(round(num, digits - int(floor(log10(abs(num)))) - 1))


# Distance result format methods
# Robinson-Foulds distance
# Keys: rf, max_rf, ref_edges_in_source, source_edges_in_ref, effective_tree_size,
#       norm_rf, treeko_dist, source_subtrees, common_edges, source_edges, ref_edges
def _rf(d):
    return f"{_most_sign(d[0] / d[1], 3)} ({int(d[0])}/{int(d[1])})"


def _rf_v(d):
    return d[0] / d[1]


# Branch score distance
def _bs(d):
    av_l = (d['stat_1']['average_length'] + d['stat_2']['average_length']) / 2
    num_edges = d['bs'] / av_l
    return f"{_most_sign(num_edges, 3)} ({_most_sign(d['bs'], 2)}/{_most_sign(av_l, 2)})"


def _bs_v(d):
    av_l = (d['stat_1']['average_length'] + d['stat_2']['average_length']) / 2
    return d['bs'] / av_l


# Kendall-Colijn distance
def _kc(d):
    return _most_sign(d, 3)


class _TreeDiffs:
    def __init__(self, norm_result, seq_type=None, data=None):
        if not seq_type:
            assert data
            self.from_dict(data)
            return

        self.seq_type = s = seq_type
        get_tree = norm_result.trees.get

        # MrBayes - RAxML; Same alignments, same partition data, diff phylogenetics methods
        self.mr_bayes_2_raxml_RF = dict(
            (on + wg, get_tree(f'{s}{on}_04_{wg}_MrBayes').distance_robinson_foulds(get_tree(f'{s}{on}_04_{wg}_RAxML')))
            for on, wg in product('on', 'WG'))
        self.mr_bayes_2_raxml_KCT = dict(
            (on + wg, get_tree(f'{s}{on}_04_{wg}_MrBayes').distance_kendall_colijn_topo(get_tree(f'{s}{on}_04_{wg}_RAxML')))
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
        self.whole_2_genes_KCT = dict(
            (on + m[0], get_tree(f'{s}{on}_04_W_{m}').distance_kendall_colijn_topo(get_tree(f'{s}{on}_04_G_{m}')))
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
        self.orig_2_norm_KCT = dict(
            (wg + m[0], get_tree(f'{s}o_04_{wg}_{m}').distance_kendall_colijn_topo(get_tree(f'{s}n_04_{wg}_{m}')))
            for wg, m in WG_MR)

    def to_dict(self):
        return self.__dict__
        # return dict((a, getattr(self, a)) for a in (
        #     'mr_bayes_2_raxml_RF', 'mr_bayes_2_raxml_KCT',
        #     'whole_2_genes_RF', 'whole_2_genes_BS', 'whole_2_genes_KC', 'whole_2_genes_KCT',
        #     'orig_2_norm_RF', 'orig_2_norm_BS', 'orig_2_norm_KC', 'orig_2_norm_KCT'))

    def from_dict(self, data):
        self.__dict__.update(data)

    def get_rows(self, *tree_diffs):
        assert all(isinstance(t, _TreeDiffs) for t in tree_diffs), tree_diffs
        tds = (self,) + tree_diffs
        nt = len(tds)
        cc = lambda x: list(chain(*x))
        lab_2_text = dict(G='GeSeq', S='Sum', A='All', N='NCBI')

        return fill_rows([
            ['', ''] + cc([lab_2_text[t.seq_type], ''] for t in tds),
            ['MrBayes-RAxML'],
            ['', 'Distance'] + ['Non-partitioned', 'Partitioned'] * nt,
            ['original', 'RF'] + cc([_rf(t.mr_bayes_2_raxml_RF['oW']), _rf(t.mr_bayes_2_raxml_RF['oG'])] for t in tds),
            ['', 'KCT'] + cc([_kc(t.mr_bayes_2_raxml_KCT['oW']), _kc(t.mr_bayes_2_raxml_KCT['oG'])] for t in tds),
            ['normalized', 'RF'] + cc([_rf(t.mr_bayes_2_raxml_RF['nW']), _rf(t.mr_bayes_2_raxml_RF['nG'])] for t in tds),
            ['', 'KCT'] + cc([_kc(t.mr_bayes_2_raxml_KCT['nW']), _kc(t.mr_bayes_2_raxml_KCT['nG'])] for t in tds),
            #
            [],
            ['Non-partitioned-Partitioned'],
            ['', 'Distance'] + ['MrBayes', 'RAxML'] * nt,
            ['original', 'RF'] + cc([_rf(t.whole_2_genes_RF['oM']), _rf(t.whole_2_genes_RF['oR'])] for t in tds),
            ['', 'BS'] + cc([_bs(t.whole_2_genes_BS['oM']), _bs(t.whole_2_genes_BS['oR'])] for t in tds),
            ['', 'KC'] + cc([_kc(t.whole_2_genes_KC['oM']), _kc(t.whole_2_genes_KC['oR'])] for t in tds),
            ['', 'KCT'] + cc([_kc(t.whole_2_genes_KCT['oM']), _kc(t.whole_2_genes_KCT['oR'])] for t in tds),
            ['normalized', 'RF'] + cc([_rf(t.whole_2_genes_RF['nM']), _rf(t.whole_2_genes_RF['nR'])] for t in tds),
            ['', 'BS'] + cc([_bs(t.whole_2_genes_BS['nM']), _bs(t.whole_2_genes_BS['nR'])] for t in tds),
            ['', 'KC'] + cc([_kc(t.whole_2_genes_KC['nM']), _kc(t.whole_2_genes_KC['nR'])] for t in tds),
            ['', 'KCT'] + cc([_kc(t.whole_2_genes_KCT['nM']), _kc(t.whole_2_genes_KCT['nR'])] for t in tds),
            #
            [],
            ['original-normalized'],
            ['', 'Distance'] + ['MrBayes', 'RAxML'] * nt,
            ['Non-partitioned', 'RF'] + cc([_rf(t.orig_2_norm_RF['WM']), _rf(t.orig_2_norm_RF['WR'])] for t in tds),
            ['', 'BS'] + cc([_bs(t.orig_2_norm_BS['WM']), _bs(t.orig_2_norm_BS['WR'])] for t in tds),
            ['', 'KC'] + cc([_kc(t.orig_2_norm_KC['WM']), _kc(t.orig_2_norm_KC['WR'])] for t in tds),
            ['', 'KCT'] + cc([_kc(t.orig_2_norm_KCT['WM']), _kc(t.orig_2_norm_KCT['WR'])] for t in tds),
            ['Partitioned', 'RF'] + cc([_rf(t.orig_2_norm_RF['GM']), _rf(t.orig_2_norm_RF['GR'])] for t in tds),
            ['', 'BS'] + cc([_bs(t.orig_2_norm_BS['GM']), _bs(t.orig_2_norm_BS['GR'])] for t in tds),
            ['', 'KC'] + cc([_kc(t.orig_2_norm_KC['GM']), _kc(t.orig_2_norm_KC['GR'])] for t in tds),
            ['', 'KCT'] + cc([_kc(t.orig_2_norm_KCT['GM']), _kc(t.orig_2_norm_KCT['GR'])] for t in tds),
        ])


class NormalizationResult:
    def __init__(self, project):
        self.project = project
        self.outgroup = None
        self.analyses_step = None
        self.tree_steps = dict()  # step_name -> step object
        self.trees = dict()       # step_name -> PhylogeneticTree object
        #
        self.G_tree_diffs = None  # Of type _TreeDiffs
        # self.N_tree_diffs = None

    def load_diffs_from_step(self, step_obj):
        data = step_obj.get_summary_data()
        assert data, step_obj.directory
        self.G_tree_diffs = _TreeDiffs(self, data=data['G'])
        # self.N_tree_diffs = _TreeDiffs(self, data=data['N'])
        return self

    def run(self, step_data):
        self._find_project_data()

        self.G_tree_diffs = _TreeDiffs(self, seq_type='G')
        # self.N_tree_diffs = _TreeDiffs(self, seq_type='N')

        rows = self.G_tree_diffs.get_rows()  # self.N_tree_diffs)
        text = str(StringColumns(rows)) + '\n'

        # Create step and collect data
        step = TableStep(self.project, step_data, remove_data=True)
        self.analyses_step.propagate_step_name_prefix(step)
        step.set_table_data(rows, [(f'c_{i}', 'str') for i in range(1, len(rows[0]) + 1)])
        step.save()  # completed=True
        #
        step.save_summary_data(dict(G=self.G_tree_diffs.to_dict()))  # , N=self.N_tree_diffs.to_dict()))
        write_str_in_file(step.step_file('summary.txt'), text)
        print(text)
        step.to_excel('normalization_result.xls', header=False)
        #
        self._create_graph(step)

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

    # One graph
    def create_graph(self, step_obj, show=False):
        if not step_obj.is_completed():
            raise ZCItoolsValueError(f"Normalization result input step is not completed!")
        self.load_diffs_from_step(step_obj)
        self._create_graph(show=show)

    def _create_graph(self, show=False):
        plt = self._create_one_graph(False)
        plt.savefig('tree_comparisons_no_labels.svg')
        plt.savefig('tree_comparisons_no_labels.png', dpi=150)
        plt.close()

        plt = self._create_one_graph(True)
        plt.savefig('tree_comparisons.svg')
        plt.savefig('tree_comparisons.png', dpi=150)

        if show:
            plt.show()

    def _create_one_graph(self, with_x_labels):
        plt = import_matplotlib_pylot()
        fig, ax = plt.subplots(figsize=(figsize_x, figsize_y), constrained_layout=True)
        # fig.subplots_adjust(right=(1 - 3 * right_axis_space), bottom=bottom_space)

        self._create_figure(fig, ax, with_x_labels, with_x_labels)
        return plt

    # More graph
    @classmethod
    def create_graphs(cls, project, step_objects, show=False):
        assert len(step_objects) > 1, len(step_objects)
        if any(not s.is_completed() for s in step_objects):
            raise ZCItoolsValueError(f"Some of normalization result input step are not completed!")
        nrs = [NormalizationResult(project).load_diffs_from_step(s) for s in step_objects]

        plt = import_matplotlib_pylot()
        # In inches. A4 is 8-1/4 * 11-3/4.
        # Note: figsize is value used to calibrate all other values
        fig, axes = plt.subplots(len(nrs), figsize=(figsize_x, figsize_y + figsize_y_no_x * (len(nrs) - 1)), constrained_layout=True)
        # # Adjust plot, because labels on X-axis and right Y axis
        # fig.subplots_adjust(right=(1 - 3 * right_axis_space), bottom=bottom_space / len(nrs), wspace=50)

        nrs[0]._create_figure(fig, axes[0], False, True)
        for nr, ax in zip(nrs[1:-1], axes[1:-1]):
            nr._create_figure(fig, ax, False, False)
        nrs[-1]._create_figure(fig, axes[-1], True, False)

        plt.savefig('all_tree_comparisons.svg')
        plt.savefig('all_tree_comparisons.png', dpi=150)

        if show:
            plt.show()

    # Figure
    def _create_figure(self, fig, ax, with_x_labels, with_legend):
        colors = dict(rf='blue', kct='orange', kc='red', bs='green')

        # Bars
        x_offset = 2  # Offset from x=0
        d_bar = 1     # Distance between neighbouring bars is 1
        d_gap_1 = 1
        d_gap_2 = 2
        d_gap_3 = 4
        bar_width = 0.9

        # Y tick
        tick_kw = dict(size=0, width=1)

        #
        _b4 = d_bar * 4
        g_starts = [0, _b4 + d_gap_1, 2 * _b4 + d_gap_1 + d_gap_2, 3 * _b4 + 2 * d_gap_1 + d_gap_2]
        group_gap = g_starts[-1] + _b4 + d_gap_3

        starts = dict()
        for group in range(3):
            gd = group * group_gap + x_offset
            for idx, x in enumerate(g_starts):
                x_group = x + gd
                for bar in range(4):
                    starts[(group, idx, bar)] = x_group + bar

        # Collect data into lists of tuples (x, y)
        gtd = self.G_tree_diffs
        rf = [(starts[(0, b, 0)], _rf_v(gtd.mr_bayes_2_raxml_RF[t])) for b, t in enumerate(('oW', 'nW', 'oG', 'nG'))]
        rf += [(starts[(1, b, 0)], _rf_v(gtd.whole_2_genes_RF[t])) for b, t in enumerate(('oM', 'oR', 'nM', 'nR'))]
        rf += [(starts[(2, b, 0)], _rf_v(gtd.orig_2_norm_RF[t])) for b, t in enumerate(('WM', 'GM', 'WR', 'GR'))]

        kct = [(starts[(0, b, 1)], gtd.mr_bayes_2_raxml_KCT[t]) for b, t in enumerate(('oW', 'nW', 'oG', 'nG'))]
        kct += [(starts[(1, b, 1)], gtd.whole_2_genes_KCT[t]) for b, t in enumerate(('oM', 'oR', 'nM', 'nR'))]
        kct += [(starts[(2, b, 1)], gtd.orig_2_norm_KCT[t]) for b, t in enumerate(('WM', 'GM', 'WR', 'GR'))]

        kc = [(starts[(1, b, 2)], gtd.whole_2_genes_KC[t]) for b, t in enumerate(('oM', 'oR', 'nM', 'nR'))]
        kc += [(starts[(2, b, 2)], gtd.orig_2_norm_KC[t]) for b, t in enumerate(('WM', 'GM', 'WR', 'GR'))]

        bs = [(starts[(1, b, 3)], _bs_v(gtd.whole_2_genes_BS[t])) for b, t in enumerate(('oM', 'oR', 'nM', 'nR'))]
        bs += [(starts[(2, b, 3)], _bs_v(gtd.orig_2_norm_BS[t])) for b, t in enumerate(('WM', 'GM', 'WR', 'GR'))]

        # Remove default labels
        ax.set_yticks([])
        ax.set_xticks([])
        for side in ('left', 'right', 'top'):
            ax.spines[side].set_visible(False)

        #
        rf_ax = ax.twinx()
        kct_ax = ax.twinx()
        kc_ax = ax.twinx()
        bs_ax = ax.twinx()

        # Offset the right spines
        kct_ax.spines['right'].set_position(('axes', 1 + right_axis_space))
        kc_ax.spines['right'].set_position(('axes', 1 + 2 * right_axis_space))
        bs_ax.spines['right'].set_position(('axes', 1 + 3 * right_axis_space))

        lines = []
        for _ax, vals, label, c in ((rf_ax, rf, 'RF', colors['rf']),
                                    (kct_ax, kct, 'KCT', colors['kct']),
                                    (kc_ax, kc, 'KC', colors['kc']),
                                    (bs_ax, bs, 'BS', colors['bs'])):
            fix_patch_spines(_ax)
            #
            mv = _max_val(max(y for x, y in vals))
            _ax.set_ylim(0, mv)
            _ax.set_yticks([mv])
            _ax.set_yticklabels(_ax.get_yticks(), rotation=90)
            _ax.set_ylabel(label, labelpad=-10)  # font size? loc='bottom',
            _ax.yaxis.label.set_color(c)
            _ax.spines['right'].set_color(c)
            _ax.tick_params(axis='y', colors=c, **tick_kw)
            #
            lines.append(_ax.bar([x for x, y in vals], [y for x, y in vals], bar_width, label=label, color=c))

        if with_legend:
            ax.legend(lines, [l.get_label() for l in lines], loc='upper left', frameon=False)

        if with_x_labels:
            # Add labels on X axis
            l_o = d_bar * 0.5
            x_labels_1 = [l_o + starts[(group, idx, 1)] for group, idx in product(range(3), range(4))]
            x_labels_2 = [(a + b) / 2 for a, b in zip(x_labels_1[::2], x_labels_1[1::2])]
            x_labels_3 = [(a + b) / 2 for a, b in zip(x_labels_2[::2], x_labels_2[1::2])]
            for idx, label in enumerate(('O', 'N', 'O', 'N', 'BI', 'ML', 'BI', 'ML', 'NP', 'P', 'NP', 'P')):
                ax.text(x_labels_1[idx], x_label_y[0], label, ha='center', va='top', fontsize=8)
            for idx, label in enumerate(('NP', 'P', 'O', 'N', 'BI', 'ML')):
                ax.text(x_labels_2[idx], x_label_y[1], label, ha='center', va='top', fontsize=10)
            for idx, label in enumerate(('BI-ML', 'NP-P', 'O-N')):
                ax.text(x_labels_3[idx], x_label_y[2], label, ha='center', va='top', fontsize=12)


def _max_val(value):
    if value == 0:
        return 1
    nd_1 = int(floor(log10(abs(value))))
    f_d = value / 10 ** nd_1
    if f_d > 5:
        return 10 ** (nd_1 + 1)
    if ceil(f_d) > 2:
        return ceil(f_d) * 10 ** nd_1
    if f_d > 1.5:
        return 2 * 10 ** nd_1
    return 1.5 * 10 ** nd_1


def fix_patch_spines(ax):
    # Having been created by twinx, par2 has its frame off, so the line of its
    # detached spine is invisible. First, activate the frame but make the patch
    # and spines invisible.
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)
    # Second, show the right spine.
    ax.spines['right'].set_visible(True)
