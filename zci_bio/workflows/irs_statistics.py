import datetime
import itertools
from collections import defaultdict
from step_project.base_workflow import BaseWorkflow
from common_utils.exceptions import ZCItoolsValueError


class IRsStatistics(BaseWorkflow):
    _WORKFLOW = 'irs_statistics'

    @staticmethod
    def required_parameters():
        # , 'plastids'
        return ('taxons', 'methods')

    @staticmethod
    def format_parameters(params):
        # Methods
        from ..chloroplast.irs.analyse_irs import METHOD_NAMES
        methods = set()
        not_known_methods = set()
        for m in params['methods'].lower().split(','):
            if m == 'all':
                methods.update(METHOD_NAMES)
            elif m in METHOD_NAMES:
                methods.add(m)
            else:
                not_known_methods.add(m)
        if not_known_methods:
            raise ZCItoolsValueError(f'Not know method(s): {", ".join(not_known_methods)}!')
        params['methods'] = methods
        #
        params['plastids'] = int(params.get('plastids', 0))
        params['remove_irl'] = int(params.get('remove_irl', 0))
        if taxa_ranks := params.get('taxa_ranks'):
            params['taxa_ranks'] = taxa_ranks.split(',')
        if 'max_update_date' in params:
            params['max_update_date'] = datetime.date.fromisoformat(params['max_update_date'])
        return params

    def _actions(self):
        from ..chloroplast.irs.analyse_irs import METHODS_USE_SEQUENCES, METHODS_SEPARATE_PATH

        taxons = ' '.join(f'-t {t}' for t in self.parameters['taxons'].split(','))
        plastids = '-P' if self.parameters['plastids'] else ''
        if max_update_date := self.parameters.get('max_update_date', ''):
            max_update_date = f'--max-update-date {max_update_date}'
        remove_irl = '--remove-irl' if int(self.parameters.get('remove_irl', '0')) else ''
        methods = self.parameters['methods']

        actions = [('01_chloroplast_list', f"ncbi_chloroplast_list {taxons} {plastids} {max_update_date} {remove_irl}")]

        # Collect data
        # Methods that use NCBI sequences from common step
        stats = [f'-m {m}' for m in methods]
        seqs_methods = [m for m in methods if m in METHODS_USE_SEQUENCES]
        seqs_methods = ' '.join(f'-s {m}' for m in seqs_methods)
        actions.append(('02_seqs', f"analyse_irs_collect_needed_data 01_chloroplast_list seqs {seqs_methods}"))

        # Methods that use separate path to collect data
        for m in methods:
            if m in METHODS_SEPARATE_PATH:
                stats.append(f'-{m[0]} 03_{m}')
                actions.append((f'02_{m}', f"analyse_irs_collect_needed_data 01_chloroplast_list {m}"))
                actions.append((f'03_{m}', f"{m} 02_{m}"))

        # Analysis
        cmd = f"analyse_irs 01_chloroplast_list 02_seqs {' '.join(stats)}"
        if ranks := self.parameters.get('taxa_ranks'):
            cmd += ' ' + ' '.join(f'-r {r}' for r in ranks)
        if names := self.parameters.get('taxons'):
            cmd += ' ' + ' '.join(f'-n {n}' for n in names.split(','))
        actions.append(('04_stats', cmd))

        # Summary, result, ...
        return actions

    def get_summary(self):
        # ---------------------------------------------------------------------
        # Collect sequence data
        # ---------------------------------------------------------------------
        if not (step_01 := self.project.read_step_if_in('01_chloroplast_list')):
            return dict(text='Project not started!')

        if not (results := self.project.read_step_if_in('04_stats')):
            return dict(text='Analysis not done!')

        #
        clade_2_num = defaultdict(int)
        clade_2_ns = defaultdict(int)
        clade_2_family = defaultdict(set)
        clade_2_genus = defaultdict(set)
        clade_2_species = defaultdict(set)
        clade_2_year = defaultdict(lambda: defaultdict(int))
        #
        methods = self.parameters['methods']
        it_2_idx = dict(exact=0, differs=1, no=2)
        m_2_it = dict((m, [0, 0, 0]) for m in methods)
        m_2_wraps = dict((m, [0, 0]) for m in methods)
        #
        m_2_dl, dl_splits = _idx_data(methods, (0, 10, 100))
        # m_2_max_dl = defaultdict(int)
        # m_2_blocks, block_splits = _idx_data(methods, (1, 10, 100))
        m_2_ddd, ddd_splits = _idx_data(methods, (5, 20, 100))
        #
        m_2_dna = dict((m, [0, 0, 0]) for m in methods)
        m_2_dna_irs = dict((m, [0, 0]) for m in methods)  # Can't be "No IRs"!
        #
        for clade, family, genus, species, published, \
                method, ir_type, ir_wraps, diff_len, \
                replace_num, replace_sum, indel_num, indel_sum, not_dna, not_dna_irs in \
                results.select(('Clade', 'family', 'genus', 'Organism', 'Published',
                                'Method', 'IR_type', 'IR_wraps', 'diff_len',
                                'replace_num', 'replace_sum', 'indel_num', 'indel_sum', 'not_dna', 'not_dna_irs')):
            if method == methods[0]:
                clade_2_num[clade] += 1
                clade_2_year[clade][published.year] += 1
                if not_dna:
                    clade_2_ns[clade] += 1
                clade_2_family[clade].add(family)
                clade_2_genus[clade].add(genus)
                clade_2_species[clade].add(species)
            ir_type_idx = it_2_idx[ir_type]
            m_2_it[method][ir_type_idx] += 1
            if ir_type != 'no':       # With IRs
                m_2_wraps[method][int(ir_wraps)] += 1
                m_2_dl[method][_idx(diff_len, dl_splits)] += 1
                # if diff_len < 20000:
                #     m_2_max_dl[method] = max(m_2_max_dl.get(method, 0), diff_len)
            if ir_type == 'differs':  # Not exact IRs
                # m_2_blocks[method][_idx(replace_num + indel_num, block_splits)] += 1
                m_2_ddd[method][_idx(replace_sum + indel_sum, ddd_splits)] += 1
            if not_dna:
                m_2_dna[method][ir_type_idx] += 1
                if not_dna_irs:
                    m_2_dna_irs[method][ir_type_idx] += 1

        #
        clades = dict()
        for c, n in sorted(clade_2_num.items()):
            clades[c] = [n, len(clade_2_family[c]), len(clade_2_genus[c]), len(clade_2_species[c]), clade_2_ns[c]]
        clades['All'] = [sum(v[i] for v in clades.values()) for i in range(5)]
        years = sorted(set(itertools.chain.from_iterable(v.keys() for v in clade_2_year.values())))
        clade_2_year['All'] = dict((y, sum(v.get(y, 0) for v in clade_2_year.values())) for y in years)
        clade_2_year = dict((k, [v[y] for y in years]) for k, v in clade_2_year.items())
        clade_2_year['Cumulative'] = list(itertools.accumulate(clade_2_year['All']))

        #
        ir_types = ('Exact IRs', 'IRs differ', 'No IRs')
        text = "Summary\n"
        # text += self._dict_table(clade_2_num, 'Number of Sequences')
        # text += self._dict_table(clade_2_family, 'Number of Families')
        text += self._methods_table(clades, 'Number of Sequences', ('Sequences', 'Families', 'Genus', 'Species', "With N's"), methods=list(clades.keys()))
        text += self._methods_table(clade_2_year, 'Number of Sequences per year', years, methods=list(clade_2_year.keys()))

        text += "\n\nAnnotated IRs, to sequence characteristics"
        text += self._methods_table(m_2_it, 'Number of Sequences with Annotated IRs', ir_types, ident='  ')
        text += self._methods_table(m_2_dna, "N's in sequence", ir_types, ident='  ')
        text += self._methods_table(m_2_dna_irs, "N's in IRs", ir_types[:-1], ident='  ')
        text += "\n\nFound IRs characteristics"
        text += self._methods_table(m_2_wraps, 'IR wraps', ('No', 'Yes'), ident='  ')
        text += self._methods_table(m_2_dl, 'IR difference in length', _idx_labels(dl_splits), ident='  ')  # , measure=' bp'
        # text += self._methods_table(dict((k, [v]) for k, v in m_2_max_dl.items()), 'Max IR difference in length', ('Max',), ident='  ')  # , measure=' bp'
        # text += self._methods_table(m_2_blocks, 'IR indel/replace number of blocks', _idx_labels(block_splits), ident='  ')
        text += self._methods_table(m_2_ddd, 'IR indel/replace number of bps', _idx_labels(ddd_splits), ident='  ', percentages=True)
        return dict(text=text)

    def _dict_table(self, _dict, title):
        lines = '\n'.join(f'  {k:<10} {n:>5}' for k, n in _dict.items())
        return f"\n{title}:\n{lines}\n  All {sum(_dict.values()):>12}\n"

    def _methods_table(self, data, title, labels, methods=None, ident='', percentages=False):
        if not methods:
            methods = self.parameters['methods']
        max_l = max(len(m) for m in methods)
        text = f"""\n{ident}{title}:
             {ident}{' '.join(m.rjust(max_l) for m in methods)}
"""
        if percentages:
            _sum = dict((m, sum(data[m][idx] for idx in range(len(labels)))) for m in methods)
            for idx, lab in enumerate(labels):
                text += f"  {ident}{lab:<10} {' '.join(str(data[m][idx]).rjust(max_l) for m in methods)}\n"
                text += f"  {ident}{' ' * 10} {' '.join(str(round(100 * data[m][idx] / _sum[m], 2)).rjust(max_l) for m in methods)}\n"
        else:
            for idx, lab in enumerate(labels):
                text += f"  {ident}{lab:<10} {' '.join(str(data[m][idx]).rjust(max_l) for m in methods)}\n"
        return text


def _idx_data(methods, lens):
    n = len(lens) + 1
    return dict((m, [0] * n) for m in methods), lens


def _idx(x, lens):
    for idx, l in enumerate(lens):
        if x <= l:
            return idx
    return len(lens)


def _idx_labels(lens, measure=''):
    lens = [str(l) for l in lens]
    max_l = len(lens[-1])
    labels = []
    for l in lens:
        if l == '0':
            labels.append(f'   {l.rjust(max_l)}')
        else:
            labels.append(f'<= {l.rjust(max_l)}')
    labels.append(f'>  {lens[-1]}')
    return [l + measure for l in labels] if measure else labels
