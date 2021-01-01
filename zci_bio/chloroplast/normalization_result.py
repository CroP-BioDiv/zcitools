from step_project.common.table.steps import TableStep
from common_utils.file_utils import get_settings
from common_utils.exceptions import ZCItoolsValueError


class NormalizationResult:
    def __init__(self, project):
        self.project = project
        self.outgroup = None
        self.analyses_step = None
        self.trees = dict()            # step_name -> step object
        self.phylos_on_same_data = []  # Pairs of strings (MrBayes step name, RAxML step name)

    def get_rooted_tree(self, step_name, library):
        tree = self.trees[step_name].get_consensus_tree(library)
        if not (leaves := tree.get_leaves_by_name(self.outgroup)):
            raise ZCItoolsValueError(f"No outgroup in step's tree! {step_name}")
        if len(leaves) != 1:
            raise ZCItoolsValueError(f"More nodes with name {outgroup} in step {step_name}!")
        tree.set_outgroup(leaves[0])
        return tree

    #
    def run(self, step_data):
        self._find_project_data()

        #
        for step_name_mr_bayes, step_name_raxml in self.phylos_on_same_data:
            tree_mr_bayes = self.get_rooted_tree(step_name_mr_bayes, 'ete')
            tree_raxml = self.get_rooted_tree(step_name_raxml, 'ete')
            d = tree_mr_bayes.compare(tree_raxml)
            # Keys: rf, max_rf, ref_edges_in_source, source_edges_in_ref, effective_tree_size,
            #       norm_rf, treeko_dist, source_subtrees, common_edges, source_edges, ref_edges

            print(step_name_mr_bayes, step_name_raxml, [d[x] for x in ('rf', 'max_rf', 'norm_rf')])
            if d['rf']:
                print('  ', d['source_edges'] - d['common_edges'])
                print('  ', d['ref_edges'] - d['common_edges'])

        # Create step and collect data
        step = TableStep(self.project, step_data, remove_data=True)
        self.analyses_step.propagate_step_name_prefix(step)
        step.save()
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
        has_A = not all(self.analyses_step.get_column_values('Part starts'))

        # Find all phylogenetic steps
        for a in ('n', 'o'):
            for b in (('S', 'A') if has_A else ('S',)):
                ab = f'{a}{b}_04'
                for parts in ('C', 'P'):
                    ab_part = f'{ab}_{parts}'
                    for phylo, dt in (('MrBayes', 'mr_bayes'), ('RAxML', 'raxml')):
                        step_name = f'{ab_part}_{phylo}'
                        self.trees[step_name] = self.project.read_step_if_in(
                            step_name, check_data_type=dt, no_check=True)
                    self.phylos_on_same_data.append((f'{ab_part}_MrBayes', f'{ab_part}_RAxML'))

        if (no_steps := [k for k, v in self.trees.items() if not v]):
            raise ZCItoolsValueError(f"No tree step(s): {', '.join(sorted(no_steps))}!")

        # o{S|A}_04_C_MrBayes, o{S|A}_04_C_RAxML
        # o{S|A}_04_P_MrBayes, o{S|A}_04_P_RAxML
        # n{S|A}_04_C_MrBayes, n{S|A}_04_C_RAxML
        # n{S|A}_04_P_MrBayes, n{S|A}_04_P_RAxML
