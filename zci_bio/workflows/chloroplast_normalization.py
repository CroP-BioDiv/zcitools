from step_project.base_workflow import BaseWorkflow
from common_utils.file_utils import read_csv
from common_utils.exceptions import ZCItoolsValueError


class ChloroplastNormalization(BaseWorkflow):
    _WORKFLOW = 'chloroplast_normalization'
    _SUMMARY_STEPS = ['04_AnalyseChloroplast', '01_chloroplast_list']

    @staticmethod
    def required_parameters():
        return ('family', 'outgroup')

    def has_A_branch(self):
        if step := self.project.read_step_if_in('04_AnalyseChloroplast', check_data_type='table'):
            if all(step.get_column_values('Part starts')):
                return False
        return True

    def _actions(self):
        family = self.parameters['family']
        outgroup = self.parameters['outgroup']
        has_A = self.has_A_branch()

        actions = [
            ('01_chloroplast_list', f"ncbi_chloroplast_list -f {family} -o {outgroup}"),
            ('02_seqs', 'fetch_seqs 01_chloroplast_list'),
            ('03_GeSeq', 'ge_seq 02_seqs'),
            ('04_AnalyseChloroplast', 'analyse_chloroplast 03_GeSeq'),
            #
            ('nS_01_seq', 'fix_by_analysis parts 04_AnalyseChloroplast'),
            ('nS_02_GeSeq', 'ge_seq nS_01_seq'),
            #
            ('oS_02_GeSeq', 'seq_subset 03_GeSeq --analyses-with-irs 04_AnalyseChloroplast'),
        ]
        if has_A:
            actions += [
                ('nA_01_seq', 'fix_by_analysis parts 04_AnalyseChloroplast -a'),
                # Note: Additional dependency added since same sequences are annotated,
                # and there is no need to have 2 steps pending for finishing.
                ('nA_02_GeSeq', 'ge_seq nA_01_seq', 'nS_02_GeSeq'),
                ('oA_02_GeSeq', 'seq_subset 03_GeSeq'),  # All
            ]

        # n{S|A}_wo_MAFFT MAFFT
        # n{S|A}_05_trees Analiza

        # o{S|A}_wo_MAFFT MAFFT
        # o{S|A}_05_trees Analiza

        tree_steps = []
        for a in ('n', 'o'):
            for b in (('S', 'A') if has_A else ('S',)):
                ab = f'{a}{b}'
                s2 = f'{ab}_02_GeSeq'
                s3 = f'{ab}_03_MAFFT'
                s4_c_b = f'{ab}_04_C_MrBayes'
                s4_c_r = f'{ab}_04_C_RAxML'
                s4_p_b = f'{ab}_04_P_MrBayes'
                s4_p_r = f'{ab}_04_P_RAxML'
                tree_steps.extend([s4_c_b, s4_p_b, s4_c_r, s4_p_r])
                actions.extend([
                    (s3, f'align_genomes {s2} w'),
                    (s4_c_b, f'mr_bayes {s3} -p'),
                    (s4_c_r, f'raxml {s3} -p'),
                    (s4_p_b, f'mr_bayes {s3} -a {s2}'),
                    (s4_p_r, f'raxml {s3} -a {s2}'),
                ])
        actions.append(('X_result', 'normalization_result', tree_steps))
        return actions

    def summary(self):
        family = self.parameters['family']
        outgroup = self.parameters['outgroup']
        text = ''

        if step := self.project.read_step_if_in('01_chloroplast_list'):
            num_acc = sum(int(s != outgroup) for s in step.get_column_values('ncbi_ident'))
            num_removed = 0
            if (ss_csv := step.step_file('same_species.csv')):
                data = read_csv(ss_csv)
                num_removed = len(data) - len(set(d[0] for d in data))
            text += f"""
Number of genomes in NCBI : {num_acc + num_removed}
Number of genomes to work : {num_acc}
"""

        return text
