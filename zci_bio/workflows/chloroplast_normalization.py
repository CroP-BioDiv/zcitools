from step_project.base_workflow import BaseWorkflow
from common_utils.cache import cache


class ChloroplastNormalization(BaseWorkflow):
    _WORKFLOW = 'chloroplast_normalization'

    @staticmethod
    def required_parameters():
        return ('family', 'outgroup')

    @cache
    def _actions(self):
        family = self.parameters['family']
        outgroup = self.parameters['outgroup']

        actions = [
            ('01_chloroplast_list', f"ncbi_chloroplast_list -f {family} -o {outgroup}"),
            ('02_seqs', 'fetch_seqs 01_chloroplast_list'),
            ('03_GeSeq', 'ge_seq 02_seqs'),
            ('04_AnalyseChloroplast', 'analyse_chloroplast 03_GeSeq'),
            #
            ('nS_01_seq', 'fix_by_analysis parts 04_AnalyseChloroplast'),
            ('nA_01_seq', 'fix_by_analysis parts 04_AnalyseChloroplast -a'),
            ('nS_02_GeSeq', 'ge_seq nS_01_seq'),
            ('nA_02_GeSeq', 'ge_seq nA_01_seq', 'nS_02_GeSeq'),  # Same sequences are annotated, no need to have 2 finishes pending!
            #
            ('oS_02_GeSeq', 'seq_subset 03_GeSeq --analyses-with-irs 04_AnalyseChloroplast'),
            ('oA_02_GeSeq', 'seq_subset 03_GeSeq'),  # All
        ]

        # +n{S|A}_03_MAFFT MAFFT
        # n{S|A}_og_MAFFT MAFFT
        # +-n{S|A}_04_MrBayes_C MrBayes
        # +-n{S|A}_04_MrBayes_P MrBayes
        # +-n{S|A}_04_RAxML_C   RAxML
        # +-n{S|A}_04_RAxML_P   RAxML
        # n{S|A}_05_trees Analiza

        # +o{S|A}_03_MAFFT MAFFT
        # o{S|A}_og_MAFFT MAFFT
        # +-o{S|A}_04_MrBayes_C MrBayes
        # +-o{S|A}_04_MrBayes_P MrBayes
        # +-o{S|A}_04_RAxML_C   RAxML
        # +-o{S|A}_04_RAxML_P   RAxML
        # o{S|A}_05_trees Analiza

        for a in ('n', 'o'):
            for b in ('S', 'A'):
                ab = f'{a}{b}'
                s2 = f'{ab}_02_GeSeq'
                s3 = f'{ab}_03_MAFFT'
                s4_b_c = f'{ab}_04_MrBayes_C'
                s4_b_p = f'{ab}_04_MrBayes_P'
                s4_r_c = f'{ab}_04_RAxML_C'
                s4_r_p = f'{ab}_04_RAxML_P'
                phylos = [s4_b_c, s4_b_p, s4_r_c, s4_r_p]
                s5 = f'{ab}_05_trees'
                actions.extend([
                    (s3, f'align_genomes {s2} w'),
                    # (s4_b_c, ['mr_bayes', s3]),  # , '?'
                    # (s4_b_p, ['mr_bayes', s3]),  # , '?'
                    # (s4_r_c, ['raxml', s3]),  # , '?'
                    # (s4_r_p, ['raxml', s3]),  # , '?'
                    # (s5, ['??'] + phylos),
                ])

        return actions
