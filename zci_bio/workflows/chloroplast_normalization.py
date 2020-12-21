from step_project.base_workflow import BaseWorkflow, WfAction
from common_utils.cache import cache


class ChloroplastNormalization(BaseWorkflow):
    _WORKFLOW = 'chloroplast_normalization'

    @staticmethod
    def required_parameters():
        return ('family', 'outgroup')

    @cache
    def actions(self):
        actions = [
            WfAction('01_chloroplast_list', None, ['ncbi_chloroplast_list', '-f', self.parameters['family']]),
            WfAction('02_seqs', '01_chloroplast_list', ['fetch_seqs', '01_chloroplast_list']),
            WfAction('03_GeSeq', '02_seqs', ['ge_seq', '02_seqs']),
            WfAction('04_AnalyseChloroplast', '03_GeSeq', ['analyse_chloroplast', '03_GeSeq']),
            #
            # WfAction('nS_01_seq', '04_AnalyseChloroplast', ['??', '04_AnalyseChloroplast']),
            # WfAction('nA_01_seq', '04_AnalyseChloroplast', ['??', '04_AnalyseChloroplast']),
            # WfAction('nS_02_GeSeq', 'nS_01_seq', ['ge_seq', 'nS_01_seq']),
            # WfAction('nA_02_GeSeq', 'nA_01_seq', ['ge_seq', 'nA_01_seq']),
            #
            # WfAction('aS_02_GeSeq', '04_AnalyseChloroplast', ['??', '04_AnalyseChloroplast']),
            # WfAction('aA_02_GeSeq', '04_AnalyseChloroplast', ['??', '04_AnalyseChloroplast']),
        ]

        # +n{S|A}_01_seqs  Obrada
        # +n{S|A}_02_GeSeq GeSeq obrada
        # +n{S|A}_03_MAFFT MAFFT
        # n{S|A}_og_MAFFT MAFFT
        # +-n{S|A}_04_MrBayes_C MrBayes
        # +-n{S|A}_04_MrBayes_P MrBayes
        # +-n{S|A}_04_RAxML_C   RAxML
        # +-n{S|A}_04_RAxML_P   RAxML
        # n{S|A}_05_trees Analiza

        # +a{S|A}_02_GeSeq Filtriranje
        # +a{S|A}_03_MAFFT MAFFT
        # a{S|A}_og_MAFFT MAFFT
        # +-a{S|A}_04_MrBayes_C MrBayes
        # +-a{S|A}_04_MrBayes_P MrBayes
        # +-a{S|A}_04_RAxML_C   RAxML
        # +-a{S|A}_04_RAxML_P   RAxML
        # a{S|A}_05_trees Analiza

        # for a in ('n', 'a'):
        #     for b in ('S', 'A'):
        #         ab = f'{a}{b}'
        #         s2 = f'{ab}_02_GeSeq'
        #         s3 = f'{ab}_03_MAFFT'
        #         s4_b_c = f'{ab}_04_MrBayes_C'
        #         s4_b_p = f'{ab}_04_MrBayes_P'
        #         s4_r_c = f'{ab}_04_RAxML_C'
        #         s4_r_p = f'{ab}_04_RAxML_P'
        #         phylos = [s4_b_c, s4_b_p, s4_r_c, s4_r_p]
        #         s5 = f'{ab}_05_trees'
        #         actions.extend([
        #             WfAction(s3, s2, ['align_genomes', s2]),
        #             # WfAction(s4_b_c, s3, ['mr_bayes', s3, '?']),
        #             # WfAction(s4_b_p, s3, ['mr_bayes', s3, '?']),
        #             # WfAction(s4_r_c, s3, ['raxml', s3, '?']),
        #             # WfAction(s4_r_p, s3, ['raxml', s3, '?']),
        #             # WfAction(s5, phylos, ['??'] + phylos),
        #         ])

        return self._check_actions(actions)
