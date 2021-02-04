import os.path
from step_project.base_workflow import BaseWorkflow
from common_utils.file_utils import read_csv
from common_utils.exceptions import ZCItoolsValueError


def workflow_branches(wf_settings, analyses_step):
    branches = ['G']
    # branches = ['G', 'N']
    if not analyses_step:
        return branches
    # parts = analyses_step.select(['GeSeq part starts', 'NCBI part starts'])
    # if any(bool(g) != bool(n) for g, n in parts):
    #     branches.append('S')
    # if int(wf_settings.get('calc_all', 0)) and not all(g or n for g, n in parts):
    #     branches.append('A')
    return branches


class ChloroplastNormalization(BaseWorkflow):
    _WORKFLOW = 'chloroplast_normalization'

    @staticmethod
    def required_parameters():
        return ('family', 'outgroup')

    def _actions(self):
        family = self.parameters['family']
        outgroup = self.parameters.get('outgroup')
        analyses_branches = workflow_branches(
            self.parameters, self.project.read_step_if_in('04_AnalyseChloroplast', check_data_type='table'))
        o_o = f' -o {outgroup}' if outgroup else ''

        actions = [
            ('01_chloroplast_list', f"ncbi_chloroplast_list -f {family}{o_o}"),
            ('02_seqs', 'fetch_seqs 01_chloroplast_list'),
            ('03_GeSeq', 'ge_seq 02_seqs'),
            ('04_AnalyseChloroplast', 'analyse_chloroplast 03_GeSeq'),

            # G(eSeq) branch
            ('Gn_01_seq', 'fix_by_analysis parts ge_seq 04_AnalyseChloroplast'),
            ('Gn_02_GeSeq', 'ge_seq Gn_01_seq'),
            #
            ('Go_02_GeSeq', 'seq_subset 03_GeSeq --analyses-with-irs 04_AnalyseChloroplast --analyses-subset ge_seq'),

            # N(CBI) branch
            ('Nn_01_seq', 'fix_by_analysis parts ncbi 04_AnalyseChloroplast'),
            ('Nn_02_GeSeq', 'ge_seq Nn_01_seq'),
            #
            ('No_02_GeSeq', 'seq_subset 03_GeSeq --analyses-with-irs 04_AnalyseChloroplast --analyses-subset ncbi'),
        ]
        if 'S' in analyses_branches:
            actions += [
                ('Sn_01_seq', 'fix_by_analysis parts sum 04_AnalyseChloroplast'),
                # Note: Additional dependency added since same sequences are annotated,
                # and there is no need to have 2 steps pending for finishing.
                ('Sn_02_GeSeq', 'ge_seq Sn_01_seq', 'Gn_02_GeSeq'),
                #
                ('So_02_GeSeq', 'seq_subset 03_GeSeq --analyses-with-irs 04_AnalyseChloroplast --analyses-subset sum'),
            ]
        if 'A' in analyses_branches:
            actions += [
                ('An_01_seq', 'fix_by_analysis parts all 04_AnalyseChloroplast'),
                # Note: Additional dependency added since same sequences are annotated,
                # and there is no need to have 2 steps pending for finishing.
                #
                ('An_02_GeSeq', 'ge_seq An_01_seq', 'Sn_02_GeSeq'),
                ('Ao_02_GeSeq', 'seq_subset 03_GeSeq'),  # All
            ]

        tree_steps = []
        for sa in analyses_branches:
            for on in 'on':
                ab = f'{sa}{on}'
                s2 = f'{ab}_02_GeSeq'
                s3 = f'{ab}_03_MAFFT'
                s4_c_b = f'{ab}_04_W_MrBayes'
                s4_c_r = f'{ab}_04_W_RAxML'
                s4_p_b = f'{ab}_04_G_MrBayes'
                s4_p_r = f'{ab}_04_G_RAxML'
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

    def get_summary(self):
        text = ''

        # ---------------------------------------------------------------------
        # Collect sequence data
        # ---------------------------------------------------------------------
        if not (step := self.project.read_step_if_in('01_chloroplast_list')):
            return text

        s = step.get_summary_data()
        text += f"""
# Input data

Family:
 * Number of family genomes in NCBI : {s['ncbi_family_num_genomes']}
 * Number of family genomes to work : {s['family_num_genomes']}
 * Same species description         : {s['same_species_description']}

Outgroup:
 * Accession number  : {s['outgroup_accession']}
 * Species           : {s['outgroup_species']}
 * Sequence length   : {s['outgroup_max_length']:,}
 * Create date       : {s['outgroup_max_create_date']}
 * Update date       : {s['outgroup_max_update_date']}

Set of genome to work:
 * Length range family : {s['family_min_length']:,} - {s['family_max_length']:,}
 * Length range all    : {s['all_min_length']:,} - {s['all_max_length']:,}
 * Date first range    : {s['all_min_create_date']} - {s['all_max_create_date']}
 * Date update range   : {s['all_min_update_date']} - {s['all_max_update_date']}
"""

        # ---------------------------------------------------------------------
        # Analyses
        # ---------------------------------------------------------------------
        if not (analyses_step := self.project.read_step_if_in('04_AnalyseChloroplast')):
            return text

        s = analyses_step.get_summary_data()
        ng = s['num_genomes']

        def _annotation_summary(title, a, genes):
            w_offset = s[f'{a}_wrong_offset_list'] or []
            if w_offset:
                w_offset_mm = f'{min(map(abs, w_offset))} - {max(map(abs, w_offset))} '
            else:
                w_offset_mm = ''

            t = f"""
{title}:
 * Number of annotations containing IRS         : {s[f'{a}_num_irs']}
 * Number of annotations without IRS            : {ng - s[f'{a}_num_irs']}
 * Number of annotations with wrong orientation : {s[f'{a}_wrong_orientations']}
 * Description wrong orientations               : {s[f'{a}_desc_wrong_orientations']}
 * Number of annotations with offset            : {s[f'{a}_wrong_offset']}
 * List of wrong offsets                        : {w_offset_mm}[{", ".join(map(str, w_offset))}]
 * Number genomes to normalize                  : {s[f'{a}_num_to_fix']}
"""
            if genes:
                t += f"""
    Range of number of genes:
     * Annotated        : {s[f'{a}_genes_min_annotated']} - {s[f'{a}_genes_max_annotated']}
     * Disjunct         : {s[f'{a}_genes_min_disjunct']} - {s[f'{a}_genes_max_disjunct']}
     * Name/strand      : {s[f'{a}_genes_min_name_strand']} - {s[f'{a}_genes_max_name_strand']}
     * Name             : {s[f'{a}_genes_min_names']} - {s[f'{a}_genes_max_names']}
     * Without location : {s[f'{a}_genes_min_without_location']} - {s[f'{a}_genes_max_without_location']}
     * Without name     : {s[f'{a}_genes_min_without_name']} - {s[f'{a}_genes_max_without_name']}

    Range of part lengths:
     * LSC : {s[f'{a}_part_min_length_lsc']} - {s[f'{a}_part_max_length_lsc']}
     * IRs : {s[f'{a}_part_min_length_ir']} - {s[f'{a}_part_max_length_ir']}
     * SSC : {s[f'{a}_part_min_length_ssc']} - {s[f'{a}_part_max_length_ssc']}
"""
            return t

        text += f"""

# Analyses

Number of genomes: {ng}

{_annotation_summary('GeSeq annotations', 'ge_seq', True)}
{_annotation_summary('NCBI annotations', 'ncbi', True)}
{_annotation_summary('Sum of annotations', 'sum', False)}
"""

        # ---------------------------------------------------------------------
        # Alignments
        # ---------------------------------------------------------------------
        analyses_branches = workflow_branches(self.parameters, analyses_step)
        gas_2_text = dict(G='GeSeq', S='Sum', A='All', N='NCBI')
        text += f"""

# Alignment

Subset        Original      Normalized
"""

        for b in analyses_branches:
            aligns = []
            for on in 'on':
                alignment_length = '-'
                if step := self.project.read_step_if_in(f'{b}{on}_03_MAFFT'):
                    if s := step.get_summary_data():
                        # ToDo: is it safe way to find number of partitions?
                        num_parts = ''
                        part_f = os.path.join(f'{b}{on}_04_G_RAxML', 'partitions.ind')
                        if os.path.isfile(part_f):
                            num_parts = f' ({sum(1 for l in open(part_f))})'
                        #
                        alignment_length = f"{s['alignment_length']:,}{num_parts}"
                aligns.append(alignment_length)
            text += f"{gas_2_text[b]:6} {aligns[0]:>15} {aligns[1]:>15}\n"

        # ToDo: Trees

        return dict(text=text)
