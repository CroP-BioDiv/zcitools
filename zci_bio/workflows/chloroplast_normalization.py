from step_project.base_workflow import BaseWorkflow
from common_utils.file_utils import read_csv
from common_utils.exceptions import ZCItoolsValueError


class ChloroplastNormalization(BaseWorkflow):
    _WORKFLOW = 'chloroplast_normalization'

    @staticmethod
    def required_parameters():
        return ('family', 'outgroup')

    def has_A_branch(self):
        if not self.parameters.get('calc_all', 0):
            return False
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
            ('Sn_01_seq', 'fix_by_analysis parts 04_AnalyseChloroplast'),
            ('Sn_02_GeSeq', 'ge_seq Sn_01_seq'),
            #
            ('So_02_GeSeq', 'seq_subset 03_GeSeq --analyses-with-irs 04_AnalyseChloroplast'),
        ]
        if has_A:
            actions += [
                ('An_01_seq', 'fix_by_analysis parts 04_AnalyseChloroplast -a'),
                # Note: Additional dependency added since same sequences are annotated,
                # and there is no need to have 2 steps pending for finishing.
                ('An_02_GeSeq', 'ge_seq An_01_seq', 'Sn_02_GeSeq'),
                ('Ao_02_GeSeq', 'seq_subset 03_GeSeq'),  # All
            ]

        tree_steps = []
        for sa in ('SA' if has_A else 'S'):
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
        family = self.parameters['family']
        outgroup = self.parameters['outgroup']
        text = ''

        if step := self.project.read_step_if_in('01_chloroplast_list'):
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
 * Length range      : {s['all_min_length']:,} - {s['all_max_length']:,}
 * Date first range  : {s['all_min_create_date']} - {s['all_max_create_date']}
 * Date update range : {s['all_min_update_date']} - {s['all_max_update_date']}
"""

        if step := self.project.read_step_if_in('04_AnalyseChloroplast'):
            # s = step.select([])
            s = step.get_summary_data()
            ng = s['num_genomes']

            def _annotation_summary(title, a, genes):
                t = f"""
{title}:
 * Number of annotations containing IRS         : {s[f'{a}_num_irs']}
 * Number of annotations without IRS            : {ng - s[f'{a}_num_irs']}
 * Number of annotations with wrong orientation : {s[f'{a}_wrong_orientations']}
 * Description wrong orientations               : {s[f'{a}_desc_wrong_orientations']}
 * Number of annotations with offset            : {s[f'{a}_wrong_offset']}
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
"""
            text += _annotation_summary('GeSeq annotation', 'ge_seq', True)
            text += _annotation_summary('NCBI annotation', 'ncbi', True)
            text += _annotation_summary('Sum annotation', 'sum', False)

        # ToDo: Alignments
        # ToDo: Trees

        return dict(text=text)
