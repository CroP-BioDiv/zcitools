import os.path
from step_project.base_workflow import BaseWorkflow
from common_utils.exceptions import ZCItoolsValueError
from ..utils.phylogenetic_tree import PhylogeneticTree


class PhylogeneticAnalysis(BaseWorkflow):
    _WORKFLOW = 'phylogenetic_analysis'

    @staticmethod
    def required_parameters():
        return ('input_format', 'input_data', 'data', 'methods')

    @staticmethod
    def format_parameters(params):
        if params['input_format'] in ('list', 'step'):
            params['input_data'] = os.path.abspath(params['input_data'])
        return params

    def _actions(self):
        params = self.parameters
        actions = []

        # Input
        input_format = params['input_format']
        input_data = params['input_data']
        if input_format == 'list':
            if not os.path.isfile(input_data):
                raise ZCItoolsValueError(f"Input data argument ({input_data}) is not a filename!")
            actions.append(
                ('01_accessions_list', f'table -f csv -c seq_ident,seq_ident {input_data}'))
        elif input_format == 'step':
            if not os.path.isfile(input_data):
                raise ZCItoolsValueError(f"Input data argument ({input_data}) is not a directory!")
            assert False, 'Not implemneted! ToDo'
            actions.append(
                ('01_accessions_list', f'copy_step {input_data} -t table'))
        elif input_format == 'accessions':
            assert False, 'Not implemneted! ToDo'
        else:
            raise ZCItoolsValueError(f"Wrong input format {input_format}! Use one of: list, step, accessions")

        # Fetch sequences
        actions.append(('02_seqs', 'fetch_seqs 01_accessions_list'))

        # Annotation
        annotation = params.get('annotate', 'ncbi').lower()
        if annotation == 'ge_seq':
            actions.append(('03_annotations', 'ge_seq 02_seqs'))
        elif annotation == 'ncbi':
            actions.append(('03_annotations', 'seq_subset 02_seqs -t annotations'))
        else:
            raise ZCItoolsValueError(f'Annotation {annotation} is not recognized! Use: ncbi (default), ge_seq.')

        # Data used for phylogeny
        on_data = params['data']
        align = params.get('align', 'mafft')
        parts = ' -w gene' if params.get('use_partitions', 1) else ''
        if on_data == 'whole_original':
            actions.append(('04_alignment', f'align_genomes 03_annotations w -p {align}{parts}'))
        elif on_data == 'whole_standardized':
            assert False, 'ToDo'
        elif on_data == 'all_genes':
            actions.append(('04_alignment', f'align_genomes 03_annotations gc -p {align}{parts}'))
        elif on_data == 'genes':
            assert False, 'ToDo'

        # Phylogenetic analyses
        methods = set(params['methods'].lower().split(','))
        if not_r := [m for m in methods if m not in ('mr_bayes', 'raxml')]:
            raise ZCItoolsValueError(f'Not known phylogenetic method(s): {", ".join(not_r)}!')
        actions.extend((f'05_{m}', f'{m} 04_alignment') for m in methods)

        return actions

    def get_summary(self):
        # ---------------------------------------------------------------------
        # Collect sequence data
        # ---------------------------------------------------------------------
        if not (step := self.project.read_step_if_in('02_seqs')):
            return dict(text='Project not started!')

        outgroup = self.parameters.get('outgroup')

        text = f"""
# Input data

Number of genomes   : {len(step.all_sequences())}
Min sequence length : {min(len(seq) for _, seq in step._iterate_records())}
Max sequence length : {max(len(seq) for _, seq in step._iterate_records())}
Outgroup            : {outgroup or '-'}
"""

        # ---------------------------------------------------------------------
        # Alignments
        # ---------------------------------------------------------------------
        if not (step := self.project.read_step_if_in('04_alignment')):
            return dict(text=text)

        alignment_length = s['alignment_length'] if (s := step.get_summary_data()) else '-'
        text += f"""
# Alignment

Length : {alignment_length}
"""

        # ---------------------------------------------------------------------
        # Trees
        # ---------------------------------------------------------------------
        methods = set(self.parameters['methods'].lower().split(','))
        if len(methods) <= 1:
            text += "\nNo trees to compare!"
        else:
            steps = [(m, self.project.read_step(f'05_{m}')) for m in sorted(methods)]
            assert len(steps) == 2
            assert outgroup
            trees = [PhylogeneticTree(s.get_consensus_file(), outgroup) for _, s in steps]
            rf = trees[0].distance_robinson_foulds(trees[1])

            text += f"""
# Phylogenetic trees

Distance (RF) : {int(rf['rf'])} / {int(rf['max_rf'])}
"""

        return dict(text=text)
