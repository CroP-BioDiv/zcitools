import os.path
import itertools
import re
from step_project.base_step import Step, StepCollection
from common_utils.file_utils import write_str_in_file, read_file_as_str
from common_utils.exceptions import ZCItoolsValueError
from ..utils.import_methods import import_bio_phylo, import_ete3_Tree, import_dendropy


class RAxMLStep(Step):
    """
Stores an RAxML calculation from alignment.
"""
    # ToDo: other alignment formats?
    _STEP_TYPE = 'raxml'
    _SEQUENCE_TYPES = ('gene', 'genes', 'whole')

    # Init object
    def _init_data(self, type_description):
        self._sequences = set()  # seq_ident
        self._seq_type = None    # None or string from _SEQUENCE_TYPES
        if type_description:
            self._sequences.update(type_description['sequences'])
            self._seq_type = type_description['sequence_type']

    def _check_data(self):
        pass

    # Set data
    def set_sequences(self, seqs):
        self._sequences.update(seqs)

    def seq_sequence_type(self, t):
        if t is not None and t not in self._SEQUENCE_TYPES:
            raise ZCItoolsValueError(f'Wrong alignment sequence type ({t})!')
        self._seq_type = t

    # Save/load data
    def save(self, create=True, completed=True):
        # Store description.yml
        self.save_description(dict(sequences=sorted(self._sequences), sequence_type=self._seq_type),
                              create=create, completed=completed)

    # Retrieve data methods
    # def all_sequences(self):
    #     return self._sequences

    # Show data
    def show_data(self, params=None):
        print('RAxML', self.directory)

    def get_consensus_file(self):
        # Returns filename of consensus tree in newick format
        return self.step_file('RAxML_bestTree.raxml_output')

    def get_consensus_tree(self, library):
        assert library in ('ete', 'dendropy'), library  # 'bio.phylo',
        if library == 'ete':
            # ete3 reads only newick format!
            return import_ete3_Tree()(self.get_consensus_file())
        if library == 'dendropy':
            return import_dendropy().Tree.get(path=self.get_consensus_file(), schema='newick')
        assert False, 'ToDo!!!'


class RAxMLSteps(StepCollection):
    _STEP_TYPE = 'raxmls'
    _SUBSTEP_CLASS = RAxMLStep

    # Show data
    def show_data(self, params=None):
        print('RAxML', self.directory)


#
class MrBayesStep(RAxMLStep):
    "Stores an MrBayes calculation on alignment."
    _STEP_TYPE = 'mr_bayes'

    def get_consensus_file(self):
        cf = self.step_file('consensus.newick')
        if not os.path.isfile(cf):
            # Note: Bio.Phylo can't handle nexus file with more than one comment in a value part!
            # Remove problematic comments from nexus file
            f = self.step_file('consensus.nexus')
            text = read_file_as_str(self.step_file('result.con.tre'))
            write_str_in_file(f, re.sub(r'\[&[^\]]*]', '', text))
            import_bio_phylo().convert(f, 'nexus', cf, 'newick')
        return cf

    def get_consensus_tree(self, library):
        assert library in ('ete', 'dendropy'), library  # 'bio.phylo',
        if library == 'ete':
            # ete3 reads only newick format!
            return import_ete3_Tree()(self.get_consensus_file())
        if library == 'dendropy':
            return import_dendropy().Tree.get(path=self.step_file('result.con.tre'), schema='nexus')
        assert False, 'ToDo!!!'

    #
    def get_p_files(self):
        return [self.step_file(f) for f in self.step_files(matches=r'.*\.run[12]\.p')]

    def get_t_files(self):
        return [self.step_file(f) for f in self.step_files(matches=r'.*\.run[12]\.t')]


class MrBayesSteps(StepCollection):
    "Stores more MrBayes calculations on one or more alignments."
    _STEP_TYPE = 'mr_bayes_s'
    _SUBSTEP_CLASS = MrBayesStep

    def get_p_files(self):
        return list(itertools.chain.from_iterable(s.get_p_files() for s in self.step_objects()))

    def get_t_files(self):
        return list(itertools.chain.from_iterable(s.get_t_files() for s in self.step_objects()))
