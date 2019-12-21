import os.path
from zcitools.base.step import Step, StepCollection


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


class RAxMLSteps(StepCollection):
    _STEP_TYPE = 'raxmls'
    _SUBSTEP_CLASS = RAxMLStep

    # Show data
    def show_data(self, params=None):
        print('MrBayes', self.directory)


#
class MrBayesStep(RAxMLStep):
    """
Stores an MrBayes calculation from alignment.
"""
    # ToDo: other alignment formats?
    _STEP_TYPE = 'mr_bayes'


class MrBayesSteps(StepCollection):
    _STEP_TYPE = 'mr_bayes_s'
    _SUBSTEP_CLASS = MrBayesStep
