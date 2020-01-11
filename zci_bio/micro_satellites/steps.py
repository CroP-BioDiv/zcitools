from step_project.base_step import Step


class QTLCartStep(Step):
    """
Stores an QTL data and calculation.
Each trait run is stored as a directory.
"""
    _STEP_TYPE = 'qtl_cart'

    # Init object
    def _init_data(self, type_description):
        self._num_traits = None
        self._permutations = None
        if type_description:
            self._num_traits = type_description['num_traits']
            self._permutations = type_description['permutations']

    def _check_data(self):
        pass

    # Set data
    def set_data(self, n, p):
        self._num_traits = n
        self._permutations = p

    # Save/load data
    def save(self, create=True, completed=True):
        # Store description.yml
        self.save_description(dict(num_traits=self._num_traits, permutations=self._permutations),
                              create=create, completed=completed)

    # Retrieve data methods
    def trait_dir(self, t_idx):
        return f'T{t_idx:02}'

    def trait_dirs(self):
        return (f'T{t_idx:02}' for t_idx in range(1, self._num_traits + 1))

    # Show data
    def show_data(self, params=None):
        print('QTLCartographer', self.directory)
