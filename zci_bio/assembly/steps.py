from step_project.base_step import Step


class AssemblyStep(Step):
    """
Stores assembly method with input data and output assembly sequence.
"""
    _STEP_TYPE = 'assembly'

    # Init object
    def _init_data(self, type_description):
        self._method = None
        if type_description:
            self._method = type_description['method']

    def _check_data(self):
        pass

    # Set data
    def set_data(self, m):
        self._method = m

    # Save/load data
    def save(self, create=True, completed=True):
        # Store description.yml
        self.save_description(dict(method=self._method),
                              create=create, completed=completed)

    # Retrieve data methods

    # Show data
    def show_data(self, params=None):
        print('Assembly', self.directory)
