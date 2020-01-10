from step_project.base_step import Step


class QTLCartStep(Step):
    """
Stores an QTL data and calculation.
Each trait run is stored as a directory.
"""
    _STEP_TYPE = 'qtl_cart'

    # Init object
    def _init_data(self, type_description):
        pass
        # self._sequences = set()  # seq_ident
        # if type_description:
        #     self._sequences.update(type_description['sequences'])

    def _check_data(self):
        pass

    # Set data

    # Save/load data
    def save(self, create=True, completed=True):
        # Store description.yml
        self.save_description(dict(),  # ToDo:
                              create=create, completed=completed)

    # Retrieve data methods

    # Show data
    def show_data(self, params=None):
        print('QTLCartographer', self.directory)
