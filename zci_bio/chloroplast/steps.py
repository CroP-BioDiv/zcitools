from step_project.base_step import Step


class ChloroplastOrientateStep(Step):
    """
Stores info about chloroplast orientation
"""
    _STEP_TYPE = 'chloroplast_orientate'

    # Init object
    def _init_data(self, type_description):
        pass

    def _check_data(self):
        pass

    # Set data

    # Save/load data
    def save(self, data, create=True, completed=True):
        # Store description.yml
        self.save_description(data, create=create, completed=completed)

    # Retrieve data methods

    # Show data
    def show_data(self, params=None):
        print('Chloroplast orientate', self.directory)
