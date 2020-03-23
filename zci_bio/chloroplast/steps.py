from step_project.base_step import Step


class ChloroplastAnalyseStep(Step):
    """
Stores info about chloroplast genomes.
"""
    _STEP_TYPE = 'chloroplast_analyse'

    # Init object
    def _init_data(self, type_description):
        pass

    def _check_data(self):
        pass

    # Set data

    # Save/load data
    def save(self, create=True, completed=True):
        # Store description.yml
        self.save_description(dict(), create=create, completed=completed)

    # Retrieve data methods

    # Show data
    def show_data(self, params=None):
        print('Chloroplast check project', self.directory)
