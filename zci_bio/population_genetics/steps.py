from step_project.base_step import Step


class NewHybridsStep(Step):
    """
Stores an NewHybrids data and calculation.
Each run is stored as a directory.
"""
    _STEP_TYPE = 'new_hybrids'

    # Init object
    def _init_data(self, type_description):
        self._run_params = None
        if type_description:
            self._run_params = type_description['run_params']

    def _check_data(self):
        pass

    # Set data
    def set_data(self, data_file, gtyp_cat_file, theta_prior, pi_prior, burn_in, num_sweeps):
        self._run_params = dict(
            data_file=data_file, gtyp_cat_file=gtyp_cat_file, theta_prior=theta_prior, pi_prior=pi_prior,
            burn_in=burn_in, num_sweeps=num_sweeps)

    # Save/load data
    def save(self, create=True, completed=True):
        # Store description.yml
        self.save_description(dict(run_params=self._run_params), create=create, completed=completed)

    # Retrieve data methods
    def seed_dir(self, seed):
        return f'seed_{seed[0]}_{seed[1]}'

    # Show data
    def show_data(self, params=None):
        print('NewHybrids', self.directory)
