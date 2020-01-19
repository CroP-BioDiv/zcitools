from step_project.base_commands import Command


class NCBITaxonomySet(Command):
    _COMMAND = 'ncbi_taxonomy_set'
    _HELP = "Download data and creates SQLite database in global cache"

    def cache_identifier(self):
        return dict(static=True, data_identifier=['NCBI', 'taxonomy'])

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('-f', '--force', action='store_true', help='Force downloading data and databas creation')

    def run(self):
        from .local_database import create_database
        return create_database(self.get_cache_object(), force=self.args.force)
