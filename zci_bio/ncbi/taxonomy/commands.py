from step_project.base_commands import NonProjectCommand


class NCBITaxonomySet(NonProjectCommand):  # Works only with cache
    _COMMAND = 'ncbi_taxonomy_set'
    _HELP = "Download data and creates SQLite database in global cache"

    def cache_identifier(self):
        return dict(static=True, data_identifier=['NCBI', 'taxonomy'])

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('-f', '--force', action='store_true', help='Force downloading data and database creation')
        parser.add_argument('-d', '--force-db', action='store_true', help='Force databas creation')

    def run(self):
        from .local_database import create_database
        return create_database(self.get_cache_object(), force=self.args.force, force_db=self.args.force_db)
