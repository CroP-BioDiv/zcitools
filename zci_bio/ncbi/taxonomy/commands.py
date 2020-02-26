from step_project.base_commands import NonProjectCommand


class NCBITaxonomySet(NonProjectCommand):  # Works only with cache
    _COMMAND = 'ncbi_taxonomy_set'
    _HELP = "Download data and creates SQLite database in global cache"
    _COMMAND_GROUP = 'NCBI'

    def db_identifier(self):
        return dict(static=True, data_identifier=['NCBI', 'taxonomy'])

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('-f', '--force', action='store_true', help='Force downloading data and database creation')
        parser.add_argument('-d', '--force-db', action='store_true', help='Force databas creation')

    def run(self):
        from .local_database import create_database
        return create_database(self.get_common_db_object(), force=self.args.force, force_db=self.args.force_db)


class NCBITaxonomyNearSpecies(NCBITaxonomySet):  # Works only with common db
    _COMMAND = 'ncbi_taxonomy_tree'
    _HELP = "Show taxonomy tree"

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('rank', help='Start tax rank')
        parser.add_argument('tax_name', help='Start tax name')
        #
        parser.add_argument('-H', '--height', help='Heighest rank to show (species, genus, family)')
        parser.add_argument(
            '-D', '--depth',
            help='Lowest rank to show (subspecies, species, genus, family). Can have more values comma separated.')

    def run(self):
        from .local_database import taxonomy_tree
        return taxonomy_tree(self.get_common_db_object(), self.args)


class NCBITaxonomyWithAssembly(NCBITaxonomyNearSpecies):  # Works only with cache
    _COMMAND = 'ncbi_taxonomy_assembly'
    _HELP = "Show taxonomy tree where species have assembly"

    def run(self):
        from .local_database import taxonomy_tree_assembly
        return taxonomy_tree_assembly(self.get_common_db_object(), self.args)
