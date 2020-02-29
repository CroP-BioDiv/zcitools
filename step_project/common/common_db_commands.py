from ..base_commands import NonProjectCommand
from common_utils.common_db import CommonDB, SEQUENCE_DBS_RELATIVE_DIR


class _CommonDBCommand(NonProjectCommand):
    _COMMAND_GROUP = 'CommonDB'

    @staticmethod
    def set_arguments(parser):
        dbs = CommonDB.get_zci_sequence_dbs()
        parser.add_argument('-S', '--sequence-db', help=f'Sequence database to use: {", ".join(dbs)}')
        parser.add_argument('path', help='CommonDB path')

    def run(self):
        idents = tuple(self.args.path.split('/'))
        db = self.args.sequence_db
        if db:
            CommonDB.get_check_zci_sequence_db(db)
            idents = (SEQUENCE_DBS_RELATIVE_DIR, db) + idents
        self._run(CommonDB.get_zci_db(idents))


class CommonDBList(_CommonDBCommand):
    _COMMAND = 'cdb_list'
    _HELP = "List CommonDB record's data"

    @staticmethod
    def set_arguments(parser):
        _CommonDBCommand.set_arguments(parser)
        parser.add_argument('-a', '--all-data', action='store_true', help='Print all data inside zip file')

    def _run(self, common_db):
        print('ToDo: list')


class CommonDBExtract(_CommonDBCommand):
    _COMMAND = 'cdb_extract'
    _HELP = "Extract CommonDB record(s) from given path"

    @staticmethod
    def set_arguments(parser):
        _CommonDBCommand.set_arguments(parser)
        parser.add_argument('-x', '--recursive', action='store_true', help='Extract recursively given path')
        # parser.add_argument('-r', '--record', help='Extract all occurences of a record in given path')

    def _run(self, common_db):
        print('ToDo: extract')


class CommonDBRemove(_CommonDBCommand):
    _COMMAND = 'cdb_remove'
    _HELP = "Remove CommonDB record(s)"

    @staticmethod
    def set_arguments(parser):
        _CommonDBCommand.set_arguments(parser)
        parser.add_argument('-x', '--recursive', action='store_true', help='Remove recursively given path')
        parser.add_argument('-r', '--record', help='Remove all occurences of a record in given path')

    def _run(self, common_db):
        print('ToDo: remove')
