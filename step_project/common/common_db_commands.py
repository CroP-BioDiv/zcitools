import os
from zipfile import ZipFile
from ..base_commands import NonProjectCommand
from common_utils.common_db import CommonDB, SEQUENCE_DBS_RELATIVE_DIR
from common_utils.terminal_layout import StringColumns
from common_utils.file_utils import read_file_as_list, silent_remove_file


def convert_bytes(num):
    for x in ('b', 'KB', 'MB', 'GB', 'TB'):
        if num < 1024.0:
            return "%3.1f %s" % (num, x)
        num /= 1024.0


class _CommonDBCommand(NonProjectCommand):
    _COMMAND_GROUP = 'CommonDB'

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('path', help='CommonDB path')
        parser.add_argument('-p', '--filter-path', help='Filter path')
        parser.add_argument('-d', '--only-directory', action='store_true', help='Search recursively given path')
        parser.add_argument(
            '-r', '--record', help='Search all occurences of a record(s) in given path. Format <rec>,<rec>,...')
        parser.add_argument(
            '--records-filename', help='Search all occurences of a records read from given file.')
        dbs = CommonDB.get_zci_sequence_dbs()
        parser.add_argument('-S', '--sequence-db', help=f'Sequence database to use: {", ".join(dbs)}')

    def _iterate_records(self):
        # Returns records as tuples (relative_dir, record ident)

        # Find records to filter
        records = set()
        if self.args.record:
            records.update(self.args.record.split(','))
        if self.args.records_filename:
            records.update(read_file_as_list(self.records_filename))

        #
        if self.args.only_directory:
            yield from ((None, f) for f in os.listdir(self._db_dir)
                        if (not records or f in records) and os.path.isfile(os.path.join(self._db_dir, f)))
        else:
            filter_path = self.args.filter_path
            n = len(self._db_dir) + 1
            for root, subdirs, files in os.walk(self._db_dir):
                rel_dir = root[n:]
                if not filter_path or filter_path in rel_dir:
                    yield from ((rel_dir, f) for f in files if not records or f in records)

    def _db_file(self, rel_dir, record):
        if rel_dir:
            return os.path.join(self._db_dir, rel_dir, record)
        return os.path.join(self._db_dir, record)

    def run(self):
        idents = tuple(self.args.path.split('/'))
        db = self.args.sequence_db
        if db:
            CommonDB.get_check_zci_sequence_db(db)
            idents = (SEQUENCE_DBS_RELATIVE_DIR, db) + idents
        self._common_db = CommonDB.get_zci_db(idents)
        self._db_dir = self._common_db.db_dir
        if not os.path.isdir(self._db_dir):
            print('No CommondDB directory!', idents, self._db_dir)
            return
        self._run()


class CommonDBList(_CommonDBCommand):
    _COMMAND = 'cdb_list'
    _HELP = "List CommonDB record's data"

    @staticmethod
    def set_arguments(parser):
        _CommonDBCommand.set_arguments(parser)
        parser.add_argument('-a', '--all-data', action='store_true', help='Print all data inside zip file')
        parser.add_argument('-1', '--only-names', action='store_true', help='Print only names one per a line')

    def _run(self):
        # How to print
        if self.args.only_names:
            header = None

            def _row(record, rel_dir, zip_filename):
                return [record]

        elif self.args.all_data:
            header = ['Record', 'Dir', 'Size', 'Data']

            def _row(record, rel_dir, zip_filename):
                with ZipFile(zip_filename, 'r') as zip_f:
                    zip_files = ', '.join(f"{z_i.filename} ({convert_bytes(z_i.file_size)})"
                                          for z_i in zip_f.infolist() if not z_i.is_dir())
                return [record, rel_dir or '', convert_bytes(os.path.getsize(zip_filename)), zip_files]

        else:
            header = ['Record', 'Dir', 'Size']

            def _row(record, rel_dir, zip_filename):
                return [record, rel_dir or '', convert_bytes(os.path.getsize(zip_filename))]

        #
        rows = [_row(record, rel_dir, self._db_file(rel_dir, record)) for rel_dir, record in self._iterate_records()]

        if rows:
            if self.args.only_names:
                for r in rows:
                    print(r[0])
            else:
                print(StringColumns(rows, header=header))
        else:
            print('No records found!')


class CommonDBGet(_CommonDBCommand):
    _COMMAND = 'cdb_get'
    _HELP = "Extract CommonDB record(s) from given path"

    def _run(self):
        for rel_dir, record in self._iterate_records():
            zip_filename = self._db_file(rel_dir, record)
            # Convert directory structure in filename
            extract_prefix = f"{'_'.join(rel_dir.split(os.sep))}_" if rel_dir else ''
            with ZipFile(zip_filename, 'r') as zip_f:
                for z_i in zip_f.infolist():
                    if not z_i.is_dir():
                        with open(extract_prefix + z_i.filename, 'wb') as save_f:
                            save_f.write(zip_f.read(z_i.filename))


class CommonDBRemove(_CommonDBCommand):
    _COMMAND = 'cdb_remove'
    _HELP = "Remove CommonDB record(s)"

    def _run(self):
        for rel_dir, record in self._iterate_records():
            silent_remove_file(self._db_file(rel_dir, record))
