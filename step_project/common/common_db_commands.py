import os
from zipfile import ZipFile
from ..base_commands import NonProjectCommand
from common_utils.common_db import CommonDB
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
        parser.add_argument(
            '-n', '--no-recursive', action='store_true', help='No recursive search, search only given path')
        parser.add_argument(
            '-r', '--record', help='Search all occurences of a record(s) in given path. Format <rec>,<rec>,...')
        parser.add_argument(
            '-R', '--records-filename', help='Search all occurences of a records read from given file.')

    def _iterate_records(self, records=None):
        # Returns records as tuples (relative_dir, record ident)

        # Find records to filter
        if records is None:
            records = set()
            if self.args.record:
                records.update(self.args.record.split(','))
            if self.args.records_filename:
                records.update(read_file_as_list(self.args.records_filename))

        #
        if self.args.no_recursive:
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

    def _db_ident(self, rel_dir, record):
        if rel_dir:
            return os.path.join(rel_dir, record)
        return record

    def _set_db_attrs(self):
        idents = tuple(self.args.path.split('/'))
        self._common_db = CommonDB.get_zci_db(idents)
        self._db_dir = self._common_db.db_dir
        if not os.path.isdir(self._db_dir):
            print('No CommondDB directory!', idents, self._db_dir)
            return False
        return True

    def run(self):
        if self._set_db_attrs():
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


class CommonDBDelete(_CommonDBCommand):
    _COMMAND = 'cdb_delete'
    _HELP = "Delete CommonDB record(s)"

    def _run(self):
        for rel_dir, record in self._iterate_records():
            silent_remove_file(self._db_file(rel_dir, record))


class CommonDBCopy(_CommonDBCommand):
    _COMMAND = 'cdb_copy'
    _HELP = "Rename CommonDB record(s) on given path"

    @staticmethod
    def set_arguments(parser):
        parser.add_argument('path', help='CommonDB path')
        parser.add_argument('from_record', help='Record current name')
        parser.add_argument('to_record', help='Record future name')
        parser.add_argument('-p', '--filter-path', help='Filter path')
        parser.add_argument(
            '-n', '--no-recursive', action='store_true', help='No recursive search, search only given path')
        parser.add_argument('-m', '--move', action='store_true', help='Move, delete original')

    def run(self):
        args = self.args
        if args.from_record == args.to_record:
            print('Current and future name are the same!')
            return
        if ',' in args.from_record:
            print('Only one record name can be renamed!!!', args.from_record)
            return
        if not self._set_db_attrs():
            return
        #
        for rel_dir, record in self._iterate_records(records=[args.from_record]):
            current_filename = self._db_file(rel_dir, record)
            #
            with ZipFile(current_filename, 'r') as zip_f:
                # Check is zip of only one file with same name as given record
                files = list(zip_f.infolist())
                if len(files) != 1:
                    print(
                        f"Warning: record {(rel_dir + '/') if rel_dir else ''}{record} has more files!",
                        [z.filename for z in files])
                    continue
                z_i = files[0]
                if z_i.is_dir():
                    print(
                        f"Warning: record's {(rel_dir + '/') if rel_dir else ''}{record} content is only a directory!",
                        z_i.filename)
                    continue
                filename_parts = z_i.filename.split('.', 1)
                if filename_parts[0] != record:
                    print(
                        f"Warning: record's {(rel_dir + '/') if rel_dir else ''}{record} filename has wrong name",
                        z_i.filename, record)
                    continue
                #
                suffix = z_i.filename
                arcname = args.to_record
                if len(filename_parts) == 2:
                    arcname += '.' + filename_parts[1]
                #
                to_record = self._db_ident(rel_dir, args.to_record)
                print(f'Storing data into {to_record}')
                data = zip_f.read(z_i.filename)

                # Replace sequence identifiers
                if filename_parts[1] in ('fa', 'gb'):
                    data = data.decode("utf-8")
                    data = data.replace(args.from_record, args.to_record)
                self._common_db.set_record_from_stream(to_record, data, arcname, force=False, info=True)

            # if args.move:
            #     silent_remove_file(current_filename)
