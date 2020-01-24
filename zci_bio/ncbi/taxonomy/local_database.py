import os.path
import sqlite3
from zipfile import ZipFile, ZIP_BZIP2
from common_utils.net_utils import download_url
from common_utils.file_utils import silent_remove_file
from common_utils.exceptions import ZCItoolsValueError

# Database is created with new_taxdump data. Check readme file:
# ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/taxdump_readme.txt
# Data to use are in ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip
# Data and database file are located in cache directory with names taxonomy.zip and taxonomy.db.
_NCBI_ZIP_FILENAME = 'taxonomy.zip'
_DB_FILENAME = 'taxonomy.db'

_gencode_columns = (
    'code_id INTEGER NOT NULL PRIMARY KEY', 'abbreviation TEXT', 'name TEXT', 'cde TEXT', 'starts TEXT')

_nodes_columns = (
    'tax_id INTEGER NOT NULL PRIMARY KEY',
    'parent_tax_id INTEGER',
    'rank TEXT', 'embl_code TEXT', 'division_id TEXT', 'inherited_div INTEGER',
    'genetic_code INTEGER', 'inherited_GC INTEGER',
    'mitochondrial_genetic_code INTEGER', 'inherited_MGC INTEGER',
    'GenBank_hidden INTEGER', 'hidden_subtree_root INTEGER',
    'comments TEXT', 'plastid_genetic_code INTEGER', 'inherited_PGC INTEGER',
    'specified_species INTEGER', 'hydrogenosome_genetic_code INTEGER', 'inherited_HGC INTEGER')
_names_columns = ('tax_id INTEGER NOT NULL', 'name TEXT', 'unique_name TEXT', 'name_class TEXT')
_rankedlineage_columns = (
    'tax_id INTEGER NOT NULL PRIMARY KEY',
    'tax_name TEXT', 'species TEXT', 'genus TEXT', 'family TEXT', 'order_ TEXT',
    'class TEXT', 'phylum TEXT', 'kingdom TEXT', 'superkingdom TEXT')
_taxidlineage_columns = ('tax_id INTEGER NOT NULL PRIMARY KEY', 'lineage TEXT')

_delnodes_columns = ('tax_id INTEGER NOT NULL PRIMARY KEY',)
_merged_columns = ('old_tax_id INTEGER NOT NULL PRIMARY KEY', 'new_tax_id INTEGER')


def create_database(cache_obj, force=False, force_db=False):
    if not cache_obj:
        raise ZCItoolsValueError('No cache specified!')

    cache_obj.ensure_location()
    zip_filename = cache_obj.get_record_filename(_NCBI_ZIP_FILENAME)
    db_filename = cache_obj.get_record_filename(_DB_FILENAME)

    if force:
        silent_remove_file(zip_filename)
        silent_remove_file(db_filename)
    elif force_db:
        silent_remove_file(db_filename)

    if os.path.isfile(db_filename):
        print('Taxonomy database already created!')
        return

    if not os.path.isfile(zip_filename):
        download_url('ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip', zip_filename)

    # Create database
    conn = sqlite3.connect(db_filename)
    cursor = conn.cursor()

    # Note: not all tables are used
    with ZipFile(zip_filename, 'r') as zip_f:
        _fill_table(conn, cursor, zip_f, 'gencode', _gencode_columns)

        # These table are merged into table taxa
        # _fill_table(conn, cursor, zip_f, 'nodes', _nodes_columns)
        # _fill_table(conn, cursor, zip_f, 'names', _names_columns)
        # _fill_table(conn, cursor, zip_f, 'rankedlineage', _rankedlineage_columns)
        # _fill_table(conn, cursor, zip_f, 'taxidlineage', _taxidlineage_columns)

        sql_cs = dict((c.split(' ', 1)[0], i) for i, c in enumerate(_rankedlineage_columns))
        # not_empty_idx = [sql_cs[c] for c in ('species', 'genus', 'family')]
        # s_idx = sql_cs['superkingdom']  # not in ('Archaea', 'Bacteria', 'Eukaryota'))
        k_idx = sql_cs['kingdom']  # 'Viridiplantae'
        rankedlineage = _read_rows_dict(zip_f, 'rankedlineage', lambda row: row[k_idx] == 'Viridiplantae')
            # lambda row: any(row[i] for i in not_empty_idx) and row[s_idx] == 'Eukaryota')
        print('  Info: read rankedlineage', len(rankedlineage))

        nodes = _read_rows_dict(zip_f, 'nodes', lambda row: row[0] in rankedlineage)
        print('  Info: read nodes', len(nodes))
        # names = _read_rows_dict(zip_f, 'names')
        taxidlineage = _read_rows_dict(zip_f, 'taxidlineage', lambda row: row[0] in rankedlineage)
        print('  Info: read taxidlineage', len(rankedlineage))

        rows = [[tax_id] for tax_id in rankedlineage.keys()]
        columns = ['tax_id INTEGER NOT NULL PRIMARY KEY']
        _add_columns(columns, rows, nodes, _nodes_columns, 'parent_tax_id', 'rank')
        _add_columns(
            columns, rows, rankedlineage, _rankedlineage_columns,
            'tax_name', 'species', 'genus', 'family')  # , 'order_', 'class', 'phylum', 'kingdom', 'superkingdom')
        _add_columns(columns, rows, taxidlineage, _taxidlineage_columns, 'lineage')

        _fill_whole_table(conn, cursor, 'taxa', columns, rows)

        cursor.execute(f"CREATE INDEX taxa_species ON taxa ('species')")
        cursor.execute(f"CREATE INDEX taxa_genus ON taxa ('genus')")
        cursor.execute(f"CREATE INDEX taxa_family ON taxa ('family')")
        # ToDo: create indices

    conn.commit()
    conn.close()


def _add_columns(columns, rows, table_data, sql_columns, *add_columns):
    # Append column definition and row data from ncbi table data into given lists.
    sql_cs = dict((c.split(' ', 1)[0], (c, i)) for i, c in enumerate(sql_columns))
    extract_row_ids = []
    for c in add_columns:
        dc, ri = sql_cs[c]
        columns.append(dc)
        extract_row_ids.append(ri)

    for row in rows:
        t_row = table_data[row[0]]
        row.extend(t_row[i] for i in extract_row_ids)


def _fill_table(conn, cursor, zip_f, table_name, columns):
    _fill_whole_table(conn, cursor, table_name, columns, _read_rows(zip_f, table_name))


def _fill_whole_table(conn, cursor, table_name, columns, rows):
    cursor.execute(f"CREATE TABLE {table_name} ({', '.join(columns)})")

    cursor.executemany(
        f'INSERT INTO {table_name} VALUES ({",".join(["?"] * len(columns))})', rows)

    conn.commit()
    print(f'  Info: created table {table_name}',
          cursor.execute(f'SELECT COUNT(*) FROM {table_name}').fetchone()[0])


def _read_rows(zip_f, table_name):
    rows = zip_f.read(f'{table_name}.dmp').decode('utf-8').splitlines()
    return [[x.strip() for x in r[:-2].split('\t|\t')] for r in rows]


def _read_rows_dict(zip_f, table_name, filter_method):
    return dict((r[0], r) for r in _read_rows(zip_f, table_name) if filter_method(r))
