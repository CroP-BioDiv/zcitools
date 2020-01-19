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


def create_database(cache_obj, force=False):
    if not cache_obj:
        raise ZCItoolsValueError('No cache specified!')

    cache_obj.ensure_location()
    zip_filename = cache_obj.get_record_filename(_NCBI_ZIP_FILENAME)
    db_filename = cache_obj.get_record_filename(_DB_FILENAME)

    if os.path.isfile(db_filename) and not force:
        print('Taxonomy database already created!')
        return

    if not os.path.isfile(zip_filename) or force:
        download_url('ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip', zip_filename)

    # Create database
    silent_remove_file(db_filename)  # Remove file to be sure
    conn = sqlite3.connect(db_filename)
    cursor = conn.cursor()

    # Note: not all tables are used
    with ZipFile(zip_filename, 'r') as zip_f:
        _fill_table(
            conn, cursor, zip_f, 'nodes',
            ('tax_id INTEGER NOT NULL PRIMARY KEY',
             'parent_tax_id INTEGER',
             'rank TEXT', 'embl_code TEXT', 'division_id TEXT', 'inherited_div INTEGER',
             'genetic_code INTEGER', 'inherited_GC INTEGER',
             'mitochondrial_genetic_code INTEGER', 'inherited_MGC INTEGER',
             'GenBank_hidden INTEGER', 'hidden_subtree_root INTEGER',
             'comments TEXT', 'plastid_genetic_code INTEGER', 'inherited_PGC INTEGER',
             'specified_species INTEGER', 'hydrogenosome_genetic_code INTEGER', 'inherited_HGC INTEGER'))

        _fill_table(
            conn, cursor, zip_f, 'names',
            ('tax_id INTEGER NOT NULL', 'name TEXT', 'unique_name TEXT', 'name_class TEXT'))

        _fill_table(
            conn, cursor, zip_f, 'gencode',
            ('code_id INTEGER NOT NULL PRIMARY KEY',
             'abbreviation TEXT', 'name TEXT', 'cde TEXT', 'starts TEXT'))

        _fill_table(
            conn, cursor, zip_f, 'rankedlineage',
            ('tax_id INTEGER NOT NULL PRIMARY KEY',
             'tax_name TEXT', 'species TEXT', 'genus TEXT', 'family TEXT', 'order_ TEXT',
             'class TEXT', 'phylum TEXT', 'kingdom TEXT', 'superkingdom TEXT'))

        _fill_table(
            conn, cursor, zip_f, 'taxidlineage',
            ('tax_id INTEGER NOT NULL PRIMARY KEY', 'lineage TEXT'))

        # ToDo: create indices

    conn.commit()
    conn.close()


def _fill_table(conn, cursor, zip_f, table_name, columns):
    cursor.execute(f"CREATE TABLE {table_name} ({', '.join(columns)})")

    rows = zip_f.read(f'{table_name}.dmp').decode('utf-8').splitlines()
    rows = [r.strip().split('\t|\t') for r in rows]
    cursor.executemany(
        f'INSERT INTO {table_name} VALUES ({",".join(["?"] * len(columns))})', rows)

    conn.commit()
    print(table_name, cursor.execute(f'SELECT COUNT(*) FROM {table_name}').fetchone()[0])
