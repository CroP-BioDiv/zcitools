import os.path
import sqlite3
import datetime
from zipfile import ZipFile, ZIP_BZIP2
from common_utils.net_utils import download_url
from common_utils.file_utils import silent_remove_file
from common_utils.exceptions import ZCItoolsValueError
from common_utils.terminal_layout import TreeBox
from common_utils.misc import split_list
from ...utils.entrez import Entrez

# Database is created with new_taxdump data. Check readme file:
# ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/taxdump_readme.txt
# Data to use are in ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip
# Data and database file are located in cache directory with names taxonomy.zip and taxonomy.db.
_NCBI_ZIP_FILENAME = 'taxonomy.zip'
_DB_FILENAME = 'taxonomy.db'
_NUM_QUERY_SPECIES = 500

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
_short_ranks = dict(
    family='f', forma='forma', genus='gen', order='o', section='section', series='series',
    species='sp', subclass='scl', subfamily='sfam', subgenus='sg', suborder='so', subsection='ssec',
    subspecies='ssp', subtribe='st', tribe='t', varietas='var')
_short_ranks['no rank'] = 'no'
_short_ranks['species group'] = 'sp g'
_ranks_have_ncbi_data = ('species', 'subspecies', 'varietas', 'no rank')

_taxa_column_names = ('tax_id', 'parent_tax_id', 'rank', 'tax_name', 'species', 'genus', 'family')


def _tax_to_dict(row):
    return dict(zip(_taxa_column_names, row))


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
        # s_idx = sql_cs['superkingdom']  # ('Archaea', 'Bacteria', 'Eukaryota'))
        # k_idx = sql_cs['kingdom']       # 'Viridiplantae'
        # p_idx = sql_cs['phylum']        # 'Streptophyta'
        c_idx = sql_cs['class']           # 'Magnoliopsida'
        rankedlineage = _read_rows_dict(zip_f, 'rankedlineage', lambda row: row[c_idx] == 'Magnoliopsida')
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
        # _add_columns(columns, rows, taxidlineage, _taxidlineage_columns, 'lineage')

        _fill_whole_table(conn, cursor, 'taxa', columns, rows)

        cursor.execute(f"CREATE INDEX taxa_parent_tax_id ON taxa ('parent_tax_id')")
        cursor.execute(f"CREATE INDEX taxa_species ON taxa ('species')")
        cursor.execute(f"CREATE INDEX taxa_genus ON taxa ('genus')")
        cursor.execute(f"CREATE INDEX taxa_family ON taxa ('family')")

        # NCBI calls cache data
        cursor.execute(
            "CREATE TABLE cached_data (tax_id INTEGER, tag TEXT, date TEXT, value TEXT, PRIMARY KEY (tax_id, tag))")

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
    _fill_whole_table(conn, cursor, table_name, columns, list(_read_rows(zip_f, table_name)))


def _fill_whole_table(conn, cursor, table_name, columns, rows):
    cursor.execute(f"CREATE TABLE {table_name} ({', '.join(columns)})")

    cursor.executemany(
        f'INSERT INTO {table_name} VALUES ({",".join(["?"] * len(columns))})', rows)

    conn.commit()
    print(f'  Info: created table {table_name}',
          cursor.execute(f'SELECT COUNT(*) FROM {table_name}').fetchone()[0])


def _read_rows(zip_f, table_name):
    rows = zip_f.read(f'{table_name}.dmp').decode('utf-8').splitlines()
    return ([x.strip() for x in r[:-2].split('\t|\t')] for r in rows)


def _read_rows_dict(zip_f, table_name, filter_method):
    return dict((r[0], r) for r in _read_rows(zip_f, table_name) if filter_method(r))


# ---------------------------------------------------------
# Query taxa
# ---------------------------------------------------------
class _TaxonomyNode:
    def __init__(self, cursor, row, rank_low):
        self.row = row  # dict
        self._data = None
        #
        if rank_low and row['rank'] in rank_low:
            self.children = []
        else:
            cursor.execute(f'SELECT * FROM taxa WHERE parent_tax_id = {row["tax_id"]}')
            # Note: fetchall() is important, since same cursor is used in this recursion!
            self.children = [_TaxonomyNode(cursor, _tax_to_dict(c_row), rank_low) for c_row in cursor.fetchall()]

    @property
    def label(self):
        if self._data:
            return f"({_short_ranks[self.row['rank']]}) {self.row['tax_name']} [{self._data}]"
        return f"({_short_ranks[self.row['rank']]}) {self.row['tax_name']}"

    def _get_nodes_of_rank(self, ranks, nodes):
        if self.row['rank'] in ranks:
            nodes[self.row['tax_id']] = self
        for c in self.children:
            c._get_nodes_of_rank(ranks, nodes)

    def _filter_with_data(self):
        self.children = [c for c in self.children if c._filter_with_data()]
        return self._data or bool(self.children)


class _TaxonomyTree:
    # Represent taxonomy tree.
    # Root is given with given name and rank, or it's ancestor.
    # Optional depth is given with rank_low.
    def __init__(self, db_filename, rank, tax_name, rank_heigh=None, rank_low=None):
        if not os.path.isfile(db_filename):
            raise ZCItoolsValueError(f"Taxonomy database {db_filename} doesn't exist!")

        self.conn = sqlite3.connect(db_filename)
        self.cursor = self.conn.cursor()

        # Find start row
        self.start_row = self._start_row(rank, tax_name)
        if self.start_row and rank_heigh:
            self.start_row = self._ancestor_start_row(self.start_row, rank_heigh)

        # Find subtree
        self.root = None
        if self.start_row:
            self.root = _TaxonomyNode(self.cursor, self.start_row, rank_low)

    def close(self):
        self.conn.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    #
    def _start_row(self, rank, tax_name):
        self.cursor.execute('SELECT * FROM taxa WHERE rank = ? and tax_name LIKE ? LIMIT 2', (rank, tax_name + '%'))
        rows = self.cursor.fetchall()
        if not rows:
            print('No rows fetched with given filter!')
            return
        if len(rows) == 2:
            print('More rows fetched with given filter!')
            return
        return _tax_to_dict(rows[0])

    def _ancestor_start_row(self, row, rank_heigh):
        start_row = row
        while True:
            ac_row = self.cursor.execute(f'SELECT * FROM taxa WHERE tax_id = {row["parent_tax_id"]}').fetchone()
            if not ac_row:
                return start_row
            ac_row = _tax_to_dict(ac_row)
            if ac_row['rank'] in rank_heigh:
                return ac_row
            #
            row = ac_row

    #
    def set_data_to_nodes_of_rank(self, tag, ranks, fetch_method):
        # Find nodes to check
        nodes = dict()
        self.root._get_nodes_of_rank(ranks, nodes)

        # What is cached
        self.cursor.execute(
            f"SELECT * FROM cached_data WHERE tag = '{tag}' and tax_id in ({','.join(map(str, nodes.keys()))})")
        for row in self.cursor:
            nodes[row[0]]._data = row[-1]  # Can be None
            nodes.pop(row[0])

        if nodes:
            f_data = fetch_method(nodes)
            for tax_id, value in f_data.items():
                nodes[tax_id]._data = value
            # Store into database
            d = str(datetime.date.today())  # ISO format
            self.cursor.executemany(f'INSERT INTO cached_data VALUES (?, ?, ?, ?)',
                                    [(tax_id, tag, d, f_data.get(tax_id)) for tax_id, n in nodes.items()])
            self.conn.commit()
        #
        self.root._filter_with_data()

    #
    def print_tree(self):
        if not self.root:
            print('No tree!')
            return
        print(TreeBox(self.root))


#
def _get_tree(cache_obj, params):
    db_filename = cache_obj.get_record_filename(_DB_FILENAME)

    if not os.path.isfile(db_filename):
        print("Error: Taxonomy database doesn't exist!")
        return

    return _TaxonomyTree(
        db_filename, params.rank, params.tax_name,
        rank_heigh=params.height.split(',') if params.height else None,
        rank_low=params.depth.split(',') if params.depth else None)


def taxonomy_tree(cache_obj, params):
    with _get_tree(cache_obj, params) as tt:
        tt.print_tree()


def taxonomy_tree_assembly(cache_obj, params):
    with _get_tree(cache_obj, params) as tt:
        tt.set_data_to_nodes_of_rank('assembly', _ranks_have_ncbi_data, _find_assemblies)
        tt.print_tree()


def _find_assemblies(sp_nodes):
    entrez = Entrez()

    # Fetch assembly ids
    print(f'Info: seaching for {len(sp_nodes)} species!', flush=True)
    assemblies = []
    for ors in split_list([n.row['tax_name'] for n in sp_nodes.values()], _NUM_QUERY_SPECIES):
        print('.', end='', flush=True)
        data = entrez.esearch(db='assembly', term=' OR '.join(f'"{o}"[ORGN]' for o in ors), retmax=_NUM_QUERY_SPECIES)
        assemblies.extend(data['IdList'])
    print()

    # Fetch assemblies to see what species are in
    print(f'Info: seaching for {len(assemblies)} assemblies!', flush=True)
    ret = dict()
    for ass_ids in split_list(list(assemblies), _NUM_QUERY_SPECIES):
        print('.', end='', flush=True)
        records = entrez.esummary(db='assembly', id=','.join(ass_ids), retmax=_NUM_QUERY_SPECIES)
        # (?) a['SpeciesTaxid']
        # for a in records['DocumentSummarySet']['DocumentSummary']:
        #     print(a['SpeciesTaxid'], a['Taxid'], a['AssemblyAccession'])
        ret.update((int(a['Taxid']), a['AssemblyAccession']) for a in records['DocumentSummarySet']['DocumentSummary']
                   if int(a['Taxid']) in sp_nodes)
    print()
    return ret
