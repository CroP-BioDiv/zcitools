from types import SimpleNamespace
from collections import defaultdict
import itertools
from .import_methods import import_ete3_NCBITaxa
from common_utils.cache import cache, cache_args
try:
    from ete3 import Tree
except ImportError:
    Tree = None

# Note: there are more ranks. These ranks have clear hierarchy.
# Ranks filtered from result of SQL query in taxa.sqlite database:
#   SELECT rank, COUNT(*) FROM species GROUP BY rank ORDER BY 1;
NCBI_RANKS = (
    'superkingdom', 'kingdom', 'subkingdom',
    'superphylum', 'phylum', 'subphylum',
    'superclass', 'class', 'subclass', 'infraclass',
    'superorder', 'order', 'suborder', 'infraorder',
    'superfamily', 'family', 'subfamily',
    'tribe', 'subtribe',
    'genus', 'subgenus',
    'species', 'subspecies',
    'varietas', 'subvariety')


def order_ranks(ranks):
    assert all(r in NCBI_RANKS for r in ranks), [r for r in ranks if r not in NCBI_RANKS]
    return sorted(ranks, key=NCBI_RANKS.index)


# NCBI taxonomy is static data (singleton)
_ncbi_taxonomy = None


def get_ncbi_taxonomy():
    global _ncbi_taxonomy
    if _ncbi_taxonomy is None:
        _ncbi_taxonomy = NCBITaxonomy()
    return _ncbi_taxonomy


class NCBITaxonomy:
    # Wrapper around ete3.NCBITaxa
    def __init__(self):
        # Taxid -> SimpleNamespace object with attr:
        #  taxid
        #  parent   : taxid
        #  spname
        #  rank
        #  parents  : list of taxids (linage)
        #  children : list of taxids
        #  down     : set of taxid down from this
        self._cache = dict()

    # Proxy methods
    get_lineage = property(lambda self: self._nt().get_lineage)
    get_rank = property(lambda self: self._nt().get_rank)
    get_taxid_translator = property(lambda self: self._nt().get_taxid_translator)

    @cache
    def _nt(self):
        return import_ete3_NCBITaxa()()  # Make an object

    def _db(self):
        return self._nt().db

    def taxa(self, taxid):
        # Find all needed taxids
        if taxid in self._cache:
            return self._cache[taxid]
        #
        db = self._db()
        all_taxids = []
        to_proc = [taxid]
        while to_proc:
            all_taxids.extend(to_proc)
            result = db.execute(f'SELECT taxid, parent FROM species WHERE parent IN ({",".join(map(str, to_proc))});')
            to_proc = [x[0] for x in result.fetchall() if x[1] not in self._cache]

        # Collect all needed data
        result = db.execute(
            f'SELECT taxid, parent, spname, rank, track FROM species WHERE taxid IN ({",".join(map(str, all_taxids))});')
        result = dict(
            (x[0], SimpleNamespace(
                taxid=x[0], parent=x[1], spname=x[2], rank=x[3], parents=tuple(map(int, x[4].split(',')[1:]))))
            for x in result.fetchall())

        return self._taxa_merge(taxid, result)

    def _taxa_merge(self, taxid, result):
        if taxid in self._cache:
            return self._cache[taxid]

        # Update data
        sn = result[taxid]
        sn.children = [t for t, s in result.items() if s.parent == taxid]
        sn.down = set(sn.children)
        for c in sn.children:
            c_sn = self._taxa_merge(c, result)
            sn.down.update(c_sn.down)
        return sn

    # Names
    def name_2_taxid(self, name):
        taxid = self.get_exact_name_translation(name)
        if taxid is None:
            taxid, _, _ = self._nt().get_fuzzy_name_translation(name)
        return taxid

    def get_fuzzy_name_translation(self, name, sim=0.9):
        name = name[0].upper() + name[1:].lower()  # Better in more cases :-)
        return self._nt().get_fuzzy_name_translation(name, sim=sim)

    def get_exact_name_translation(self, name):
        # Check NCBITaxa.get_fuzzy_name_translation() implementation
        print(f"Trying exact search for {name}")
        # Note: spname column is COLLATE NOCASE!
        result = self._db().execute(f'SELECT taxid FROM species WHERE spname = "{name}" LIMIT 1;')
        taxid = result.fetchone()
        if taxid:
            taxid = int(taxid[0])
            print(f"   FOUND!  taxid {taxid}")
        else:
            print(f"   Didn't find taxid for '{name}' with exact search")
            print(f"Trying synonym search for {name}")
            result = self._db().execute(f'SELECT taxid FROM synonym WHERE spname = "{name}" LIMIT 1;')
            if taxid := result.fetchone():
                taxid = int(taxid[0])
                print(f"   FOUND!  taxid {taxid}")
            else:
                print(f"   Didn't find taxid for '{name}' with synonym search")

        return taxid

    # Group in species
    def group_taxids_in_species(self, taxids):
        nt = self._nt()
        parents = nt.get_lineage_translator(taxids)
        query = ','.join(map(str, itertools.chain(*parents.values())))
        result = self._db().execute(f"SELECT taxid, spname, rank FROM species WHERE taxid IN ({query});")
        id2data = dict((tax_id, (spname, rank)) for tax_id, spname, rank in result.fetchall())
        without_sp = []  # [(taxid, sp_name, rank), ]
        species = defaultdict(list)  # (taxid, sp_name) -> [(taxid, sp_name, rank), ]
        for tax_id, ps in parents.items():
            b_spname, b_rank = id2data[tax_id]
            for p_id in ps[::-1]:
                spname, rank = id2data[p_id]
                if rank == 'species':
                    species[(p_id, spname)].append((tax_id, b_spname, b_rank))
                    break
                if rank == 'genus' or 'family' in rank:
                    without_sp.append((tax_id, b_spname, b_rank))
                    break
        return species, without_sp

    def num_species_below(self, taxid):
        result = self._db().execute(f"SELECT COUNT(*) FROM species WHERE INSTR(track, ',{taxid},') > 0;")
        return int(result.fetchone()[0])

    #
    @cache_args
    def species_tree(self, taxid):
        sn = self.taxa(taxid)
        tree = Tree(name=sn.taxid)  # PhyloTree? Has annotate_ncbi_taxa() method.
        tree.add_feature(spname=sn.spname, rank=sn.rank)
        tree.children.extend(self.species_tree(taxid) for taxid in sn.children)
        return tree

    def find_close_taxids(self, taxid, max_taxid, from_taxids):
        parents = self._nt().get_lineage_translator([taxid] + list(from_taxids))
        for p_id in parents[taxid][-2::-1]:
            try:
                if f_taxids := [t for t in from_taxids if p_id in parents[t]]:
                    return f_taxids
            except KeyError as err:
                print("""
---------------------------------------------------------------------
Parent taxa wasn't found!
Try to upgrade ETE taxonomy database.
In python console run these lines:

from ete3 import NCBITaxa
ncbi = NCBITaxa()
ncbi.update_taxonomy_database()

---------------------------------------------------------------------
""")
                raise err
            if p_id == max_taxid:
                break


class GroupOfSpecies:
    def __init__(self, taxids):  # , ranks=None, names=None):
        # Find all taxids, species and parent clades
        all_taxids = set(taxids)
        ncbi_taxonomy = get_ncbi_taxonomy()
        gl = ncbi_taxonomy.get_lineage
        self.taxid_2_lineage = dict((t, gl(t)) for t in all_taxids)
        all_taxids.update(itertools.chain.from_iterable(self.taxid_2_lineage.values()))
        self.taxid_2_rank = ncbi_taxonomy.get_rank(all_taxids)
        # self.taxid_2_name = ncbi_taxonomy.get_taxid_translator(all_taxids)

    def group_by_rank(self, rank, return_species=False):
        group_taxids = set(t for t, r in self.taxid_2_rank.items() if r == rank)
        groups = defaultdict(list)
        for taxid, linage in self.taxid_2_lineage.items():
            group = list(group_taxids & set(linage))
            assert len(group) == 1, group
            groups[group[0]].append(taxid)
        if not return_species:
            return groups
        #
        species_2_group = dict()
        for group, taxids in groups.items():
            for t in taxids:
                species_2_group[t] = group
        return groups, species_2_group
