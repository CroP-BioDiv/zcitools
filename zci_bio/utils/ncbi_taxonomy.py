from types import SimpleNamespace
from collections import defaultdict
import itertools
from .import_methods import import_ete3_NCBITaxa
from common_utils.cache import cache, cache_args
try:
    from ete3 import Tree
except ImportError:
    Tree = None

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
        result = self._db().execute(f'SELECT taxid FROM species WHERE lower(spname) = "{name.lower()}" LIMIT 1;')
        taxid = result.fetchone()
        if taxid:
            taxid = int(taxid[0])
            print(f"   FOUND!  taxid {taxid}")
        else:
            print(f"   Didn't find taxid for '{name}'")

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
