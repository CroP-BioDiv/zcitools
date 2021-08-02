import os.path
from functools import wraps
import difflib
from datetime import datetime
from Bio.SeqFeature import FeatureLocation, CompoundLocation
from .entrez import Entrez
from .diff_sequences import diff_check_memory
from ..chloroplast.utils import find_chloroplast_irs


def with_seq(func):
    @wraps(func)
    def func_wrapper(self, seq_ident=None, seq=None, key=None):
        if seq is None:
            assert seq_ident
            seq = self.sequences_step.get_sequence_record(seq_ident, cache=False)
        if not seq_ident:
            seq_ident = seq.name
        return func(self, seq, seq_ident, key)
    return func_wrapper


def with_seq_ident(func):
    @wraps(func)
    def func_wrapper(self, seq_ident=None, seq=None, key=None):
        if not seq_ident:
            seq_ident = seq.name
        return func(self, seq_ident, key)
    return func_wrapper


class ExtractData:
    # Extract/fetch additional data from NCBI GenBank file or by Entrez searches.
    # Idea is to cache that data
    # Two set of methods:
    #  - set of methods that extract or fetch data
    #  - set of methods use upper methods and do caching.
    def __init__(self, properties_db=None, sequences_step=None):
        self.properties_db = properties_db
        self.sequences_step = sequences_step

    def _seq_filename(self, seq_ident):
        return os.path.abspath(self.sequences_step.get_sequence_filename(seq_ident))

    # Extract and fetch methods
    @with_seq
    def genbank_data(self, seq, seq_ident, key):
        annotations = seq.annotations

        vals = dict(length=len(seq.seq))
        if not_dna := [i for i, c in enumerate(str(seq.seq)) if c not in 'ATCG']:
            vals['not_dna'] = not_dna

        vals.update((k, v) for k in ('organism', 'sequence_version') if (v := annotations.get(k)))
        if v := annotations.get('date'):
            vals['update_date'] = datetime.strptime(v, '%d-%b-%Y').date()

        refs = annotations['references']
        if refs[0].title != 'Direct Submission':
            vals['article_title'] = refs[0].title
            vals['journal'] = refs[0].journal
            if refs[0].pubmed_id:
                vals['pubmed_id'] = int(refs[0].pubmed_id)
        if refs[-1].title == 'Direct Submission':
            # ToDo: re ...
            vals['first_date'] = datetime.strptime(
                refs[-1].journal.split('(', 1)[1].split(')', 1)[0], '%d-%b-%Y').date()

        if (sc := annotations.get('structured_comment')) and \
           (ad := sc.get('Assembly-Data')):
            vals['assembly_method'] = ad.get('Assembly Method')
            vals['sequencing_technology'] = ad.get('Sequencing Technology')
        return vals

    @with_seq
    def sra_count(self, seq, seq_ident, key):
        vals_sra = dict()
        for x in seq.dbxrefs:  # format ['BioProject:PRJNA400982', 'BioSample:SAMN07225454'
            if x.startswith('BioProject:'):
                if bp := x.split(':', 1)[1]:
                    vals_sra['bio_project'] = bp
                    sra_count = Entrez().search_count('sra', term=f"{bp}[BioProject]")
                    vals_sra['sra_count'] = sra_count or None  # None means empty cell :-)
        return vals_sra

    @with_seq
    def annotation(self, seq, seq_ident, key):
        if irs := find_chloroplast_irs(seq, check_length=False):
            ira, irb = irs
            d = dict(length=len(seq.seq),
                     ira=self._loc(ira.location.parts),
                     irb=self._loc(irb.location.parts))
            #
            ira_s = ira.extract(seq)
            irb_s = irb.extract(seq)
            if ira.strand == irb.strand:
                irb_s = irb_s.reverse_complement()
            if desc := self._irs_desc(seq, seq_ident, key, ira_s, irb_s, d):
                d.update(desc)
                return d
            return None
        return dict(length=len(seq.seq))

    @staticmethod
    def _loc(ir_parts):
        assert 1 <= len(ir_parts) <= 2, ir_parts
        p1 = ir_parts[0]
        if len(ir_parts) == 1:
            return [int(p1.start), int(p1.end)]
        if p1.start == 0:
            return [int(ir_parts[1].start), int(p1.end)]
        assert ir_parts[1].start == 0
        return [int(p1.start), int(ir_parts[1].end)]

    @with_seq
    def small_d(self, seq, seq_ident, key):
        return self._small_d_annotation(seq, seq_ident, key, no_prepend_workaround=True, no_dna_fix=True)

    @with_seq
    def small_d_P(self, seq, seq_ident, key):
        return self._small_d_annotation(seq, seq_ident, key, no_prepend_workaround=False, no_dna_fix=True)

    @with_seq
    def small_d_D(self, seq, seq_ident, key):
        return self._small_d_annotation(seq, seq_ident, key, no_prepend_workaround=True, no_dna_fix=False)

    @with_seq
    def small_d_all(self, seq, seq_ident, key):
        return self._small_d_annotation(seq, seq_ident, key, no_prepend_workaround=False, no_dna_fix=False)

    @with_seq_ident
    def chloroplot(self, seq_ident, key):
        from ..chloroplast.irs.chloroplot import chloroplot as chloroplot_ann
        return self._from_indices(
            self.sequences_step.get_sequence_record(seq_ident, cache=False), seq_ident, key,
            chloroplot_ann(self._seq_filename(seq_ident)))

    @with_seq_ident
    def pga(self, seq_ident, key):
        from ..chloroplast.irs.pga import pga
        return self._from_indices(
            self.sequences_step.get_sequence_record(seq_ident, cache=False), seq_ident, key,
            pga(self._seq_filename(seq_ident)))

    @with_seq_ident
    def pga_sb(self, seq_ident, key):
        return self._self_blast('pga', seq_ident, key)

    @with_seq_ident
    def plann(self, seq_ident, key):
        from ..chloroplast.irs.plann import plann
        return self._from_indices(
            self.sequences_step.get_sequence_record(seq_ident, cache=False), seq_ident, key,
            plann(self._seq_filename(seq_ident)))

    @with_seq_ident
    def plann_sb(self, seq_ident, key):
        return self._self_blast('plann', seq_ident, key)

    @with_seq_ident
    def org_annotate(self, seq_ident, key):
        from ..chloroplast.irs.org_annotate import org_annotate
        return self._from_indices(
            self.sequences_step.get_sequence_record(seq_ident, cache=False), seq_ident, key,
            org_annotate(self._seq_filename(seq_ident)))

    #
    def _small_d_annotation(self, seq, seq_ident, key, no_prepend_workaround=True, no_dna_fix=True):
        from ..chloroplast.irs.small_d import small_d
        return self._from_indices(
            seq, seq_ident, key,
            small_d(seq, no_prepend_workaround=no_prepend_workaround, no_dna_fix=no_dna_fix))

    def _self_blast(self, variant, seq_ident, key):
        from ..chloroplast.irs.self_blast import self_blast
        return self._from_indices(
            self.sequences_step.get_sequence_record(seq_ident, cache=False), seq_ident, key,
            self_blast(variant, self._seq_filename(seq_ident)))

    #
    def _from_indices(self, seq, seq_ident, key, irs):
        if irs:
            ira, irb = irs
            seq_len = len(seq.seq)
            d = dict(length=len(seq.seq),
                     ira=[ira[0], ira[1]],
                     irb=[irb[0], irb[1]])
            #
            ira = self._feature(seq, *ira, 1)
            irb = self._feature(seq, *irb, -1)
            if desc := self._irs_desc(seq, seq_ident, key, ira.extract(seq), irb.extract(seq), d):
                d.update(desc)
                return d
            return None
        print(f'{key.split(" ", 1)[1]} {seq_ident}: no IRs')
        return dict(length=len(seq.seq))

    def _feature(self, seq, s, e, strand):
        if s < e:
            return FeatureLocation(s, e, strand=strand)
        return CompoundLocation([FeatureLocation(s, len(seq.seq), strand=strand),
                                 FeatureLocation(0, e, strand=strand)])

    def _irs_desc(self, seq, seq_ident, key, ira, irb, irs_d):
        ira = str(ira.seq)
        irb = str(irb.seq)
        print(f'{key.split(" ", 1)[1]} {seq_ident}: {len(ira)}, {len(irb)}: {irs_d["ira"]}, {irs_d["irb"]}')
        # if len(ira) > 30000 or len(irb) > 30000:
        #     return  # Hack: for now
        if ira == irb:
            return dict(type='+')

        # Check is same IRs region already inspected.
        seq_ident = seq.name.split('.')[0]
        for _, data in self.properties_db.get_properties_key2_like(seq_ident, 'annotation %').items():
            if irs_d['ira'] == data.get('ira') and irs_d['irb'] == data.get('irb'):
                print(f'  found diff: {data["type"]}')
                return dict(type=data['type'], diff=data['diff']) if 'diff' in data else dict(type=data['type'])

        # Problem:
        # Original NCBI annotation of NC_031898, has IR lengths 23595 and 74817.
        # IRn is annotated as complement(10932..85748). It should be 85751..109342.
        # Copy/paste error?
        if len(irb) > 3 * len(ira) / 2 or len(ira) > 3 * len(irb) / 2:
            return dict(type='???')
        print(f'  calculate diff')
        diff = diff_check_memory(ira, irb)
        return dict(type=diff.in_short(), diff=diff.get_opcodes())


# Add cache methods into ExtractData class
# Note: these methods are not decorators, but quite similar
def cache_fetch(key, method):
    def func_wrapper(self, seq_ident=None, seq=None):
        if not seq_ident:
            seq_ident = seq.name
        return self.properties_db.fetch_property(seq_ident, key, method(self), seq_ident=seq_ident, seq=seq, key=key)
    return func_wrapper


def cache_fetch_keys1(key, method):
    def func_wrapper(self, seq_idents, seq_step=None):
        if not seq_step:
            seq_step = self.sequences_step
        sr = seq_step.get_sequence_record
        return self.properties_db.fetch_properties_keys1(
            seq_idents, key, lambda s: method(self, seq_ident=s, seq=sr(s, cache=False), key=key))
    return func_wrapper


for key, m_sufix, cls_method in (
        ('NCBI GenBank data', 'genbank_data', ExtractData.genbank_data),
        ('NCBI SRA count', 'sra_count', ExtractData.sra_count),
        ('annotation ncbi', 'annotation_ncbi', ExtractData.annotation),
        ('annotation ge_seq', 'annotation_ge_seq', ExtractData.annotation),
        ('annotation small_d', 'annotation_small_d', ExtractData.small_d),
        ('annotation small_d_P', 'annotation_small_d_P', ExtractData.small_d_P),
        ('annotation small_d_D', 'annotation_small_d_D', ExtractData.small_d_D),
        ('annotation small_d_all', 'annotation_small_d_all', ExtractData.small_d_all),
        ('annotation chloe', 'annotation_chloe', ExtractData.annotation),
        ('annotation chloroplot', 'annotation_chloroplot', ExtractData.chloroplot),
        ('annotation pga', 'annotation_pga', ExtractData.pga),
        ('annotation pga_sb', 'annotation_pga_sb', ExtractData.pga_sb),
        ('annotation plann', 'annotation_plann', ExtractData.plann),
        ('annotation plann_sb', 'annotation_plann_sb', ExtractData.plann_sb),
        ('annotation org_annotate', 'annotation_org_annotate', ExtractData.org_annotate),
        ):
    # Caching
    # Interface: method_name(seq_ident=None, seq=None)
    setattr(ExtractData, f'cache_{m_sufix}', cache_fetch(key, cls_method))
    # Bulk fetch
    # Interface: method_name(seq_idents, seq_step=None)
    setattr(ExtractData, f'cache_keys1_{m_sufix}', cache_fetch_keys1(key, cls_method))
