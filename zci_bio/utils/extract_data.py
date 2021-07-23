import os.path
from functools import wraps
import difflib
from datetime import datetime
from Bio.SeqFeature import FeatureLocation, CompoundLocation
from .entrez import Entrez
from .diff_sequences import Diff
from ..chloroplast.utils import find_chloroplast_irs
from ..chloroplast.irs.small_d import small_d
from ..chloroplast.irs.chloroplot import chloroplot as chloroplot_ann
from ..chloroplast.irs.self_blast import self_blast


def with_seq(func):
    @wraps(func)
    def func_wrapper(self, seq_ident=None, seq=None):
        if seq is None:
            assert seq_ident
            seq = self.sequences_step.get_sequence_record(seq_ident)
        return func(self, seq)
    return func_wrapper


def with_seq_ident(func):
    @wraps(func)
    def func_wrapper(self, seq_ident=None, seq=None):
        if not seq_ident:
            seq_ident = seq.name
        return func(self, seq_ident)
    return func_wrapper


# Note: these methods are not decorators, since wanted functionality depends on two parameters, not on any 'code'
def cache_fetch(key, method):
    def func_wrapper(self, seq_ident=None, seq=None):
        if not seq_ident:
            seq_ident = seq.name
        return self.properties_db.fetch_property(seq_ident, key, method(self), seq_ident=seq_ident, seq=seq)
    return func_wrapper


def cache_fetch_keys1(key, method):
    def func_wrapper(self, seq_idents, seq_step=None):
        if not seq_step:
            seq_step = self.sequences_step
        sr = seq_step.get_sequence_record
        return self.properties_db.fetch_properties_keys1(
            seq_idents, key, lambda s: method(self, seq_ident=s, seq=sr(s)))
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
    def genbank_data(self, seq):
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
    def sra_count(self, seq):
        vals_sra = dict()
        for x in seq.dbxrefs:  # format ['BioProject:PRJNA400982', 'BioSample:SAMN07225454'
            if x.startswith('BioProject:'):
                if bp := x.split(':', 1)[1]:
                    vals_sra['bio_project'] = bp
                    sra_count = Entrez().search_count('sra', term=f"{bp}[BioProject]")
                    vals_sra['sra_count'] = sra_count or None  # None means empty cell :-)
        return vals_sra

    @with_seq
    def annotation(self, seq):
        return self._seq_annotation(seq)

    @with_seq
    def small_d(self, seq):
        return self._small_d_annotation(seq, no_prepend_workaround=True, no_dna_fix=True)

    @with_seq
    def small_d_P(self, seq):
        return self._small_d_annotation(seq, no_prepend_workaround=False, no_dna_fix=True)

    @with_seq
    def small_d_D(self, seq):
        return self._small_d_annotation(seq, no_prepend_workaround=True, no_dna_fix=False)

    @with_seq
    def small_d_all(self, seq):
        return self._small_d_annotation(seq, no_prepend_workaround=False, no_dna_fix=False)

    @with_seq_ident
    def chloroplot(self, seq_ident):
        return chloroplot_ann(self._seq_filename(seq_ident))

    @with_seq_ident
    def pga_sb(self, seq_ident):
        return self._self_blast('pga', seq_ident)

    @with_seq_ident
    def plann_sb(self, seq_ident):
        return self._self_blast('plann', seq_ident)

    # Caching
    # Interface: method_name(seq_ident=None, seq=None)
    cache_genbank_data = cache_fetch('NCBI GenBank data', genbank_data)
    cache_sra_count = cache_fetch('NCBI SRA count', sra_count)
    cache_annotation_ncbi = cache_fetch('annotation ncbi', annotation)
    cache_annotation_ge_seq = cache_fetch('annotation ge_seq', annotation)
    cache_annotation_small_d = cache_fetch('annotation small_d', small_d)
    cache_annotation_small_d_P = cache_fetch('annotation small_d_P', small_d_P)
    cache_annotation_small_d_D = cache_fetch('annotation small_d_D', small_d_D)
    cache_annotation_small_d_all = cache_fetch('annotation small_d_all', small_d_all)
    cache_annotation_chloe = cache_fetch('annotation chloe', annotation)
    cache_annotation_chloroplot = cache_fetch('annotation chloroplot', chloroplot)
    cache_annotation_pga_sb = cache_fetch('annotation pga_sb', pga_sb)
    cache_annotation_plann_sb = cache_fetch('annotation plann_sb', plann_sb)

    # Bulk fetch
    # Interface: method_name(seq_idents, seq_step=None)
    cache_keys1_genbank_data = cache_fetch_keys1('NCBI GenBank data', genbank_data)
    cache_keys1_sra_count = cache_fetch_keys1('NCBI SRA count', sra_count)
    cache_keys1_annotation_ncbi = cache_fetch_keys1('annotation ncbi', annotation)
    cache_keys1_annotation_ge_seq = cache_fetch_keys1('annotation ge_seq', annotation)
    cache_keys1_annotation_small_d = cache_fetch_keys1('annotation small_d', small_d)
    cache_keys1_annotation_small_d_P = cache_fetch_keys1('annotation small_d_P', small_d_P)
    cache_keys1_annotation_small_d_D = cache_fetch_keys1('annotation small_d_D', small_d_D)
    cache_keys1_annotation_small_d_all = cache_fetch_keys1('annotation small_d_all', small_d_all)
    cache_keys1_annotation_chloe = cache_fetch_keys1('annotation chloe', annotation)
    cache_keys1_annotation_chloroplot = cache_fetch_keys1('annotation chloroplot', annotation)
    cache_keys1_annotation_pga_sb = cache_fetch_keys1('annotation pga_sb', pga_sb)
    cache_keys1_annotation_plann_sb = cache_fetch_keys1('annotation plann_sb', plann_sb)

    #
    def _seq_annotation(self, seq):
        if irs := find_chloroplast_irs(seq, check_length=False):
            ira, irb = irs
            # To be sure!
            ira_p = ira.location.parts
            irb_p = irb.location.parts
            d = dict(length=len(seq.seq),
                     ira=[int(ira_p[0].start), int(ira_p[-1].end)],
                     irb=[int(irb_p[0].start), int(irb_p[-1].end)])
            #
            ira_s = ira.extract(seq)
            irb_s = irb.extract(seq)
            if ira.strand == irb.strand:
                irb_s = irb_s.reverse_complement()
            d.update(self._irs_desc(seq, ira_s, irb_s, d))
            return d
        return dict(length=len(seq.seq))

    def _small_d_annotation(self, seq, no_prepend_workaround=True, no_dna_fix=True):
        if res := small_d(seq, no_prepend_workaround=no_prepend_workaround, no_dna_fix=no_dna_fix):
            return self._from_indices(seq, *res)
        return dict(length=len(seq.seq))

    def _self_blast(self, variant, seq_ident):
        seq = self.sequences_step.get_sequence_record(seq_ident)
        if irs := self_blast(variant, self._seq_filename(seq_ident)):
            return self._from_indices(seq, *irs)
        return dict(length=len(seq.seq))

    def _from_indices(self, seq, ira, irb):
        seq_len = len(seq.seq)
        d = dict(length=len(seq.seq),
                 ira=[ira[0], ira[1]],
                 irb=[irb[0], irb[1]])
        #
        ira = self._feature(seq, *ira, 1)
        irb = self._feature(seq, *irb, -1)
        d.update(self._irs_desc(seq, ira.extract(seq), irb.extract(seq), d))
        return d

    def _feature(self, seq, s, e, strand):
        if s < e:
            return FeatureLocation(s, e, strand=strand)
        return CompoundLocation([FeatureLocation(s, len(seq.seq), strand=strand),
                                 FeatureLocation(0, e, strand=strand)])

    def _irs_desc(self, seq, ira, irb, irs_d):
        ira = str(ira.seq)
        irb = str(irb.seq)
        # print(len(ira), len(irb))
        # print('  ira', ira[:20], ira[-20:])
        # print('  irb', irb[:20], irb[-20:])
        if ira == irb:
            return dict(type='+')

        # Check is same IRs region already inspected.
        seq_ident = seq.name.split('.')[0]
        for _, data in self.properties_db.get_properties_key2_like(seq_ident, 'annotation %').items():
            if irs_d['ira'] == data.get('ira') and irs_d['irb'] == data.get('irb'):
                return dict(type=data['type'], diff=data['diff'])

        print(f'diff {seq_ident}: lengths {len(ira)} and {len(irb)}')
        diff = Diff(ira, irb)
        return dict(type=diff.in_short(), diff=diff.get_opcodes())
