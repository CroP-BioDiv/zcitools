from functools import wraps
import difflib
from datetime import datetime
from Bio.SeqFeature import FeatureLocation, CompoundLocation
from .entrez import Entrez
from ..chloroplast.utils import find_chloroplast_irs
from ..chloroplast.irs.small_d import small_d


def with_seq(func):
    @wraps(func)
    def func_wrapper(self, seq_ident=None, seq=None):
        if seq is None:
            assert seq_ident
            seq = self.sequences_step.get_sequence_record(seq_ident)
        # if not seq_ident:
        #     seq_ident = seq.name
        return func(self, seq_ident=seq_ident, seq=seq)
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

    # Extract and fetch methods
    @with_seq
    def genbank_data(self, seq_ident=None, seq=None):
        annotations = seq.annotations

        vals = dict(length=len(seq.seq))

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
    def sra_count(self, seq_ident=None, seq=None):
        vals_sra = dict()
        for x in seq.dbxrefs:  # format ['BioProject:PRJNA400982', 'BioSample:SAMN07225454'
            if x.startswith('BioProject:'):
                if bp := x.split(':', 1)[1]:
                    vals_sra['bio_project'] = bp
                    sra_count = Entrez().search_count('sra', term=f"{bp}[BioProject]")
                    vals_sra['sra_count'] = sra_count or None  # None means empty cell :-)
        return vals_sra

    @with_seq
    def annotation(self, seq_ident=None, seq=None):
        return _seq_annotation(seq)

    @with_seq
    def small_d(self, seq_ident=None, seq=None):
        return _small_d_annotation(seq)

    # Caching
    # Interface: method_name(seq_ident=None, seq=None)
    cache_genbank_data = cache_fetch('NCBI GenBank data', genbank_data)
    cache_sra_count = cache_fetch('NCBI SRA count', sra_count)
    cache_annotation_ncbi = cache_fetch('annotation ncbi', annotation)
    cache_annotation_ge_seq = cache_fetch('annotation ge_seq', annotation)
    cache_annotation_small_d = cache_fetch('annotation small_d', small_d)

    # Bulk fetch
    # Interface: method_name(seq_idents, seq_step=None)
    cache_keys1_genbank_data = cache_fetch_keys1('NCBI GenBank data', genbank_data)
    cache_keys1_sra_count = cache_fetch_keys1('NCBI SRA count', sra_count)
    cache_keys1_annotation_ncbi = cache_fetch_keys1('annotation ncbi', annotation)
    cache_keys1_annotation_ge_seq = cache_fetch_keys1('annotation ge_seq', annotation)
    cache_keys1_annotation_small_d = cache_fetch_keys1('annotation small_d', small_d)


#
def _seq_annotation(seq):
    if irs := find_chloroplast_irs(seq, check_length=False):
        ira, irb = irs
        # To be sure!
        ira_p = ira.location.parts
        irb_p = irb.location.parts
        d = dict(length=len(seq.seq),
                 ira=[int(ira_p[0].start) - 1, int(ira_p[-1].end)],
                 irb=[int(irb_p[0].start) - 1, int(irb_p[-1].end)])
        #
        ira_s = ira.extract(seq)
        irb_s = irb.extract(seq)
        if ira.strand == irb.strand:
            irb_s = irb_s.reverse_complement()
        d.update(_irs_desc(ira_s, irb_s))
        return d
    return dict(length=len(seq.seq))


def _small_d_annotation(seq):
    # Pokupit rezultate, i ako ih ima: kreirati feature, napraviti extract, pa isto ko gore
    if res := small_d(seq):
        ira, irb = res
        seq_len = len(seq.seq)
        d = dict(length=len(seq.seq),
                 ira=[ira[0] - 1, ira[1]],
                 irb=[irb[0] - 1, irb[1]])
        #
        ira = _feature(seq, *ira, 1)
        irb = _feature(seq, *irb, -1)
        d.update(_irs_desc(ira.extract(seq), irb.extract(seq)))
        return d
    return dict(length=len(seq.seq))


def _feature(seq, s, e, strand):
    if s < e:
        return FeatureLocation(s, e, strand=strand)
    return CompoundLocation([FeatureLocation(s, len(seq.seq), strand=strand),
                             FeatureLocation(0, e, strand=strand)])


def _irs_desc(ira, irb):
    ira = str(ira.seq)
    irb = str(irb.seq)
    # print(len(ira), len(irb))
    # print('  ira', ira[:20], ira[-20:])
    # print('  irb', irb[:20], irb[-20:])
    if ira == irb:
        return dict(type='+')

    # Quite slow method. Not unbearable, but slow.
    # ToDo: possible speedup is to remove equal ends of sequences.
    print('eva', len(ira), len(irb))
    diff = difflib.SequenceMatcher(a=ira, b=irb, autojunk=False)
    opcodes = diff.get_opcodes()
    # R, I -> [num blocks, num bps]
    RM = [0, 0]
    IN = [0, 0]
    for x in opcodes:
        print(x)
        if x[0] == 'equal':
            continue
        if x[0] == 'replace':
            c, bp = RM, (x[2] - x[1])
        elif x[0] == 'delete':
            c, bp = IN, (x[2] - x[1])
        elif x[0] == 'insert':
            c, bp = IN, (x[4] - x[3])
        c[0] += 1
        c[1] += bp

    return dict(type=';'.join(f'{label}:{d[0]},{d[1]}' for d, label in ((RM, 'R'), (IN, 'I')) if d[0]),
                diff=opcodes)
