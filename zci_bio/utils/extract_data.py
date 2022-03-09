import os
from functools import wraps
import time
from datetime import datetime
from Bio.SeqFeature import FeatureLocation, CompoundLocation
from Bio import SeqIO
from Bio.SeqUtils import GC
from zci_bio.utils.entrez import Entrez
from zci_bio.utils.diff_sequences import diff_check_memory
from zci_bio.chloroplast.utils import find_chloroplast_irs, ir_loc, cycle_distance_lt


def with_seq(func):
    @wraps(func)
    def func_wrapper(self, seq_ident=None, seq=None, key=None, seq_filename=None):
        if seq_filename:
            seq = SeqIO.read(seq_filename, 'genbank')  # ToDo: Type?
        if seq is None:
            assert seq_ident
            seq = self.sequences_step.get_sequence_record(seq_ident, cache=False)
        if not seq_ident:
            seq_ident = seq.name.split('.')[0]
        return func(self, seq, seq_ident, key)
    return func_wrapper


def with_seq_ident(func):
    @wraps(func)
    def func_wrapper(self, seq_ident=None, seq=None, key=None, seq_filename=None):
        if not seq_ident:
            if seq_filename:
                seq_ident = os.path.splitext(os.path.basename(seq_filename))[0]
            else:
                seq_ident = seq.name
        return func(self, seq_ident, key, seq_filename)
    return func_wrapper


class ExtractData:
    # Extract/fetch additional data from NCBI GenBank file or by Entrez searches.
    # Idea is to cache that data
    # Two set of methods:
    #  - set of methods that extract or fetch data
    #  - set of methods use upper methods and do caching.
    def __init__(self, properties_db=None, sequences_step=None, look_for_diff=True):
        self.properties_db = properties_db
        self.sequences_step = sequences_step
        self.look_for_diff = look_for_diff

    def _seq_filename(self, seq_ident, seq_filename):
        if seq_filename:
            return os.path.abspath(seq_filename)
        return os.path.abspath(self.sequences_step.get_sequence_filename(seq_ident))

    def _seq_record(self, seq_ident, seq_filename):
        if seq_filename:
            return SeqIO.read(seq_filename, 'genbank')  # ToDo: Type?
        return self.sequences_step.get_sequence_record(seq_ident, cache=False)

    # Extract and fetch methods
    @with_seq
    def genbank_data(self, seq, seq_ident, key):
        annotations = seq.annotations

        genes = [f for f in seq.features if f.type == 'gene' and 'gene' in f.qualifiers]
        vals = dict(length=len(seq.seq),
                    gc=GC(seq.seq),
                    num_genes=len(genes),
                    num_unique_genes=len(set(f.qualifiers['gene'][0] for f in genes)))
        if not_dna := [i for i, c in enumerate(str(seq.seq)) if c not in 'ATCG']:
            vals['not_dna'] = not_dna

        vals.update((k, v) for k in ('organism', 'sequence_version') if (v := annotations.get(k)))
        if v := annotations.get('date'):
            vals['update_date'] = datetime.strptime(v, '%d-%b-%Y').date()

        refs = annotations['references']
        if refs[0].title != 'Direct Submission':
            vals['article_title'] = refs[0].title
            vals['journal'] = refs[0].journal
            vals['authors'] = refs[0].authors
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
        return self._from_features(seq, seq_ident, key, find_chloroplast_irs(seq, check_length=False))

    @with_seq
    def airpg(self, seq, seq_ident, key):
        from zci_bio.chloroplast.irs.airpg import airpg
        return self._from_features(seq, seq_ident, key, airpg(seq, ret_features=True))

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
    def chloroplot(self, seq_ident, key, seq_filename):
        from zci_bio.chloroplast.irs.chloroplot import chloroplot as chloroplot_ann
        return self._from_indices(
            self._seq_record(seq_ident, seq_filename), seq_ident, key,
            lambda: chloroplot_ann(self._seq_filename(seq_ident, seq_filename)))

    @with_seq_ident
    def pga(self, seq_ident, key, seq_filename):
        from zci_bio.chloroplast.irs.pga import pga
        return self._from_indices(
            self._seq_record(seq_ident, seq_filename), seq_ident, key,
            lambda: pga(self._seq_filename(seq_ident, seq_filename)))

    @with_seq_ident
    def pga_sb(self, seq_ident, key, seq_filename):
        return self._self_blast('pga', seq_ident, key, seq_filename)

    @with_seq_ident
    def plann(self, seq_ident, key, seq_filename):
        from zci_bio.chloroplast.irs.plann import plann
        return self._from_indices(
            self._seq_record(seq_ident, seq_filename), seq_ident, key,
            lambda: plann(self._seq_filename(seq_ident, seq_filename)))

    @with_seq_ident
    def plann_sb(self, seq_ident, key, seq_filename):
        return self._self_blast('plann', seq_ident, key, seq_filename)

    @with_seq_ident
    def org_annotate(self, seq_ident, key, seq_filename):
        from zci_bio.chloroplast.irs.org_annotate import org_annotate
        return self._from_indices(
            self._seq_record(seq_ident, seq_filename), seq_ident, key,
            lambda: org_annotate(self._seq_filename(seq_ident, seq_filename)))

    @with_seq_ident
    def zci_method(self, seq_ident, key, seq_filename):
        try:
            from irs_annotate import annotate_filename
        except ImportError:
            print('\nError: irs_annotate.py is not part of ZCI tools project.\nSet PYTHONPATH!!!\n\n')
            raise
        return self._from_indices(
            self._seq_record(seq_ident, seq_filename), seq_ident, key,
            lambda: annotate_filename(self._seq_filename(seq_ident, seq_filename)))

    #
    def _small_d_annotation(self, seq, seq_ident, key, no_prepend_workaround=True, no_dna_fix=True):
        from zci_bio.chloroplast.irs.small_d import small_d
        return self._from_indices(
            seq, seq_ident, key,
            lambda: small_d(seq, no_prepend_workaround=no_prepend_workaround, no_dna_fix=no_dna_fix))

    def _self_blast(self, variant, seq_ident, key, seq_filename):
        from zci_bio.chloroplast.irs.self_blast import self_blast
        return self._from_indices(
            self._seq_record(seq_ident, seq_filename), seq_ident, key,
            lambda: self_blast(variant, self._seq_filename(seq_ident, seq_filename)))

    #
    def _from_features(self, seq, seq_ident, key, irs):
        if irs:
            ira, irb = irs
            d = dict(length=len(seq.seq),
                     ira=ir_loc(ira.location.parts),
                     irb=ir_loc(irb.location.parts))
            # Sequences NC_035878 and NC_045527 have molecule_type = RNA?!?!
            seq.annotations['molecule_type'] = 'DNA'
            ira_s = ira.extract(seq)
            irb_s = irb.extract(seq)
            if ira.strand == irb.strand:
                irb_s = irb_s.reverse_complement()
            if desc := self._irs_desc(seq, seq_ident, key, ira_s, irb_s, d):
                d.update(desc)
                return d
            return None
        print(f'{key.split(" ", 1)[1]} {seq_ident}: no IRs')
        return dict(length=len(seq.seq))

    def _from_indices(self, seq, seq_ident, key, irs_method):
        start = time.perf_counter()
        irs = irs_method()
        d = dict(length=len(seq.seq), calculation_time=_format_time(time.perf_counter() - start))
        if irs:
            ira, irb = irs
            d['ira'] = [ira[0], ira[1]]
            d['irb'] = [irb[0], irb[1]]
            # Sequences NC_035878 and NC_045527 have molecule_type = RNA?!?!
            seq.annotations['molecule_type'] = 'DNA'
            ira = self._feature(seq, *ira, 1)
            irb = self._feature(seq, *irb, -1)
            if desc := self._irs_desc(seq, seq_ident, key, ira.extract(seq), irb.extract(seq), d):
                d.update(desc)
                return d
            return None
        print(f'{key.split(" ", 1)[1]} {seq_ident}: no IRs')
        return d

    def _feature(self, seq, s, e, strand):
        if s < e:
            return FeatureLocation(s, e, strand=strand)
        if strand == 1:
            return CompoundLocation([FeatureLocation(s, len(seq.seq), strand=1),
                                     FeatureLocation(0, e, strand=1)])
        return CompoundLocation([FeatureLocation(0, e, strand=-1),
                                 FeatureLocation(s, len(seq.seq), strand=-1)])

    def _irs_desc(self, seq, seq_ident, key, ira, irb, irs_d):
        ira = str(ira.seq)
        irb = str(irb.seq)
        print(f'{key.split(" ", 1)[1]} {seq_ident}: {len(ira)}, {len(irb)}: {irs_d["ira"]}, {irs_d["irb"]}')
        # if len(ira) > 30000 or len(irb) > 30000:
        #     return  # Hack: for now
        if ira == irb:
            return self._irs_desc_add_1_store(irs_d, dict(type='+'), seq)

        # Check is same IRs region already inspected.
        if self.look_for_diff:
            seq_ident = seq.name.split('.')[0]
            for _, data in self.properties_db.get_properties_key2_like(seq_ident, 'annotation %').items():
                if irs_d['ira'] == data.get('ira') and irs_d['irb'] == data.get('irb'):
                    print(f'  found diff: {data["type"]}')
                    return dict((k, v) for k, v in data.items()
                                if k in ('type', 'diff', 'ir_lengths', 'diff_len', 'max_indel_length', 'not_dna',
                                         'ssc_length', 'lsc_length'))

        # Problem:
        # Original NCBI annotation of NC_031898, has IR lengths 23595 and 74817.
        # IRn is annotated as complement(10932..85748). It should be 85751..109342.
        # Copy/paste error?
        if len(irb) > 3 * len(ira) / 2 or len(ira) > 3 * len(irb) / 2:
            return self._irs_desc_add_1_store(irs_d, dict(type='???'), seq)
        print(f'  calculate diff')
        diff = diff_check_memory(ira, irb)
        return self._irs_desc_add_1_store(irs_d, dict(type=diff.in_short(), diff=diff.get_opcodes()), seq)

    #
    def _irs_desc_add_1_store(self, irs_d, d, seq):
        tmp_d = dict(d)
        tmp_d.update(irs_d)
        not_dna = sum(1 for c in str(seq.seq) if c not in 'ATCG')
        d_1 = self._irs_desc_add_1(tmp_d, seq, not_dna, store=True)
        d.update(d_1)
        return d

    def _irs_desc_add_1(self, irs_d, seq, not_dna, store=False, ira=None, irb=None):
        # Addition to method _irs_desc().
        seq_length = irs_d['length']
        ira_l = cycle_distance_lt(*irs_d['ira'], seq_length)
        irb_l = cycle_distance_lt(*irs_d['irb'], seq_length)
        d = dict(ir_lengths=[ira_l, irb_l], diff_len=abs(ira_l - irb_l))

        a_s, a_e = irs_d['ira']
        b_s, b_e = irs_d['irb']
        d['ssc_length'], d['lsc_length'] = sorted([(b_s - a_e) % seq_length, (a_s - b_e) % seq_length])

        if diff := irs_d.get('diff'):
            # ToDo: <a>,equal,<a>,equal by few bases
            max_diff = max(max(src_e - src_s, tgt_e - tgt_s)
                           for action, src_s, src_e, tgt_s, tgt_e in diff if action != 'equal')
            d['max_indel_length'] = max_diff

        if not_dna:
            if ira is None and seq:
                ira = self._feature(seq, *irs_d['ira'], 1)
                irb = self._feature(seq, *irs_d['irb'], -1)
            if ira and irb:
                d['not_dna'] = [sum(1 for c in str(ira.extract(seq).seq) if c not in 'ATCG'),
                                sum(1 for c in str(irb.extract(seq).seq) if c not in 'ATCG')]

        if store:
            irs_d.update(d)
        return d

    def ir_data(self, cdb):
        # Applies _irs_desc_add_1() to stored IR data
        like = 'annotation %'
        seq_idents = self.properties_db.get_keys1_key2_like(like)
        for seq_ident in sorted(seq_idents):
            annots = self.properties_db.get_properties_key2_like(seq_ident, like)

            # If all annotations are without IRs, than there is nothing to do
            if all('ira' not in data for data in annots.values()):
                print(f'{seq_ident}: no method located IRs')
                continue

            # Note: sequence is needed only for finding do IRs contain not DNA bases.
            # First check number of not DNA bases in sequence (NCBI GenBank data) than load sequence if needed
            seq = None

            # First check stored data, and if there is no one, than calculate and cache it
            if not (gb_data := self.properties_db.get_property(seq_ident, 'NCBI GenBank data')):
                # Fetch sequence
                if seq_io := cdb.get_record_stringIO(seq_ident):
                    seq = SeqIO.read(seq_io, 'genbank')
                    gb_data = self.cache_genbank_data(seq=seq, seq_ident=seq_ident)
                else:
                    print(f'Warning: no sequence {seq_ident}!!!')

            not_dna = gb_data.get('not_dna', 0) if gb_data else 0
            if not_dna and not seq:
                if seq_io := cdb.get_record_stringIO(seq_ident):
                    seq = SeqIO.read(seq_io, 'genbank')
                else:
                    not_dna = 0  # Not possible to check not DNA bases inside IRs

            # Proci po anotacijama i ekstrahirati feature IR-ova
            # _irs_desc(seq, seq_ident, key, ira, irb, irs_d)
            cached_data = dict()  # (ira, irb) -> found data
            print(f'{seq_ident}: ', end='', flush=True)
            for a_method, irs_d in annots.items():
                print('.', end='', flush=True)
                if 'ira' not in irs_d:
                    continue

                # Additional data already stored!
                if 'irs_lengths' in irs_d:
                    continue

                c_key = (tuple(irs_d['ira']), tuple(irs_d['irb']))
                if not (new_data := cached_data.get(c_key)):
                    new_data = self._irs_desc_add_1(irs_d, seq, not_dna, store=False)
                    cached_data[c_key] = new_data

                # Store data
                if new_data:
                    irs_d.update(new_data)
                    self.properties_db.set_property(seq_ident, a_method, irs_d)
            print()

    def diff_stats(self, seq_idents):
        from zci_bio.chloroplast.irs.diff_stats import analyse_diff
        like = 'annotation %'
        for seq_ident in sorted(seq_idents):
            print(seq_ident)
            # print(f'{seq_ident}: ', end='', flush=True)
            annots = self.properties_db.get_properties_key2_like(seq_ident, like)
            for a_method, irs_d in sorted(annots.items(), key=lambda x: max(x[1].get('ir_lengths', [0]))):
                # Test
                if a_method.endswith('_sb') or 'small_d' in a_method:
                    continue
                # print('.', end='', flush=True)
                if 'diff' in irs_d:  # 'max_indel' not in irs_d and
                    print('  ', a_method.split()[1], max(irs_d['ir_lengths']), irs_d['ira'], irs_d['irb'])
                    new_data = analyse_diff(irs_d)
                    print('    ', ' '.join(f'{n}{m.lower()}' for m, n in new_data['diff_simple']))
                    # print('  ', ' '.join(map(str, new_data['equal_lengths'])))
                    self.properties_db.set_property(seq_ident, a_method, {**irs_d, **new_data})
                elif 'ira' in irs_d:
                    print('  *', a_method.split()[1], max(irs_d['ir_lengths']), irs_d['ira'], irs_d['irb'])
                else:
                    print('  ', a_method.split()[1], '-')
            # print()

    def lsc_ssc_lengths(self, seq_idents):
        like = 'annotation %'
        for seq_ident in sorted(seq_idents):
            print(f'{seq_ident}: ', end='', flush=True)
            annots = self.properties_db.get_properties_key2_like(seq_ident, like)
            for a_method, irs_d in sorted(annots.items(), key=lambda x: max(x[1].get('ir_lengths', [0]))):
                if 'ira' in irs_d and 'lsc_length' not in irs_d:
                    print('.', end='', flush=True)
                    seq_l = irs_d['length']
                    a_s, a_e = irs_d['ira']
                    b_s, b_e = irs_d['irb']
                    irs_d['ssc_length'], irs_d['lsc_length'] = sorted([(b_s - a_e) % seq_l, (a_s - b_e) % seq_l])
                    self.properties_db.set_property(seq_ident, a_method, irs_d)
            print()

    @with_seq
    def not_ns(self, seq, seq_ident, key):
        if not_ns := [i for i, c in enumerate(str(seq.seq)) if c not in 'ATCGN']:
            print(seq_ident, len(not_ns))


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
        ('annotation airpg', 'annotation_airpg', ExtractData.airpg),
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
        ('annotation zci', 'annotation_zci', ExtractData.zci_method),
        ):
    # Caching
    # Interface: method_name(seq_ident=None, seq=None)
    setattr(ExtractData, f'cache_{m_sufix}', cache_fetch(key, cls_method))
    # Bulk fetch
    # Interface: method_name(seq_idents, seq_step=None)
    setattr(ExtractData, f'cache_keys1_{m_sufix}', cache_fetch_keys1(key, cls_method))


def _format_time(seconds):
    if seconds < 0.001:
        return f'{round(seconds * 1000000)}ns'
    if seconds < 0.2:
        return f'{round(seconds * 1000)}ms'
    if seconds < 10:
        return f'{round(seconds, 2)}s'
    if seconds > 3600:
        hours = int(seconds // 3600)
        rest = seconds - hours * 3600
        minutes = int(rest // 60)
        return f'{hours}h {minutes}m'
    if seconds > 60:
        minutes = int(seconds // 60)
        rest = seconds - minutes * 60
        return f'{minutes}m {int(rest)}s'
    return f'{round(seconds, 1)}s'


if __name__ == '__main__':
    import argparse
    from common_utils.properties_db import PropertiesDB
    from common_utils.exceptions import ZCItoolsValueError

    _arg_2_method_key = dict(
        genbank_data=('genbank_data', 'NCBI GenBank data'),
        ncbi=('annotation', 'annotation ncbi'),
        ge_seq=('annotation', 'annotation ge_seq'),
        chloe=('annotation', 'annotation chloe'),
        airpg=('airpg', 'annotation airpg'),
        small_d=('small_d', 'annotation small_d'),
        small_d_P=('small_d_P', 'annotation small_d_P'),
        small_d_D=('small_d_D', 'annotation small_d_D'),
        small_d_all=('small_d_all', 'annotation small_d_all'),
        chloroplot=('chloroplot', 'annotation chloroplot'),
        pga=('pga', 'annotation pga'),
        pga_sb=('pga_sb', 'annotation pga_sb'),
        plann=('plann', 'annotation plann'),
        plann_sb=('plann_sb', 'annotation plann_sb'),
        org_annotate=('org_annotate', 'annotation org_annotate'),
        zci=('zci_method', 'annotation zci'),
        ir_data=None,
        diff_stats=None,
        lsc_ssc_lengths=None,
        not_ns=('not_ns', None),
    )
    parser = argparse.ArgumentParser(description="""
Calls ExtractData method on given sequence

Examples:
python3 extract_data.py pga <path_to_sequence_filename>
Make PGA IRs location on sequence given with a filename.

python3 extract_data.py plann -c sequences <sequence_accession_number>
Make Plann IRs location on sequence in stored in Common DB.
File is searched with path 'sequences' and given ident (accession number).
""", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('method_name', choices=sorted(_arg_2_method_key.keys()), help='Method name')

    # Process one sequence (file or cached)
    parser.add_argument(
        '-c', '--common-db',
        help='Common DB path. If specified, seq_filename is used as ident inside Common DB specified with given path.')
    parser.add_argument('-f', '--seq-filename', help='Sequence filename')
    parser.add_argument('-l', '--look-for-diff', action='store_true', help='Search for diff in properties DB.')

    # Process whole
    parser.add_argument('-C', '--process-common-db', action='store_true', help='Process Common DB files')
    parser.add_argument('-s', '--starts-with', help='Common DB filename starts with')
    parser.add_argument('-F', '--force-update', action='store_true', help='Update all common DB entities, not only missing in properties')
    params = parser.parse_args()

    #
    if params.method_name == 'ir_data':
        from common_utils.common_db import CommonDB
        ed = ExtractData(properties_db=PropertiesDB())
        ed.ir_data(CommonDB.get_zci_db(tuple(params.common_db.split('/'))))

    elif params.method_name in ('diff_stats', 'lsc_ssc_lengths'):
        from common_utils.common_db import CommonDB
        seq_idents = set()
        if params.common_db:
            cdb = CommonDB.get_zci_db(tuple(params.common_db.split('/')))
            seq_idents.update(p[-1] for p in cdb.get_all_record_ident(startswith=params.starts_with))
        if params.seq_filename:
            seq_idents.add(params.seq_filename)
        ed = ExtractData(properties_db=PropertiesDB())
        getattr(ed, params.method_name)(seq_idents)

    elif params.process_common_db:
        if not params.common_db:
            raise ZCItoolsValueError('Common DB is not set!!!')

        import tempfile
        from common_utils.common_db import CommonDB
        cdb = CommonDB.get_zci_db(tuple(params.common_db.split('/')))
        properties_db = PropertiesDB()
        ed = ExtractData(properties_db=properties_db, look_for_diff=params.look_for_diff)
        method_name, key2 = _arg_2_method_key[params.method_name]

        # Find idents to process
        seq_idents_2_paths = dict((p[-1], p) for p in cdb.get_all_record_ident(startswith=params.starts_with))
        if not params.force_update:
            use_idents = properties_db.not_stored_keys1(set(seq_idents_2_paths.keys()), key2)
            seq_idents_2_paths = dict((k, v) for k, v in seq_idents_2_paths.items() if k in use_idents)

        # Process sequences
        tmp_dir = tempfile.gettempdir()
        for seq_ident, cds_path in sorted(seq_idents_2_paths.items()):
            if key2:
                print(seq_ident)
            seq_filename = cdb.get_record(cds_path, tmp_dir, info=False)
            data = getattr(ed, method_name)(seq_ident=seq_ident, key=key2, seq_filename=seq_filename)
            if data and key2:
                properties_db.set_property(seq_ident, key2, data)
            os.remove(seq_filename)

    else:
        seq_filename = params.seq_filename
        tmp_file = None

        if params.common_db:
            import tempfile
            from common_utils.common_db import CommonDB
            cdb = CommonDB.get_zci_db(tuple(params.common_db.split('/')))
            seq_filename = cdb.get_record(tuple(params.seq_filename.split('/')), tempfile.gettempdir(), info=True)
            tmp_file = seq_filename

        #
        if seq_filename:
            ed = ExtractData(properties_db=PropertiesDB(), look_for_diff=params.look_for_diff)
            method_name, key = _arg_2_method_key[params.method_name]
            ir_desc = getattr(ed, method_name)(key=key, seq_filename=seq_filename)
            print(ir_desc)

            if tmp_file:
                os.remove(tmp_file)
        else:
            print('No sequence file set!!!')
