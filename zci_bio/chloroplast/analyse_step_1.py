from collections import defaultdict
import itertools
from datetime import datetime
from .utils import find_chloroplast_partition, create_chloroplast_partition, \
    chloroplast_parts_orientation, trnF_GAA_start  # , trnH_GUG_start
from ..utils.entrez import Entrez
from ..utils.helpers import fetch_from_properties_db
from ..utils.features import Feature
from common_utils.cache import cache


class SequenceDesc:
    def __init__(self, seq_ident, seq, analyse):  # table_data, sequence_step,
        self.seq_ident = seq_ident
        self._seq = seq
        self._analyse = analyse
        self._table_data = analyse.table_data
        #
        self.length = len(seq)
        self._genes = find_uniq_features(seq, 'gene')
        self._cds = find_uniq_features(seq, 'CDS')
        self._partition = find_chloroplast_partition(seq)
        self._parts_data = _PartsDesc(self._partition, self._genes, self.length) if self._partition else None
        # Set in step 2. (Calculate missing data)
        self.irs_took_from = None
        self.irs_took_reason = None
        self.part_orientation = None

        # Extract data from NCBI GenBank files (comments)
        self.extract_ncbi_comments(analyse.properties_db)

    # Environment
    step = property(lambda self: self._analyse.step)
    annotations_step = property(lambda self: self._analyse.annotations_step)
    sequences_step = property(lambda self: self._analyse.sequences_step)
    table_step = property(lambda self: self._analyse.table_step)
    #
    num_genes = property(lambda self: len(self._genes))
    num_cds = property(lambda self: len(self._cds))
    taxid = property(lambda self: self._table_data.get_cell(self.seq_ident, 'tax_id'))
    title = property(lambda self: self._table_data.get_cell(self.seq_ident, 'title'))
    created_date = property(lambda self: self._table_data.get_cell(self.seq_ident, 'create_date'))
    parts_desc = property(lambda self: self._parts_data.desc_str() if self._parts_data else None)
    part_starts = property(lambda self: self._parts_data.starts_str() if self._parts_data else None)
    part_lengths = property(lambda self: self._parts_data.lengths_str() if self._parts_data else None)
    part_num_genes = property(lambda self: self._parts_data.num_genes_str() if self._parts_data else None)
    part_offset = property(lambda self: self._parts_data.offset if self._parts_data else None)
    part_trnH_GUG = property(lambda self: self._parts_data.trnH_GUG if self._parts_data else None)

    @property
    @cache
    def trnF_GAA(self):
        return trnF_GAA_start(self._seq, self._partition)

    @property
    def trnH_GUG(self):
        if ret := trnH_GUG_start(self._seq, self._partition):
            if 'lsc_offset' in ret:
                return f"{ret['lsc_offset']} : {ret['strategy']}"
            return f"{ret['zero_offset']} : rc" if ret['reverse'] else str(ret['zero_offset'])

    def part_lengths_all(self):
        return self._parts_data.parts_length() if self._parts_data else None

    _ncbi_comment_fields = dict(
        (x, None) for x in ('artcle_title', 'journal', 'pubmed_id', 'first_date',
                            'assembly_method', 'sequencing_technology', 'bio_project', 'sra_count'))
    _key_genbank_data = 'NCBI GenBank data'
    _key_sra_count = 'NCBI SRA count'

    def extract_ncbi_comments(self, properties_db):
        # Notes:
        #  - sequences_step contains GenBank files collected from NCBI which contains genome info
        #  - properties_db is used to cache results dince they are static
        # Set empty fields
        self.__dict__.update(self._ncbi_comment_fields)

        self.__dict__.update(fetch_from_properties_db(
            properties_db, self.seq_ident, self._key_genbank_data,
            self._genbank_data, self.seq_ident, properties_db))

    def _genbank_data(self, seq_ident, properties_db):
        vals = dict()
        seq = self.sequences_step.get_sequence_record(seq_ident)

        refs = seq.annotations['references']
        if refs[0].title != 'Direct Submission':
            vals['artcle_title'] = refs[0].title
            vals['journal'] = refs[0].journal
            if refs[0].pubmed_id:
                vals['pubmed_id'] = int(refs[0].pubmed_id)
        if refs[-1].title == 'Direct Submission':
            # ToDo: re ...
            vals['first_date'] = datetime.strptime(
                refs[-1].journal.split('(', 1)[1].split(')', 1)[0], '%d-%b-%Y').date()

        if (sc := seq.annotations.get('structured_comment')) and \
           (ad := sc.get('Assembly-Data')):
            vals['assembly_method'] = ad.get('Assembly Method')
            vals['sequencing_technology'] = ad.get('Sequencing Technology')

        #
        vals.update(fetch_from_properties_db(
            properties_db, self.seq_ident, self._key_sra_count,
            self._sra_count_data, seq))
        return vals

    def _sra_count_data(self, seq):
        vals_sra = dict()
        for x in seq.dbxrefs:  # format ['BioProject:PRJNA400982', 'BioSample:SAMN07225454'
            if x.startswith('BioProject:'):
                if bp := x.split(':', 1)[1]:
                    vals_sra['bio_project'] = bp
                    sra_count = Entrez().search_count('sra', term=f"{bp}[BioProject]")
                    vals_sra['sra_count'] = sra_count or None  # None means empty cell :-)
        return vals_sra

    #
    def set_parts_data(self):
        if self._partition:  # Take better one
            # Part orientation
            orient = chloroplast_parts_orientation(self._seq, self._partition, self._genes)
            if ppp := [p for p in ('lsc', 'ira', 'ssc') if not orient[p]]:
                self.part_orientation = ','.join(ppp)

    def set_took_part(self, ira, irb, transfer_from, reason):
        assert not self._partition, "For now there should not be original IRs!"
        self._partition = create_chloroplast_partition(self.length, ira, irb, in_interval=True)
        self._parts_data = _PartsDesc(self._partition, self._genes, self.length)
        self.irs_took_from = transfer_from
        self.irs_took_reason = reason


def _parts_desc(parts, genes, length):
    # Formats:
    #  part starts     : '<lsc start>, <ira start>, <ssc strart>, <irb start>'
    #  part lengths    : '<lsc_length>, <ir length>, <ssc length>'
    #  number of genes : '<lsc start>, <ira start>, <ssc strart>, <irb start>'
    lsc, ira, ssc, irb = [parts[p] for p in ('lsc', 'ira', 'ssc', 'irb')]
    lsc_s = seq_offset(length, lsc.real_start)
    part_starts = ', '.join([str(lsc_s)] + [str(p.real_start) for p in (ira, ssc, irb)])
    part_lens = ', '.join(str(len(p)) for p in (lsc, ira, ssc))
    part_starts = f'{part_starts}\n{part_lens}'
    if genes is None:
        return part_starts, part_lens

    # Number of genes in parts
    part_genes = parts.put_features_in_parts([Feature(length, feature=f) for f in genes])
    gs = ', '.join(str(len(part_genes[p])) for p in ('lsc', 'ssc', 'ira', 'irb'))
    return part_starts, part_lens, gs


class _PartsDesc:
    def __init__(self, parts, genes, seq_length):
        self.seq_length = seq_length
        self.oriented = [parts[p] for p in ('lsc', 'ira', 'ssc', 'irb')]
        self.lsc, self.ira, self.ssc, self.irb = self.oriented
        self.offset = seq_offset(seq_length, self.lsc.real_start)
        if genes is None:
            self.part_genes = None
            self.offset = None
            self.trnH_GUG = None
        else:
            part_genes = parts.put_features_in_parts([Feature(seq_length, feature=f) for f in genes])
            self.part_genes = [part_genes[p.name] for p in self.oriented]
            self.offset = seq_offset(seq_length, self.lsc.real_start)
            self.trnH_GUG = trnH_GUG_offset(seq_length, self.offset, genes)

    def _in_k(self, num):
        s = str(round(num / 1000, 1))
        return s[:-2] if s.endswith('.0') else s

    @cache
    def starts_str(self):
        return ', '.join([str(self.offset)] + [str(p.real_start) for p in self.oriented[1:]])

    @cache
    def lengths_str(self):
        return ', '.join(self._in_k(len(p)) for p in self.oriented[:-1])

    @cache
    def num_genes_str(self):
        if self.part_genes:
            return ', '.join(str(len(p)) for p in self.part_genes)

    def parts_length(self):
        return [len(p) for p in self.oriented]


def seq_offset(seq_length, offset):
    o2 = offset - seq_length
    return offset if offset <= abs(o2) else o2


def trnH_GUG_offset(seq_length, lsc_start, genes):
    features = [f for f in genes if f.qualifiers['gene'][0] == 'trnH-GUG']
    if features:
        rel_to_lsc = [((f.location.start - lsc_start) % seq_length) for f in features]
        return min(seq_offset(seq_length, d) for d in rel_to_lsc)


def find_uniq_features(seq, _type):
    fs = defaultdict(list)
    for idx, f in enumerate(seq.features):
        if f.type == _type:
            name = f.qualifiers['gene'][0]
            if not fs[name] or \
               ((s_f := set(f.location)) and not any(set(x.location).intersection(s_f) for _, x in fs[name])):
                fs[name].append((idx, f))
    return [y[1] for y in sorted(itertools.chain(*fs.values()))]
