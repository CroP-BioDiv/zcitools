from collections import defaultdict
import itertools
from ..utils.entrez import Entrez


def find_uniq_features(seq, _type):
    fs = defaultdict(list)
    for idx, f in enumerate(seq.features):
        if f.type == _type and f.location:
            name = f.qualifiers['gene'][0]
            if not fs[name] or \
               ((s_f := set(f.location)) and not any(set(x.location).intersection(s_f) for _, x in fs[name])):
                fs[name].append((idx, f))
    return [y[1] for y in sorted(itertools.chain(*fs.values()))]


def extract_ncbi_comments(data, sequences_step, properties_db):
    # Notes:
    #  - sequences_step contains GenBank files collected from NCBI which contains genome info
    #  - properties_db is used to cache results dince they are static
    # Set empty fields
    n = dict((x, None) for x in ('artcle_title', 'journal', 'pubmed_id', 'first_date',
                                 'assembly_method', 'sequencing_technology', 'bio_project', 'sra_count'))
    for d in data.values():
        d.update(n)
    if not sequences_step:
        return

    prop_key = 'NCBI GenBank data'
    for seq_ident in sequences_step.all_sequences():
        if not properties_db or not (vals := properties_db.get_property(seq_ident, prop_key)):
            vals = dict()
            seq = sequences_step.get_sequence_record(seq_ident)

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
            prop_sra = 'NCBI SRA count'
            if not properties_db or not (vals_sra := properties_db.get_property(seq_ident, prop_sra)):
                vals_sra = dict()
                for x in seq.dbxrefs:  # format ['BioProject:PRJNA400982', 'BioSample:SAMN07225454'
                    if x.startswith('BioProject:'):
                        if bp := x.split(':', 1)[1]:
                            vals_sra['bio_project'] = bp
                            sra_count = Entrez().search_count('sra', term=f"{bp}[BioProject]")
                            vals_sra['sra_count'] = sra_count or None  # None means empty cell :-)
                if properties_db:
                    properties_db.set_property(seq_ident, prop_sra, vals_sra)

            vals.update(vals_sra)
            if properties_db:
                properties_db.set_property(seq_ident, prop_key, vals)

        data[seq_ident].update(vals)


