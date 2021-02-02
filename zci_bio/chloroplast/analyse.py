import os.path
from collections import defaultdict
from functools import reduce
# import multiprocessing
# from concurrent.futures import ThreadPoolExecutor
from step_project.common.table.steps import TableStep
from common_utils.file_utils import ensure_directory, write_str_in_file, write_fasta
from common_utils.properties_db import PropertiesDB
from common_utils.cache import cache
from .analyse_step_1 import SequenceDesc
from .analyse_step_calc import run_align_cmd  # , find_missing_partitions
from .utils import find_chloroplast_irs
from .constants import DEFAULT_KEEP_OFFSET
from ..utils.features import Feature
from ..utils.ncbi_taxonomy import get_ncbi_taxonomy


# ---------------------------------------------------------
# Analyse genomes
# Steps:
#  1. Collect data from steps (input annotations, previous table)
#  2. Calculate 'missing data'
#  3. Set table from collected and evaluated data
# ---------------------------------------------------------
def analyse_genomes(step_data, annotations_step):
    step = TableStep(annotations_step.project, step_data, remove_data=True)
    annotations_step.propagate_step_name_prefix(step)
    AnalyseGenomes(step, annotations_step).run()
    step.save()
    step.to_excel('chloroplast_analysis.xls')  # Store data for a check
    return step


def _get_attr(d, c):
    if '.' in c:
        for c in c.split('.'):
            d = getattr(d, c)
        return d
    return getattr(d, c)


class AnalyseGenomes:
    def __init__(self, step, annotations_step):
        self.step = step
        self.annotations_step = annotations_step
        project = step.project
        self.table_step = project.find_previous_step_of_type(self.annotations_step, 'table')
        self.sequences_step = project.find_previous_step_of_type(self.annotations_step, 'sequences')

    @property
    @cache
    def ncbi_2_max_taxid(self):
        return self.table_step.mapping_between_columns('ncbi_ident', 'max_taxid')

    @property
    @cache
    def taxid_2_ncbi(self):
        return self.table_step.mapping_between_columns('tax_id', 'ncbi_ident')

    def run(self):
        self.properties_db = PropertiesDB()

        # 1. Collect data from steps
        self.table_data = self.table_step.index_on_table('ncbi_ident')

        ncbi_rec = self.sequences_step.get_sequence_record
        data = dict((seq_ident, SequenceDesc(seq_ident, seq, ncbi_rec(seq_ident), self))
                    for seq_ident, seq in self.annotations_step._iterate_records())
        self.seq_descs = data

        # #  2. Calculate 'missing data'
        # find_missing_partitions(self.seq_descs)

        # # Set partition data
        # for seq_ident, d in data.items():
        #     d.set_parts_data()

        #  3. Set table from collected and evaluated data
        columns = [
            # tuples (dict's attribute, column name, column type)
            ('seq_ident', 'AccesionNumber', 'seq_ident'),
            ('bio_project', 'BioProject', 'str'),
            ('title', 'Title', 'str'),
            ('created_date', 'Date', 'date'),
            ('first_date', 'First date', 'date'),
            ('length', 'Length', 'int'),

            # GeSeq annotation
            ('ge_seq.num_genes_stat', 'GeSeq genes', 'str'),
            ('ge_seq.part_starts', 'GeSeq part starts', 'str'),
            ('ge_seq.part_lengths', 'GeSeq part lengths', 'str'),
            ('ge_seq.part_num_genes', 'GeSeq part genes', 'str'),
            ('ge_seq.part_orientation', 'GeSeq orientation', 'str'),
            ('ge_seq.trnF_GAA', 'GeSeq trnF-GAA', 'int'),
            # ('ge_seq.trnH_GUG', 'trnH-GUG', 'str'),

            # NCBI annotation
            ('ncbi.num_genes_stat', 'NCBI genes', 'str'),
            ('ncbi.part_starts', 'NCBI part starts', 'str'),
            ('ncbi.part_lengths', 'NCBI part lengths', 'str'),
            ('ncbi.part_num_genes', 'NCBI part genes', 'str'),
            ('ncbi.part_orientation', 'NCBI orientation', 'str'),
            ('ncbi.trnF_GAA', 'NCBI trnF-GAA', 'int'),

            #
            ('artcle_title', 'Article', 'str'),
            ('journal', 'Journal', 'str'),
            ('pubmed_id', 'PubMed', 'int'),
            ('assembly_method', 'Assembly Method', 'str'),
            ('sequencing_technology', 'Sequencing Technology', 'str'),
            ('sra_count', 'SRA count', 'int'),
        ]
        ga = _get_attr
        self.step.set_table_data(
            [[ga(d, c) for c, _, _ in columns] for seq_ident, d in sorted(data.items())],
            [(n, t) for _, n, t in columns])

        # Summary data
        some_data = next(iter(data.values()))
        sequences_step = some_data.sequences_step

        # Corrections
        offset_off = lambda o: (abs(o or 0) > DEFAULT_KEEP_OFFSET)

        def _desc_wrong_orientations(annot_attr):
            _po_count = defaultdict(int)
            for d in data.values():
                if o := getattr(d, annot_attr).part_orientation:
                    for p in o.split(','):
                        _po_count[p] += 1
            return '; '.join(f'{n}-{p}' for p, n in sorted(_po_count.items())) if _po_count else None

        summary_data = dict(
            num_genomes=len(data),
            # min_sequence_length=min(len(seq) for _, seq in sequences_step._iterate_records()),
            # max_sequence_length=max(len(seq) for _, seq in sequences_step._iterate_records()),
        )
        for a_attr in ('ge_seq', 'ncbi', 'sum'):
            summary_data[f'{a_attr}_num_irs'] = sum(1 for d in data.values() if getattr(d, a_attr).has_irs)
            summary_data[f'{a_attr}_wrong_orientations'] = sum(1 for d in data.values() if getattr(d, a_attr).part_orientation)
            summary_data[f'{a_attr}_desc_wrong_orientations'] = _desc_wrong_orientations(a_attr)
            summary_data[f'{a_attr}_wrong_offset'] = sum(1 for d in data.values() if offset_off(getattr(d, a_attr).part_offset))
            if summary_data[f'{a_attr}_wrong_offset']:
                summary_data[f'{a_attr}_wrong_offset_list'] = sorted(o for d in data.values() if offset_off(o := getattr(d, a_attr).part_offset))
            else:
                summary_data[f'{a_attr}_wrong_offset_list'] = None
            summary_data[f'{a_attr}_num_to_fix'] = sum(
                1 for d in data.values() if getattr(d, a_attr).part_orientation or offset_off(getattr(d, a_attr).part_offset))

        for a_attr in ('ge_seq', 'ncbi'):
            gks = list(getattr(some_data, a_attr).genes_stat.keys())
            min_genes = reduce(lambda x, y: dict((k, min(x[k], y[k])) for k in gks),
                               (getattr(d, a_attr).genes_stat for d in data.values()))
            max_genes = reduce(lambda x, y: dict((k, max(x[k], y[k])) for k in gks),
                               (getattr(d, a_attr).genes_stat for d in data.values()))

            for k, v in min_genes.items():
                summary_data[f'{a_attr}_genes_min_{k}'] = v
                summary_data[f'{a_attr}_genes_max_{k}'] = max_genes[k]

            p_lengths = [l for d in data.values() if (l := getattr(d, a_attr).part_lengths_all())]
            for p_idx, p_name in enumerate(('lsc', 'ir', 'ssc')):
                summary_data[f'{a_attr}_part_min_length_{p_name}'] = min(l[p_idx] for l in p_lengths)
                summary_data[f'{a_attr}_part_max_length_{p_name}'] = max(l[p_idx] for l in p_lengths)
        #
        self.step.save_summary_data(summary_data)


# ---------------------------------------------------------
# Ns in sequences
# ---------------------------------------------------------
def _ns_parts(sequences):
    # ToDo: faster
    ns = []
    for i, c in enumerate(sequences):
        if c == 'N':
            if ns and ns[-1][1] == i:
                ns[-1][1] += 1
            else:
                ns.append([i, i + 1])
    return ns


def analyse_ns(step_data, sequences_step):
    project = sequences_step.project
    step = TableStep(project, step_data, remove_data=True)

    table_step = project.find_previous_step_of_type(sequences_step, 'table')
    ncbi_2_taxid = table_step.mapping_between_columns('ncbi_ident', 'tax_id')
    taxid_2_ncbi = dict((v, k) for k, v in ncbi_2_taxid.items())
    ncbi_2_max_taxid = table_step.mapping_between_columns('ncbi_ident', 'max_taxid')
    # ncbi_2_title = table_step.mapping_between_columns('ncbi_ident', 'title')

    sequences = dict(sequences_step._iterate_records())  # Copy sequences for matching
    seq_ident_2_ns = dict((seq_ident, ns) for seq_ident, seq in sequences.items()
                          if (ns := _ns_parts(seq.seq)))
    all_taxids = set(ncbi_2_taxid[s] for s in sequences_step.all_sequences())

    if seq_ident_2_ns:
        rows = []
        ncbi_tax = get_ncbi_taxonomy()
        # with ThreadPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
        for seq_ident, ns in seq_ident_2_ns.items():
            row = [seq_ident, len(sequences[seq_ident]), len(ns), sum(b - a for a, b in ns), None, None]
            rows.append(row)

            taxid = ncbi_2_taxid[seq_ident]
            search_in = set(all_taxids)
            search_in.discard(taxid)
            close_taxids = ncbi_tax.find_close_taxids(taxid, ncbi_2_max_taxid[seq_ident], search_in)
            if not close_taxids:
                print(f"Warning: sequence {seq_ident} doesn't have close relative in accession set!")
                continue

            f_dir = step.step_file('repair_ns', seq_ident)
            ensure_directory(f_dir)
            # seq_fasta = os.path.join(f_dir, f"{seq_ident}.fa")
            # write_fasta(seq_fasta, [(seq_ident, seq_data['_seq'].seq)])

            # executor.submit(_run_manage_ns, seq_ident, sequences, ns, f_dir, [taxid_2_ncbi[t] for t in close_taxids])
            _run_manage_ns(seq_ident, sequences, ns, f_dir, [taxid_2_ncbi[t] for t in close_taxids])
        #
        columns = [('seq_ident', 'seq_ident'), ('length', 'int'),
                   ('num_ns_parts', 'int'), ('ns_length', 'int'),
                   ('close_seq_idents', 'str'), ('fix', 'str')]
        step.set_table_data(rows, columns)
    else:
        step.set_columns([('seq_ident', 'seq_ident')])  # Dummy table

    step.save()
    return step


def _run_manage_ns(seq_ident, sequences, ns, f_dir, close_seq_idents, match_sides=30):
    # Store close sequences to search in
    close_filename = os.path.join(f_dir, 'close_sequences.fa')
    query_filename = os.path.join(f_dir, 'query_sequences.fa')
    write_fasta(close_filename, [(c, sequences[c].seq) for c in close_seq_idents[:5]])

    # Find parts to extract
    seq = sequences[seq_ident]
    seq_length = len(seq)
    a, b = ns[0]
    m = max(match_sides, (b - a) * 2)
    parts = [Feature(seq_length, interval=(a - m, b + m))]
    for i, (a, b) in enumerate(ns[1:]):
        m = max(match_sides, (b - a) * 2)
        f = Feature(seq_length, interval=(a - m, b + m))
        if parts[-1].intersects(f):
            parts[-1] = Feature(seq_length, interval=(parts[-1].real_start, f.real_end))
        else:
            parts.append(Feature(seq_length, interval=(a - m, b + m)))

    # Extract query data
    write_fasta(query_filename, [(f"{f.real_start}_{f.real_end}", f.extract(seq).seq) for f in parts])
    align = run_align_cmd(close_filename, query_filename, 'result')
    # ToDo: ...
