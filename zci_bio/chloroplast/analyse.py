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
from .analyse_step_calc import run_align_cmd, find_missing_partitions
from .utils import find_chloroplast_partition, find_chloroplast_irs
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

        data = dict((seq_ident, SequenceDesc(seq_ident, seq, self))
                    for seq_ident, seq in self.annotations_step._iterate_records())
        self.seq_descs = data

        #  2. Calculate 'missing data'
        find_missing_partitions(self.seq_descs)

        # Set partition data
        for seq_ident, d in data.items():
            d.set_parts_data()

        #  3. Set table from collected and evaluated data
        columns = [
            # tuples (dict's attribute, column name, column type)
            ('seq_ident', 'AccesionNumber', 'seq_ident'),
            ('bio_project', 'BioProject', 'str'),
            ('title', 'Title', 'str'),
            ('created_date', 'Date', 'date'),
            ('first_date', 'First date', 'date'),
            ('length', 'Length', 'int'),
            # ('num_genes', 'Genes', 'int'),
            ('num_genes_stat', 'Genes', 'str'),
            # ('num_cds', 'CDS', 'int'),
            ('irs_took_from', 'IRS took', 'seq_ident'),
            ('part_starts', 'Part starts', 'str'),
            ('part_lengths', 'Part lengths', 'str'),
            ('part_num_genes', 'Part genes', 'str'),
            ('part_orientation', 'Orientation', 'str'),
            ('trnF_GAA', 'trnF-GAA', 'int'),
            # ('trnH_GUG', 'trnH-GUG', 'str'),
            ('artcle_title', 'Article', 'str'),
            ('journal', 'Journal', 'str'),
            ('pubmed_id', 'PubMed', 'int'),
            ('assembly_method', 'Assembly Method', 'str'),
            ('sequencing_technology', 'Sequencing Technology', 'str'),
            ('sra_count', 'SRA count', 'int'),
        ]
        self.step.set_table_data(
            [[getattr(d, c) for c, _, _ in columns] for seq_ident, d in sorted(data.items())],
            [(n, t) for _, n, t in columns])

        # Statistics
        some_data = next(iter(data.values()))
        sequences_step = some_data.sequences_step
        ncbi_with = [seq_ident for seq_ident, seq in sequences_step._iterate_records()
                     if find_chloroplast_irs(seq, check_size=False)]  # Note: MCBI IRs are used without length check!
        ge_seq_with = [seq_ident for seq_ident, seq in self.annotations_step._iterate_records()
                       if find_chloroplast_irs(seq, check_size=True)]
        p_lengths = [l for d in data.values() if (l := d.part_lengths_all())]
        #
        gks = list(some_data._genes_stat.keys())
        min_genes = reduce(lambda x, y: dict((k, min(x[k], y[k])) for k in gks), (d._genes_stat for d in data.values()))
        max_genes = reduce(lambda x, y: dict((k, max(x[k], y[k])) for k in gks), (d._genes_stat for d in data.values()))

        # Corrections
        offset_off = lambda o: (abs(o or 0) > DEFAULT_KEEP_OFFSET)
        parts_orient_count = defaultdict(int)
        for d in self.seq_descs.values():
            if d.part_orientation:
                for p in d.part_orientation.split(','):
                    parts_orient_count[p] += 1
        parts_orient_count = '; '.join(f'{n} ({p})' for p, n in sorted(parts_orient_count.items()))
        if parts_orient_count:
            parts_orient_count = f' ({parts_orient_count})'

        summary = f"""Statistics:

Minimum sequence length: {min(len(seq) for _, seq in sequences_step._iterate_records())}
Maximum sequence length: {max(len(seq) for _, seq in sequences_step._iterate_records())}
Number of genomes containing IRS by NCBI annotation  : {len(ncbi_with)}
Number of genomes containing IRS by GeSeq annotation : {len(ge_seq_with)}
Number of genomes without IRS                        : {len(data) - len(set(ge_seq_with + ncbi_with))}

Range of dates:
  First date     : {min(d.first_date for d in data.values())} - {max(d.first_date for d in data.values())}
  Published date : {min(d.created_date for d in data.values())} - {max(d.created_date for d in data.values())}

Range of number of genes:
  Annotated        : {min_genes['annotated']} - {max_genes['annotated']}
  Disjunct         : {min_genes['disjunct']} - {max_genes['disjunct']}
  Name/strand      : {min_genes['name_strand']} - {max_genes['name_strand']}
  Name             : {min_genes['names']} - {max_genes['names']}
  Without location : {min_genes['without_location']} - {max_genes['without_location']}
  Without name     : {min_genes['without_name']} - {max_genes['without_name']}

Range of part lengths:
  LSC : {min(l[0] for l in p_lengths)} - {max(l[0] for l in p_lengths)}
  IRs : {min(l[1] for l in p_lengths)} - {max(l[1] for l in p_lengths)}
  SSC : {min(l[2] for l in p_lengths)} - {max(l[2] for l in p_lengths)}

Corrections by parts:
  Wrong orientation : {sum(1 for d in self.seq_descs.values() if d.part_orientation)}{parts_orient_count}
  With offset       : {sum(1 for d in self.seq_descs.values() if offset_off(d.part_offset))}
  Fixed             : {sum(1 for d in self.seq_descs.values() if d.part_orientation or offset_off(d.part_offset))}

{self._find_species_stats()}
"""

        print(summary)
        write_str_in_file(self.step.step_file('summary.txt'), summary)

        #
        if errors := self._find_errors():
            r_file = self.step.step_file('errors.txt')
            write_str_in_file(r_file, '\n'.join(errors) + '\n')
            print(f'Problems found in sequences! Check file {r_file}.')

    def _find_species_stats(self):
        ncbi_tax = get_ncbi_taxonomy()
        ident_2_taxid = dict(self.table_data.iterate_column('tax_id'))
        taxid_2_ident = dict((v, k) for k, v in ident_2_taxid.items())
        species, without_sp = ncbi_tax.group_taxids_in_species(set(ident_2_taxid.values()))

        # Species with more genomes
        ret = []
        for (sp_taxid, sp_name), data in species.items():
            if len(data) > 1:
                ret.append(f'  Genomes for species {sp_name} ({sp_taxid}):')
                ret.extend(f'   - {taxid_2_ident[taxid]} {name} ({taxid}) of rank {rank}' for taxid, name, rank in data)
        if ret:
            ret = ['Species with more accessions:'] + ret
        else:
            ret.append('No species found with more genomes!')

        # Without species
        if without_sp:
            ret.extend(['', "No 'base' species in NCBI taxa database found for:"])
            ret.extend(f' - {taxid_2_ident[taxid]} {name} ({taxid}) of rank {rank}' for taxid, name, rank in without_sp)
        return '\n'.join(ret)

    def _find_errors(self):
        offset_off = lambda o: (abs(o or 0) > DEFAULT_KEEP_OFFSET)
        data = self.seq_descs
        errors = []
        l_start = '\n - '
        if f_ps := [(s, d.irs_took_from) for s, d in data.items() if d.irs_took_from]:
            errors.append(f'Found partitions for sequences:{l_start}{l_start.join(sorted(f"{s} ({m})" for s, m in f_ps))}')
        if without_parts := [seq_ident for seq_ident, d in data.items() if not d._partition]:
            errors.append(f"No partitions for sequences:{l_start}{l_start.join(sorted(without_parts))}")
        if wrong_oriented_parts := [seq_ident for seq_ident, d in data.items() if d.part_orientation]:
            errors.append(f"Partitions with wrong orientation in sequences:{l_start}{l_start.join(sorted(wrong_oriented_parts))}")
        if with_offset := [seq_ident for seq_ident, d in data.items() if offset_off(d.part_offset)]:
            errors.append(f"Sequences with offset:{l_start}{l_start.join(sorted(with_offset))}")
        # if with_trnH_offset := [seq_ident for seq_ident, d in data.items() if offset_off(d.part_trnH_GUG)]:
        #     errors.append(f"Sequneces with trnH-GUG offset:{l_start}{l_start.join(sorted(with_trnH_offset))}")
        return errors


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
