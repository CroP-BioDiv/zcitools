import os.path
# import multiprocessing
# from concurrent.futures import ThreadPoolExecutor
from step_project.common.table.steps import TableStep
from common_utils.file_utils import ensure_directory, write_str_in_file, write_fasta
from common_utils.properties_db import PropertiesDB
from common_utils.cache import cache
from .analyse_step_1 import SequenceDesc
from .analyse_step_calc import run_align_cmd, find_missing_partitions
from .utils import find_chloroplast_partition
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
    step.to_excel('analyse_chloroplast.xls')  # Store data for a check
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
            ('num_genes', 'Genes', 'int'),
            ('num_cds', 'CDS', 'int'),
            ('irs_took_from', 'IRS took', 'seq_ident'),
            ('part_starts', 'Part starts', 'str'),
            ('part_lengths', 'Part lengths', 'str'),
            ('part_num_genes', 'Part genes', 'str'),
            ('part_orientation', 'Orientation', 'str'),
            ('trnH_GUG', 'trnH-GUG', 'str'),
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
        sequences_step = next(iter(data.values())).sequences_step
        sp_stats = self._find_species_stats()
        num_ge_seq_irs = sum(int(bool(seq_data._partition)) for seq_data in data.values())
        num_ncbi_irs = sum(int(bool(find_chloroplast_partition(seq))) for _, seq in sequences_step._iterate_records())
        summary = f"""Statistics:

Minimum sequence length: {min(len(seq) for _, seq in sequences_step._iterate_records())}
Maximum sequence length: {max(len(seq) for _, seq in sequences_step._iterate_records())}
Number of genomes containing IRS by NCBI annotation  : {num_ncbi_irs}
Number of genomes containing IRS by GeSeq annotation : {num_ge_seq_irs}
"""
        if sp_stats := self._find_species_stats():
            summary += '\nSpecies with more accessions:\n\n' + '\n'.join(sp_stats) + '\n'

        if corrections := self._find_corrections_stat():
            summary += '\nCorrections:\n\n' + '\n'.join(corrections) + '\n'
        else:
            summary += '\nNo corrections were done!'

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
                ret.append(f'Species {sp_name} ({sp_taxid}) has more genomes:')
                ret.extend(f' - {taxid_2_ident[taxid]} {name} ({taxid}) of rank {rank}' for taxid, name, rank in data)
        if not ret:
            ret.append('No species found with more genomes!')

        # Wihtout species
        if without_sp:
            ret.append("No 'base' species in NCBI taxa database found for:")
            ret.extend(f' - {taxid_2_ident[taxid]} {name} ({taxid}) of rank {rank}' for taxid, name, rank in without_sp)
        return ret

    def _find_corrections_stat(self):
        corrs = []
        if num := sum(1 for d in self.seq_descs.values() if d.part_orientation):
            corrs.append(f'Found partitions: {num}')
        if num := sum(1 for d in self.seq_descs.values() if not d._partition):
            corrs.append(f'Without partitions: {num}')
        if num := sum(1 for d in self.seq_descs.values() if d.part_orientation):
            corrs.append(f'Wrong partition orientation: {num}')
        if num := sum(1 for d in self.seq_descs.values() if abs(d.part_offset or 0) > 50):
            corrs.append(f'Wrong LSC offset: {num}')
        if num := sum(1 for d in self.seq_descs.values() if abs(d.part_trnH_GUG or 0) > 50):
            corrs.append(f'Wrong trnH-GUG offset: {num}')
        #
        if num := sum(1 for d in self.seq_descs.values()
                      if d.part_orientation or
                      not d._partition or
                      d.part_orientation or
                      abs(d.part_offset or 0) > 50 or
                      abs(d.part_trnH_GUG or 0) > 50):
            corrs.append(f'Fixed: {num}')
        return corrs

    def _find_errors(self):
        data = self.seq_descs
        errors = []
        l_start = '\n - '
        if f_ps := [(s, d.irs_took_from) for s, d in data.items() if d.irs_took_from]:
            errors.append(f'Found partitions for sequences:{l_start}{l_start.join(sorted(f"{s} ({m})" for s, m in f_ps))}')
        if without_parts := [seq_ident for seq_ident, d in data.items() if not d._partition]:
            errors.append(f"No partitions for sequences:{l_start}{l_start.join(sorted(without_parts))}")
        if wrong_oriented_parts := [seq_ident for seq_ident, d in data.items() if d.part_orientation]:
            errors.append(f"Partitions with wrong orientation in sequences:{l_start}{l_start.join(sorted(wrong_oriented_parts))}")
        if with_offset := [seq_ident for seq_ident, d in data.items() if abs(d.part_offset or 0) > 50]:
            errors.append(f"Sequences with offset:{l_start}{l_start.join(sorted(with_offset))}")
        if with_trnH_offset := [seq_ident for seq_ident, d in data.items() if abs(d.part_trnH_GUG or 0) > 50]:
            errors.append(f"Sequneces with trnH-GUG offset:{l_start}{l_start.join(sorted(with_trnH_offset))}")
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
