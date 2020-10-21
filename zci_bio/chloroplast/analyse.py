import os.path
from datetime import datetime
# import multiprocessing
# from concurrent.futures import ThreadPoolExecutor
from step_project.common.table.steps import TableStep
from common_utils.file_utils import ensure_directory, write_fasta
from common_utils.properties_db import PropertiesDB
from .analyse_step_1 import SequenceDesc
from .analyse_step_2 import evaluate_credibility
from .analyse_step_3 import run_align_cmd, find_missing_partitions
from ..utils.features import Feature
from ..utils.ncbi_taxonomy import ncbi_taxonomy


# ---------------------------------------------------------
# Analyse genomes
# Steps:
#  1. Collect data from steps (input annotations, previous table)
#  2. Evaluate credibility of collected data
#  3. Find missing or more credible data
#  4. Set table from collected and evaluated data
# ---------------------------------------------------------
def analyse_genomes(step_data, annotations_step):
    project = annotations_step.project
    step = TableStep(project, step_data, remove_data=True)
    annotations_step.propagate_step_name_prefix(step)
    properties_db = PropertiesDB()

    # 1. Collect data from steps
    table_step = project.find_previous_step_of_type(annotations_step, 'table')
    sequence_step = project.find_previous_step_of_type(annotations_step, 'sequences')
    table_data = table_step.index_on_table('ncbi_ident')

    data = dict((seq_ident, SequenceDesc(seq_ident, seq, table_data, sequence_step, properties_db))
                for seq_ident, seq in annotations_step._iterate_records())

    #  2. Evaluate credibility of collected data
    evaluate_credibility(data, annotations_step)

    #  3. Find missing or more credible data
    find_missing_partitions(step, data, table_step)

    # Set partition data
    for seq_ident, d in data.items():
        d.set_parts_data()

    #  4. Set table from collected and evaluated data
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
        ('part_starts', 'Part starts', 'str'),
        ('part_lengths', 'Part lengths', 'str'),
        ('part_num_genes', 'Part genes', 'str'),
        ('irs_took_from', 'IRS took', 'seq_ident'),
        ('took_part_starts', 'Took part starts', 'str'),
        ('offset', 'Offset', 'int'),
        ('trnH_GUG', 'trnH-GUG', 'int'),
        ('part_orientation', 'Orientation', 'str'),
        ('artcle_title', 'Article', 'str'),
        ('journal', 'Journal', 'str'),
        ('pubmed_id', 'PubMed', 'int'),
        ('assembly_method', 'Assembly Method', 'str'),
        ('sequencing_technology', 'Sequencing Technology', 'str'),
        ('sra_count', 'SRA count', 'int'),
    ]
    step.set_table_data(
        [[getattr(d, c) for c, _, _ in columns] for seq_ident, d in sorted(data.items())],
        [(n, t) for _, n, t in columns])
    step.save()

    #
    errors = []
    l_start = '\n - '
    if f_ps := [(s, d.irs_took_from) for s, d in data.items() if d.irs_took_from]:
        errors.append(f'Found partitions for sequences:{l_start}{l_start.join(sorted(f"{s} ({m})" for s, m in f_ps))}')
    if without_parts := [seq_ident for seq_ident, d in data.items() if not d._parts]:
        errors.append(f"No partitions for sequences:{l_start}{l_start.join(sorted(without_parts))}")
    if wrong_oriented_parts := [seq_ident for seq_ident, d in data.items() if d.part_orientation]:
        errors.append(f"Partitions with wrong orientation in sequences:{l_start}{l_start.join(sorted(wrong_oriented_parts))}")
    if with_offset := [seq_ident for seq_ident, d in data.items() if abs(d.offset) > 50]:
        errors.append(f"Sequences with offset:{l_start}{l_start.join(sorted(with_offset))}")
    if with_trnH_offset := [seq_ident for seq_ident, d in data.items() if abs(d.trnH_GUG) > 50]:
        errors.append(f"Sequneces with trnH-GUG offset:{l_start}{l_start.join(sorted(with_trnH_offset))}")
    if errors:
        with open((r_file := step.step_file('README.txt')), 'w') as _out:
            _out.write('\n'.join(errors))
            _out.write('\n')
        print(f'Problems found in sequences! Check file {r_file}.')
    step.to_excel('analyse.xls')  # Test

    return step


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
        ncbi_tax = ncbi_taxonomy()
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

    print(len(parts), len(ns))
    # Extract query data
    write_fasta(query_filename, [(f"{f.real_start}_{f.real_end}", f.extract(seq).seq) for f in parts])
    align = run_align_cmd(close_filename, query_filename, 'result')
    # ToDo: ...
