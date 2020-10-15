import os.path
import multiprocessing
from concurrent.futures import ThreadPoolExecutor
from step_project.common.table.steps import TableStep
from common_utils.file_utils import ensure_directory, write_fasta
from .utils import find_chloroplast_partition, create_chloroplast_partition
from ..utils.features import Feature
from ..utils.mummer import MummerDelta
from ..utils.ncbi_taxonomy import ncbi_taxonomy


def _run_align_cmd(seq_fasta, qry_fasta, out_prefix):
    # Blast version
    # cmd = f"blastn -subject {seq_fasta} -query {qry_fasta} " + \
    #     f"-perc_identity 40 -num_alignments 2 -max_hsps 2 -outfmt 5 > {out_prefix}.xml"
    # ...

    # Mummer version
    out_prefix = os.path.join(os.path.dirname(seq_fasta), out_prefix)
    cmd = f"nucmer -p {out_prefix} {seq_fasta} {qry_fasta}"
    print(f"Command: {cmd}")
    os.system(cmd)

    # Read output file
    return MummerDelta(f'{out_prefix}.delta')


# ---------------------------------------------------------
# Analyse genomes
# ---------------------------------------------------------
def analyse_genomes(step_data, annotations_step):
    project = annotations_step.project
    step = TableStep(project, step_data, remove_data=True)

    #
    table_step = project.find_previous_step_of_type(annotations_step, 'table')
    ncbi_2_taxid = table_step.mapping_between_columns('ncbi_ident', 'tax_id')
    taxid_2_ncbi = dict((v, k) for k, v in ncbi_2_taxid.items())
    ncbi_2_title = table_step.mapping_between_columns('ncbi_ident', 'title')
    ncbi_2_max_taxid = table_step.mapping_between_columns('ncbi_ident', 'max_taxid')

    data = dict((seq_ident, dict(
            _seq=seq,
            _parts=find_chloroplast_partition(seq),
            seq_ident=seq_ident,
            length=len(seq),
            genes=sum(1 for f in seq.features if f.type == 'gene' and f.location),
            cds=sum(1 for f in seq.features if f.type == 'CDS' and f.location),
            taxid=ncbi_2_taxid[seq_ident],
            title=ncbi_2_title[seq_ident],
            irs_transfered_from=None,
        )) for seq_ident, seq in annotations_step._iterate_records())

    # Find missing IR partitions
    find_missing_partitions(step, data, ncbi_2_max_taxid, taxid_2_ncbi)

    # Set partition data
    cs = ('lsc', 'ssc', 'ira')
    for seq_ident, d in data.items():
        if parts := d['_parts']:
            # Part lengths
            for part in cs:
                d[part] = len(parts.get_part_by_name(part))

            # SSC ends
            ssc_start, ssc_end = parts.get_part_by_name('ssc').ends()
            d['ssc_ends'] = f'{ssc_start}-{ssc_end}'

            # Number of genes in parts
            seq = d['_seq']
            length = len(seq)
            part_genes = parts.put_features_in_parts(
                [Feature(length, feature=f) for f in seq.features if f.type == 'gene' and f.location])
            for p in ('lsc', 'ssc', 'ira', 'irb'):
                d[f'{p}_genes'] = len(part_genes[p])

            # Offset
            d['offset'] = parts.get_part_by_name('lsc').real_start

            # Part orientation
            orient = chloroplast_parts_orientation(d['_seq'], parts)
            d['part_orientation'] = ''.join('P' if orient[p] else 'N' for p in ('lsc', 'ira', 'ssc'))

        else:
            for x in cs + ('ssc_ends', 'lsc_genes', 'ssc_genes', 'ira_genes', 'irb_genes',
                           'offset', 'part_orientation'):
                d[x] = None

    # Extract data, and set table properties
    columns = [
        # tuples (dict's attribute, column name, column type)
        ('seq_ident', 'AccesionNumber', 'seq_ident'),
        ('title', 'Title', 'str'),
        ('length', 'Length', 'int'),
        ('genes', 'Genes', 'int'),
        ('irs_transfered_from', 'IRS took', 'seq_ident'),
        ('cds', 'CDS', 'int'),
        ('lsc', 'LSC', 'int'),
        ('ssc', 'SSC', 'int'),
        ('ira', 'IR', 'int'),
        ('ssc_ends', 'SSC ends', 'str'),
        ('lsc_genes', 'LSC genes', 'int'),
        ('ssc_genes', 'SSC genes', 'int'),
        ('ira_genes', 'IRA genes', 'int'),
        ('irb_genes', 'IRB genes', 'int'),
        ('offset', 'Offset', 'int'),
        ('part_orientation', 'Orientation', 'str'),
    ]
    step.set_table_data(
        [[d[c] for c, _, _ in columns] for seq_ident, d in sorted(data.items())],
        [(n, t) for _, n, t in columns])
    step.save()
    step.to_excel('analyse.xls')  # Test
    return step


def find_missing_partitions(step, data, ncbi_2_max_taxid, taxid_2_ncbi):
    with_irs = set(d['taxid'] for d in data.values() if d['_parts'])
    if len(with_irs) == len(data):  # All in
        return
    if not with_irs:
        print('Warning: all genomes miss IR parts!!!')
        return

    # Run alignment of related species IR ends onto sequences without IRs
    ncbi_tax = ncbi_taxonomy()
    seq_2_result_object = dict()
    with ThreadPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
        for seq_ident, seq_data in data.items():
            if seq_data['_parts']:
                continue

            close_taxids = ncbi_tax.find_close_taxids(seq_data['taxid'], ncbi_2_max_taxid[seq_ident], with_irs)
            if not close_taxids:
                print(f"Warning: sequence {seq_ident} doesn't have close relative with IR partition!")
                continue

            # _run_align(step, seq_ident, seq_data, [data[taxid_2_ncbi[t]] for t in close_taxids],
            #            seq_2_result_object, taxid_2_ncbi)
            executor.submit(_run_align, step, seq_ident, seq_data, [data[taxid_2_ncbi[t]] for t in close_taxids],
                            seq_2_result_object, taxid_2_ncbi)

    #
    for seq_ident, (ira, irb, transfer_from) in seq_2_result_object.items():
        d = data[seq_ident]
        d['_parts'] = create_chloroplast_partition(len(d['_seq']), ira, irb, in_interval=True)
        d['irs_transfered_from'] = transfer_from


def _run_align(step, seq_ident, seq_data, close_data, seq_2_result_object, taxid_2_ncbi, match_length=100):
    # Store fasta
    f_dir = step.step_file('find_irs', seq_ident)
    ensure_directory(f_dir)
    seq_fasta = os.path.join(f_dir, f"{seq_ident}.fa")
    write_fasta(seq_fasta, [(seq_ident, seq_data['_seq'].seq)])

    all_aligns = []
    for d in close_data:
        ira = d['_parts'].get_part_by_name('ira')
        rec = ira.extract(d['_seq'])
        qry_fasta = os.path.join(f_dir, f"qry_{d['seq_ident']}.fa")
        write_fasta(qry_fasta, [('end1', rec.seq[:match_length]), ('end2', rec.seq[-match_length:])])
        align = _run_align_cmd(seq_fasta, qry_fasta, f"res_{d['seq_ident']}")
        if 2 == len(align.aligns(seq_ident, 'end1')) == len(align.aligns(seq_ident, 'end2')):  # All sides
            ira_1, irb_2 = align.aligns(seq_ident, 'end1')
            ira_2, irb_1 = align.aligns(seq_ident, 'end2')
            seq_2_result_object[seq_ident] = \
                ((ira_1.sequence_interval[0], ira_2.sequence_interval[1]),
                 (irb_1.sequence_interval[0], irb_2.sequence_interval[1]),
                 d['seq_ident'])
            break
        all_aligns.append(align)
    else:
        # No alignment matched all ends.
        # Try partial (2+1)
        for align in all_aligns:
            if len(align.aligns(seq_ident, 'end1')) == 2 and len(align.aligns(seq_ident, 'end2')) == 1:
                ira_1, irb_2 = align.aligns(seq_ident, 'end1')
                ir = align.aligns(seq_ident, 'end2')[0]
                if ir.positive:  # In positive direction (IRA matched)
                    x = ir.sequence_interval[1]
                    y = irb_2.sequence_interval[1] - (x - ira_1.sequence_interval[0])
                else:
                    y = ir.sequence_interval[0]
                    x = ira_1.sequence_interval[0] + (irb_2.sequence_interval[1] - y)
                seq_2_result_object[seq_ident] = \
                    ((ira_1.sequence_interval[0], x),
                     (y, irb_2.sequence_interval[1]),
                     d['seq_ident'])
                return
        # Try partial (1+2)
        for align in all_aligns:
            if len(align.aligns(seq_ident, 'end1')) == 1 and len(align.aligns(seq_ident, 'end2')) == 2:
                ir = align.aligns(seq_ident, 'end1')[0]
                ira_2, irb_1 = align.aligns(seq_ident, 'end2')
                if ir.positive:  # In positive direction (IRA matched)
                    x = ir.sequence_interval[0]
                    y = irb_1.sequence_interval[0] + (ira_2.sequence_interval[1] - x)
                else:
                    y = ir.sequence_interval[1]
                    x = ira_2.sequence_interval[1] - (y - irb_1.sequence_interval[0])
                seq_2_result_object[seq_ident] = \
                    ((x, ira_2.sequence_interval[1]),
                     (irb_1.sequence_interval[0], y),
                     d['seq_ident'])
                return


def chloroplast_parts_orientation(seq_rec, partition):
    # Check chloroplast sequence part orientation.
    # Default orientation is same as one uses in Fast-Plast. Check:
    #  - source file orientate_plastome_v.2.0.pl
    #    (https://github.com/mrmckain/Fast-Plast/blob/master/bin/orientate_plastome_v.2.0.pl)
    #  - explanation https://github.com/mrmckain/Fast-Plast/issues/22
    # Consitent with Wikipedia image:
    #  - https://en.wikipedia.org/wiki/Chloroplast_DNA#/media/File:Plastomap_of_Arabidopsis_thaliana.svg

    l_seq = len(seq_rec)
    in_parts = partition.put_features_in_parts(
        Feature(l_seq, feature=f) for f in seq_rec.features if f.type == 'gene' and f.location)

    lsc_count = sum(f.feature.strand if any(x in f.name for x in ('rpl', 'rps')) else 0
                    for f in in_parts.get('lsc', []))
    ssc_count = sum(f.feature.strand for f in in_parts.get('ssc', []))
    ira_count = sum(f.feature.strand if 'rrn' in f.name else 0 for f in in_parts.get('ira', []))

    return dict(lsc=(lsc_count <= 0),
                ssc=(ssc_count <= 0),
                ira=(ira_count >= 0))


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
    align = _run_align_cmd(close_filename, query_filename, 'result')
    # ToDo: ...
