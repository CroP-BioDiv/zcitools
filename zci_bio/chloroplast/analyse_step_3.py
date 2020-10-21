import os
import multiprocessing
from concurrent.futures import ThreadPoolExecutor
from common_utils.file_utils import ensure_directory, write_fasta
from .utils import create_chloroplast_partition
from ..utils.mummer import MummerDelta
from ..utils.ncbi_taxonomy import ncbi_taxonomy


def run_align_cmd(seq_fasta, qry_fasta, out_prefix):
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


def find_missing_partitions(step, data, table_step):
    with_irs = set(d.taxid for d in data.values() if d._parts)
    if len(with_irs) == len(data):  # All in
        return
    if not with_irs:
        print('Warning: all genomes miss IR parts!!!')
        return

    ncbi_2_max_taxid = table_step.mapping_between_columns('ncbi_ident', 'max_taxid')
    taxid_2_ncbi = table_step.mapping_between_columns('tax_id', 'ncbi_ident')

    # Run alignment of related species IR ends onto sequences without IRs
    ncbi_tax = ncbi_taxonomy()
    seq_2_result_object = dict()  # dict seq_ident -> tuple (ira interval, irb interval, matche seq_ident)
    with ThreadPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
        for seq_ident, seq_data in data.items():
            if seq_data._parts:
                continue

            close_taxids = ncbi_tax.find_close_taxids(seq_data.taxid, ncbi_2_max_taxid[seq_ident], with_irs)
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
        d._took_parts = create_chloroplast_partition(d.length, ira, irb, in_interval=True)
        d.irs_took_from = transfer_from

    return seq_2_result_object


def _run_align(step, seq_ident, seq_data, close_data, seq_2_result_object, taxid_2_ncbi, match_length=100):
    # Store fasta
    f_dir = step.step_file('find_irs', seq_ident)
    ensure_directory(f_dir)
    seq_fasta = os.path.join(f_dir, f"{seq_ident}.fa")
    write_fasta(seq_fasta, [(seq_ident, seq_data._seq.seq)])

    all_aligns = []
    for d in close_data:
        ira = d._parts.get_part_by_name('ira')
        rec = ira.extract(d._seq)
        qry_fasta = os.path.join(f_dir, f"qry_{d.seq_ident}.fa")
        write_fasta(qry_fasta, [('end1', rec.seq[:match_length]), ('end2', rec.seq[-match_length:])])
        align = run_align_cmd(seq_fasta, qry_fasta, f"res_{d.seq_ident}")
        if 2 == len(align.aligns(seq_ident, 'end1')) == len(align.aligns(seq_ident, 'end2')):  # All sides
            ira_1, irb_2 = align.aligns(seq_ident, 'end1')
            ira_2, irb_1 = align.aligns(seq_ident, 'end2')
            seq_2_result_object[seq_ident] = \
                ((ira_1.sequence_interval[0], ira_2.sequence_interval[1]),
                 (irb_1.sequence_interval[0], irb_2.sequence_interval[1]),
                 d.seq_ident)
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
                     d.seq_ident)
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
                     d.seq_ident)
                return
