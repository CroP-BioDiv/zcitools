import os
from concurrent.futures import ThreadPoolExecutor
import multiprocessing
from common_utils.exceptions import ZCItoolsValueError
from common_utils.file_utils import ensure_directory, write_fasta, write_yaml, run_module_script, set_run_instructions
from common_utils.value_data_types import rows_2_excel
from zci_bio.annotations.steps import AnnotationsStep
from ..utils.import_methods import import_bio_seq_io, import_bio_align_io
from . import run_mafft_irs

_instructions = ''


# Mummer methods
def _run_single(mummer_exe, n, input_filename, output_filename):
    cmd = f"{mummer_exe} -n {n} {input_filename} > {output_filename}"
    print(f"Command: {cmd}")
    os.system(cmd)


def _read_mummer_repeat(result_filename):
    max_r = None
    with open(result_filename, 'r') as output:
        read = False
        for line in output:
            fields = line.strip().split()
            if read:
                if fields[1][-1] == 'r':  # Inverted
                    length = int(fields[2])
                    if max_r is None or length > max_r[0]:
                        # Mummer indexed from 1?
                        start_1 = int(fields[0]) - 1
                        start_2 = int(fields[1][:-1]) - 1
                        # start_2 = start_2 - length + 1
                        max_r = [length, start_1, start_2]  # Return list, because of yaml serialization
            elif fields[0] == 'Start1':
                read = True
    return max_r


# Mafft methods
def _extract_subseq_plus(seq, from_ind, length):
    if from_ind + length > len(seq):
        return seq[from_ind:] + seq[:length - (len(seq) - from_ind)]
    return seq[from_ind: from_ind + length]


def _extract_subseq_minus(seq, from_ind, length):
    if from_ind - length < 0:
        return seq[from_ind::-1] + seq[:-(length - from_ind - 1):-1]
    return seq[from_ind:from_ind - length:-1]


#
def create_irs_data(step_data, input_step, params):
    # Creates Annotations step from input sequences/annotations
    # Steps subdirectory 'run_dir' contains input and output calculation files
    SeqIO = import_bio_seq_io()
    seq_idents = input_step.all_sequences()

    step = input_step.project.new_step(AnnotationsStep, step_data)
    step.set_sequences(seq_idents)
    # seq_ident -> mummer data ([length, start_1, start_2])
    mummer_results = step.get_type_description_elem('mummer_results', default=dict())
    #
    ensure_directory(step.step_file('run_dir'))
    calc_mummer = []  # tuples (seq_ident, fasta file, mummer output file)

    # Mummer
    for seq_ident in sorted(seq_idents - set(mummer_results)):
        fa_file = step.step_file('run_dir', f'{seq_ident}.fa')
        mummer_res_file = step.step_file('run_dir', f'{seq_ident}.out')
        if not os.path.isfile(fa_file):
            seq_rec = input_step.get_sequence_record(seq_ident)
            SeqIO.write([seq_rec], fa_file, 'fasta')
            calc_mummer.append((seq_ident, fa_file, mummer_res_file))
        elif not os.path.isfile(mummer_res_file):
            calc_mummer.append((seq_ident, fa_file, mummer_res_file))

    # Run mummer
    if calc_mummer:
        mummer_exe = 'repeat-match'  # ToDo:
        n = 3000
        threads = multiprocessing.cpu_count()
        with ThreadPoolExecutor(max_workers=threads) as executor:
            for seq_ident, fa_file, mummer_res_file in calc_mummer:
                executor.submit(_run_single, mummer_exe, n, fa_file, mummer_res_file)

        for seq_ident, _, mummer_res_file in calc_mummer:
            rep = _read_mummer_repeat(mummer_res_file)
            if not rep:
                raise ZCItoolsValueError(f'No repeat for sequence {seq_ident}!')
            mummer_results[seq_ident] = rep

    # Find sequences extend with alignment
    files_to_zip = []
    calc_mafft = []
    for seq_ident in sorted(seq_idents):
        length, s1, s2 = mummer_results[seq_ident]
        if length >= 23000:
            continue

        if step.is_file('run_dir', f'{seq_ident}_right_align.fa') and \
           step.is_file('run_dir', f'{seq_ident}_right_align.fa'):
            continue

        #
        calc_mafft.append(seq_ident)
        _seq = input_step.get_sequence_record(seq_ident).seq
        seq = str(_seq)
        comp_seq = str(_seq.complement())
        missing = 26000 - length

        # Right side
        p1 = _extract_subseq_plus(seq, s1 + length, missing)
        p2 = _extract_subseq_minus(comp_seq, s2 - length, missing)
        assert len(p1) == len(p2), (length, s1, s2, missing, (len(p1), len(p2)))
        files_to_zip.append(step.step_file('run_dir', f'{seq_ident}_right.fa'))
        write_fasta(files_to_zip[-1], [('p1', p1), ('p2', p2)])

        # Left side
        p1 = _extract_subseq_minus(comp_seq, s1 - 1, missing)
        p2 = _extract_subseq_plus(seq, s2 + 1, missing)
        assert len(p1) == len(p2), (length, s1, s2, missing, (len(p1), len(p2)))
        files_to_zip.append(step.step_file('run_dir', f'{seq_ident}_left.fa'))
        write_fasta(files_to_zip[-1], [('p1', p1), ('p2', p2)])

    # Mafft
    if calc_mafft:
        finish_f = step.step_file('finish.yml')
        write_yaml(dict(calc_seq_idents=calc_mafft), finish_f)

        run = True  # ToDo: ...
        step.save(additional_data=dict(mummer_results=mummer_results), completed=False)
        if run:
            run_module_script(run_mafft_irs, step)
            finish_irs_data(step)
        else:
            files_to_zip.append(finish_f)
            set_run_instructions(run_mafft_irs, step, files_to_zip, _instructions)
    #
    elif params.force_parse:
        finish_irs_data(step)

    return step


def finish_irs_data(step_object):
    mummer_results = step_object.get_type_description_elem('mummer_results')
    assert mummer_results

    rows = []
    for seq_ident in sorted(step_object.all_sequences()):
        length, s1, s2 = mummer_results[seq_ident]
        # print(length, s1, s2)
        if length >= 23000:
            rows.append([seq_ident, 'mummer', length, s1 + length, s2 - length])
            continue
        #
        d_left = _mafft_length(step_object.step_file('run_dir', f'{seq_ident}_left_align.fa'))
        s1 -= d_left
        d_right = _mafft_length(step_object.step_file('run_dir', f'{seq_ident}_right_align.fa'))
        s2 += d_right
        length += d_left + d_right
        rows.append([seq_ident, 'MAFFT', length, s1 + length, s2 - length])
        print(seq_ident, length, s1, s2)

    # Create table
    output_file = f'chloroplast_blast_ssc.xlsx'
    columns = ['Seq', 'Method', 'Length', 'SSC start', 'SSC end']
    rows_2_excel(output_file, columns, rows)


def _mafft_length(align_file):
    a_m = import_bio_align_io().read(align_file, 'fasta')
    return min(map(_mafft_l, a_m))


def _mafft_l(seq_rec):
    # States: align, gap, check_gap
    state = 'align'
    last_gap = 0
    last_check_gap = 0
    num_chars = 0
    for c in str(seq_rec.seq):
        if c == '-':
            if state == 'align':
                state = 'gap'
                last_gap = 1
            elif state == 'gap':
                last_gap += 1
            else:  # check_gap
                # Next gap starts before check length passed
                return num_chars
        else:
            if state == 'align':
                num_chars += 1
            elif state == 'gap':
                state = 'check_gap'
                last_check_gap = 1
                last_gap *= 10  # How much next alignment part should be longer than last gap
            else:  # check_gap
                last_check_gap += 1
                if last_check_gap >= last_gap:
                    state = 'align'
                    num_chars += last_check_gap
    return num_chars
