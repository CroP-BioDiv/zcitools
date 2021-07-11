import os
import tempfile
import subprocess
from Bio import SeqIO


def small_d(seq_rec, working_dir=None):
    if not working_dir:
        working_dir = tempfile.gettempdir()
    fasta_filename = os.path.join(working_dir, f'{seq_rec.name}.fa')
    irscan_exe = os.environ.get('IRSCAN_EXE', 'irscan')

    for offset in (0, 500, 5000, 25000, 30000):
        s_rec = (seq_rec[-offset:] + seq_rec) if offset else seq_rec
        SeqIO.write([s_rec], fasta_filename, 'fasta')
        try:
            result = subprocess.run([irscan_exe, '-f', fasta_filename], check=True, stdout=subprocess.PIPE)
        except subprocess.CalledProcessError:
            print(f'\nWarning: sequence {seq_rec.name} has IR on the start, with offset {offset}!')
            continue
        #
        os.remove(fasta_filename)
        seq_length = len(seq_rec.seq)
        irs = tuple(int(x) - offset for x in result.stdout.decode('utf-8').split(';')[:4])
        ira, irb = _ir(seq_length, *irs[:2]), _ir(seq_length, *irs[2:])

        # Check IR order
        if (irb[0] - ira[1]) % seq_length < (ira[0] - irb[1]) % seq_length:
            return ira, irb
        return irb, ira


def _ir(seq_length, ira, irb):
    # Note: irscan output is 1-based, and range means [start..end]
    return ((ira - 1) % seq_length), (irb if irb > 0 else (irb % seq_length))


def small_d_on_file(seq_filename):
    _ext_2_bio_io_type = dict(
        gb='genbank', gbff='genbank',
        fa='fasta',  fas='fasta',
        fastq='fastq',
    )

    base_filename, file_extension = os.path.splitext(seq_filename)
    in_format = _ext_2_bio_io_type[file_extension[1:]]
    return small_d(SeqIO.read(seq_filename, in_format))


if __name__ == '__main__':
    import sys
    print(small_d_on_file(sys.argv[1]))
