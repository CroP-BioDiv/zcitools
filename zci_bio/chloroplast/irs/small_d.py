#!/usr/bin/env python3

import os
import re
import tempfile
import subprocess
from Bio import SeqIO, Seq

"""
Wrapper around irscan program. irscan program is used by OGDRAW to locate IRs.

'For quick and precise detection of inverted repeat regions in organellar genomes,
a small D program (http://www.digitalmars.com/d/index.html) was developed and included in the package.'
From paper:
(2007) OrganellarGenomeDRAW (OGDRAW): a tool for the easy generation of high-quality custom graphical
maps of plastid and mitochondrial genomes.
https://link.springer.com/article/10.1007/s00294-007-0161-y

'Small D' executables and source can be found on OGDRAW download page:
 - https://chlorobox.mpimp-golm.mpg.de/OGDraw-Downloads.html
 - File is GeneMap-1.1.1.tar.gz.
Subproject irscan contains directories:
 - bin : Win, Linux and Mac executables.
 - src : D code.

Note: it is not possible to compile code with current D version.

Note: irscan crashes on some 'proper' input.
This script handle crashes and do workarounds to get meaningfull result.
These workarounds can be turned off with script switches.

I encountered two types of crashes, which stderr outputs are:
Error: ArrayBoundsError irscan(161)
 - IR starts at sequence start. Method screen_from_to_pos() moves pos into negative
     - Workaround is to run irscan on sequence prepended by sequence end. Like: seq[-100:] + seq

Error: ArrayBoundsError irscan(354)
 - IRs were not located. Printing of results crashes.
    - No workaround since there are no IRs. At least not in irscan definition of IRs.

 - Can be caused if sequence contains character different from ATCG.
   irscan.d, return in line 82 is dangerous.
"""


def small_d(seq_rec, working_dir=None, leave_tmp_file=False, no_prepend_workaround=False, no_dna_fix=False):
    if not working_dir:
        working_dir = tempfile.gettempdir()
    fasta_filename = os.path.join(working_dir, f'{seq_rec.name}.fa')
    irscan_exe = os.environ.get('IRSCAN_EXE', 'irscan')
    offsets = (0,) if no_prepend_workaround else (0, 500, 5000, 25000, 30000)
    res_irs = None

    if not no_dna_fix and any(c not in 'ATCG' for c in str(seq_rec.seq)):
        dna = re.sub(r'[^ATCG]', 'A', str(seq_rec.seq))
        seq_rec.seq = Seq.Seq(dna)

    for offset in offsets:
        s_rec = (seq_rec[-offset:] + seq_rec) if offset else seq_rec
        SeqIO.write([s_rec], fasta_filename, 'fasta')
        with tempfile.TemporaryFile() as stderr:
            try:
                result = subprocess.run([irscan_exe, '-f', fasta_filename],
                                        check=True,
                                        stdout=subprocess.PIPE, stderr=stderr)
            except subprocess.CalledProcessError:
                stderr.seek(0)
                err = stderr.read().decode('utf-8')
                print(f'\nWarning: sequence {seq_rec.name} has IR on the start, with offset {offset}!\n{err}\n')
                if '161' in err:
                    # Try with longer offset
                    continue
                # continue
                return  # Nothing to do more!

        seq_length = len(seq_rec.seq)
        output = result.stdout.decode('utf-8')
        # Note: output can contain lines of type
        #   found illegal character N at position 1944
        for line in output.splitlines():
            if ';' in line:
                irs = tuple(int(x) - offset for x in line.split(';')[:4])
                ira, irb = _ir(seq_length, *irs[:2]), _ir(seq_length, *irs[2:])

                # Check IR order
                res_irs = (ira, irb) if (irb[0] - ira[1]) % seq_length < (ira[0] - irb[1]) % seq_length else (irb, ira)
                break

    #
    if leave_tmp_file:
        print(f'Note: tmp file {fasta_filename} is not removed!')
    else:
        os.remove(fasta_filename)

    return res_irs


def _ir(seq_length, start, end):
    # Note: irscan output is 1-based, and range means [start..end]
    if end == 0:
        end = seq_length
    return ((start - 1) % seq_length), (end if (end > 0 or end == seq_length) else (end % seq_length))


def small_d_on_file(seq_filename, leave_tmp_file=False, no_prepend_workaround=False, no_dna_fix=False):
    _ext_2_bio_io_type = dict(
        gb='genbank', gbff='genbank',
        fa='fasta',  fas='fasta',
        fastq='fastq',
    )

    base_filename, file_extension = os.path.splitext(seq_filename)
    in_format = _ext_2_bio_io_type[file_extension[1:]]
    return small_d(SeqIO.read(seq_filename, in_format),
                   leave_tmp_file=leave_tmp_file,
                   no_prepend_workaround=no_prepend_workaround,
                   no_dna_fix=no_dna_fix)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Run irscan on a sequence.")
    parser.add_argument('seq_filename', help='Sequence filename')
    parser.add_argument('-T', '--leave-tmp-file', action='store_true', help='Leave temporary file. For testing.')
    parser.add_argument('-P', '--no-prepend-workaround', action='store_true', help='Omit workaround by prepending')
    parser.add_argument('-D', '--no-dna-fix', action='store_true', help='Omit fixing DNA bases')
    params = parser.parse_args()

    print(small_d_on_file(params.seq_filename,
                          leave_tmp_file=params.leave_tmp_file,
                          no_prepend_workaround=params.no_prepend_workaround,
                          no_dna_fix=params.no_dna_fix))
