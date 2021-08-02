#!/usr/bin/env python3

import sys
import os.path
import shutil
import tempfile
import subprocess
from Bio import SeqIO


def remove_directory(_dir, create):
    # assert not _dir.endswith("/") or _dir.endswith("\\"), _dir
    _dir = os.path.normpath(_dir)

    if os.path.isdir(_dir):
        if sys.platform == "win32":
            temp_path = _dir + "_"

            if os.path.exists(temp_path):
                remove_directory(temp_path, False)

            try:
                os.renames(_dir, temp_path)
            except OSError as exception:
                if exception.errno != errno.ENOENT:
                    raise
            else:
                shutil.rmtree(temp_path)
        else:
            shutil.rmtree(_dir)

    if create:
        os.makedirs(_dir)


def plann(seq_filename, leave_tmp_file=False):
    if not (plann_script := os.environ.get('PLANN_SCRIPT')):
        # Chech is plann.pl executable and on the PATH
        if p_pl := shutil.which('plann.pl'):
            plann_script = os.path.realpath(p_pl)

    if not plann_script:
        raise OSError('Plann script not found! Set PLANN_SCRIPT environment command.')

    # Tmp directories
    tmp_d = tempfile.gettempdir()
    tmp_plann = os.path.join(tmp_d, 'plann')
    remove_directory(tmp_plann, True)

    fa_filename = os.path.join(tmp_plann, 'x.fa')
    SeqIO.convert(seq_filename, 'genbank', fa_filename, 'fasta')

    res_irs = None
    try:
        cmd = ['perl', plann_script, '-reference', seq_filename, '-fasta', fa_filename, '-out', 'plann_result']
        result = subprocess.run(cmd, cwd=tmp_plann, check=True,
                                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        tbl_filename = os.path.join(tmp_plann, 'plann_result.tbl')
        if os.path.isfile(tbl_filename):
            with open(tbl_filename, 'r') as _in:
                irs = [(int(fs[0]), int(fs[1]))
                       for line in _in.readlines()
                       if 'repeat_region' in line and (fs := line.split()) and len(fs) == 3 and fs[2] == 'repeat_region']
                if len(irs) == 2:
                    ira, irb = irs
                    # Our indexing
                    ira = (ira[0] - 1, ira[1])
                    irb = (irb[0] - 1, irb[1])
                    # Find sequence length
                    with open(os.path.join(tmp_plann, 'plann_result.fsa'), 'r') as _fsa:
                        seq_length = len(list(_fsa.readlines())[1]) - 1
                    # Check IR order
                    res_irs = (ira, irb) if (irb[0] - ira[1]) % seq_length < (ira[0] - irb[1]) % seq_length else (irb, ira)

    except subprocess.CalledProcessError:
        print('?' * 100)
        print('PLANN', seq_filename)
        print('?' * 100)

    if not leave_tmp_file:
        remove_directory(tmp_plann, False)

    return res_irs


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Run plann on given GenBank file.")
    parser.add_argument('filename', help='Sequence GenBank filename')
    parser.add_argument('-T', '--leave-tmp-file', action='store_true', help='Leave temporary file. For testing.')

    params = parser.parse_args()
    print(plann(params.filename, leave_tmp_file=params.leave_tmp_file))
