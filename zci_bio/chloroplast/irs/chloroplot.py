#!/usr/bin/env python3

import os.path
import subprocess

"""
Python wrapper method around chloroplot.R wrapper :-)
"""


def chloroplot(genbank_file, print_chloroplot_output=False):
    _dir = os.path.dirname(os.path.abspath(__file__))
    r_script = os.path.join(_dir, 'chloroplot.R')
    try:
        result = subprocess.run(['Rscript', r_script, genbank_file],
                                check=True,
                                stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        return

    output = result.stdout.decode('utf-8')
    if print_chloroplot_output:
        print(output)
    return parse_output(output.split('\n'))


def parse_output(lines):
    # Format of output is:
    # $ir_table
    #    chr  start    end name       text center
    # 1 chr1      0  81884  LSC LSC: 81934  40917
    # 2 chr1  81884 107546  IRA IRA: 25662  94715
    # 3 chr1 107546 124801  SSC SSC: 17255 116174
    # 4 chr1 124801 150460  IRB IRB: 25659 137630
    # 5 chr1 150460 150510  LSC            150485
    #
    # $indel_table
    #    mismatch_type position string    col
    # 1         insert    82161      C  green
    # 2         insert    82162      C  green
    # 3         insert    82163      C  green
    # 4         delete   150187      D yellow
    # 5        replace   104414      G    red
    # 6        replace   104420      G    red
    # 7        replace   104421      C    red
    # ...

    iras, irbs = [], []
    seq_length = 0
    for line in lines:
        fields = line.split()
        if len(fields) < 5:
            continue
        if fields[4] in ('IRA', 'IRB', 'LSC', 'SSC'):
            seq_length = int(fields[3])
            if fields[4] == 'IRA':
                iras.append((int(fields[2]), int(fields[3])))
            if fields[4] == 'IRB':
                irbs.append((int(fields[2]), int(fields[3])))
    if iras and irbs:
        ira, irb = _loc(iras), _loc(irbs)
        return (ira, irb) if ((irb[0] - ira[1]) % seq_length < (ira[0] - irb[1]) % seq_length) else (irb, ira)


def _loc(ir_parts):
    assert 1 <= len(ir_parts) <= 2, ir_parts
    if len(ir_parts) == 1:
        return tuple(ir_parts[0])
    # Of type:
    #    chr  start    end name       text center
    # 1 chr1      0      3  IRB IRB: 24645 139014
    # 2 chr1      3  83717  LSC LSC: 83714  41860
    # 3 chr1  83717 108362  IRA IRA: 24645  96040
    # 4 chr1 108362 126691  SSC SSC: 18329 117526
    # 5 chr1 126691 151333  IRB            139012
    return ir_parts[1][0], ir_parts[0][1]


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Run Chlorolpot IR detection on sequence stored as genbank file.")
    parser.add_argument('genbank_file', help='Sequence filename')
    parser.add_argument('-C', '--print-chloroplot-output', action='store_true', help='Print Chloroplot R output.')
    params = parser.parse_args()
    print(chloroplot(params.genbank_file, print_chloroplot_output=params.print_chloroplot_output))

    # # Parse test
    # import sys
    # print(parse_output(list(open(sys.argv[1]).readlines())))
