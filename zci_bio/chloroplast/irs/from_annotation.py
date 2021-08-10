#!/usr/bin/env python3

from Bio import SeqIO
from zci_bio.chloroplast.utils import find_chloroplast_irs, ir_loc

"""
Python wrapper around method for finding IRs from annotation.
Script produces same output as other IRs location scripts from this directory.
"""


def from_annotation(filename):
    seq = SeqIO.read(filename, 'genbank')
    if irs := find_chloroplast_irs(seq, check_length=False):
        ira, irb = irs
        return tuple(ir_loc(ira.location.parts)), tuple(ir_loc(irb.location.parts))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Extract IR location from annotation")
    parser.add_argument('filename', help='Sequence GenBank filename')
    params = parser.parse_args()
    print(from_annotation(params.filename))
