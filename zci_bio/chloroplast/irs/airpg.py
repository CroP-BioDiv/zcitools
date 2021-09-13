#!/usr/bin/env python3

from Bio import SeqIO

"""
Wrapper around airpg module ir_operations.py.
Repository: https://github.com/michaelgruenstaeudl/airpg
File: airpg/ir_operations.py
      https://github.com/michaelgruenstaeudl/airpg/blob/main/airpg/ir_operations.py
"""


def airpg_filename(gb_filename, min_ir_length=1000, ret_features=False):
    return airpg(SeqIO.read(gb_filename, 'genbank'), min_ir_length=min_ir_length, ret_features=ret_features)


def airpg(seq_rec, min_ir_length=1000, ret_features=False):
    try:
        from .ir_operations import IROperations
    except ImportError:
        from ir_operations import IROperations
    iro = IROperations()
    try:
        ira, irb = iro.identify_inverted_repeats(seq_rec, min_IR_len=min_ir_length)
    except Exception as r:
        return None

    # Note: NC_039346 has one repeat of length 2bp, and airpg doesn't locate it good.
    if not ira or not irb:
        return

    # Check orientation
    diff_1 = (irb.location.parts[0].start - ira.location.parts[-1].end) % len(seq_rec)
    diff_2 = (ira.location.parts[0].start - irb.location.parts[-1].end) % len(seq_rec)
    if diff_1 > diff_2:
        ira, irb = irb, ira

    if ret_features:
        return ira, irb
    return _loc(ira), _loc(irb)


def _loc(ir):
    ir = ir.location.parts
    assert len(ir) <= 2, ir
    if ir[0].strand > 0:
        return (int(ir[0].start), int(ir[-1].end))
    return (int(ir[-1].start), int(ir[0].end))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Run airpg on given GenBank file.")
    parser.add_argument('filename', help='Sequence GenBank filename')
    parser.add_argument('-l', '--min-ir-length', default=1000, type=int, help='Minimal IR length')
    params = parser.parse_args()
    print(airpg_filename(params.filename, min_ir_length=params.min_ir_length))
