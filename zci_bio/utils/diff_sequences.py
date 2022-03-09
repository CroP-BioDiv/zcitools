import difflib
try:
    import edlib
    import re
except ImportError:
    edlib = None
try:
    import Levenshtein
except ImportError:
    Levenshtein = None

"""
Implements classes with same interface for comparing sequences.
Methods:
 - get_opcodes(self)
   Returns list of operations needed to convert first sequence into the second.
   Format is same as in method get_opcode() from standard implemented in difflib library.
   https://docs.python.org/3/library/difflib.html#difflib.SequenceMatcher.get_opcodes
   List of 5-tuples describing how to turn a into b.
   Each tuple is of the form (tag, i1, i2, j1, j2).
   The first tuple has i1 == j1 == 0, and remaining tuples have i1 equal to the i2 from the preceding tuple,
   and, likewise, j1 equal to the previous j2.
 - in_short(self)
   Returns string of format [R:<num_blocks>,<num_characters>][I:<num_blocks>,<num_characters>].
   Note: Inserts and deletes are threated as one (I).

"""


class Diff_difflib:
    def __init__(self, a, b):
        diff = difflib.SequenceMatcher(a=a, b=b, autojunk=False)
        self.opcodes = diff.get_opcodes()

    def get_opcodes(self):
        return self.opcodes

    def in_short(self):
        # R, I -> [num blocks, num bps]
        RM = [0, 0]
        IN = [0, 0]
        for x in self.opcodes:
            print('  ', x)
            if x[0] == 'equal':
                continue
            if x[0] == 'replace':
                c, bp = RM, (x[2] - x[1])
            elif x[0] == 'delete':
                c, bp = IN, (x[2] - x[1])
            elif x[0] == 'insert':
                c, bp = IN, (x[4] - x[3])
            c[0] += 1
            c[1] += bp
        return ';'.join(f'{label}:{d[0]},{d[1]}' for d, label in ((RM, 'R'), (IN, 'I')) if d[0])


class Diff_Levenshtein(Diff_difflib):
    def __init__(self, a, b):
        self.opcodes = Levenshtein.opcodes(a, b)


class Diff_edlib(Diff_difflib):
    # Note: edlib has no gap penalty :-/
    # That means it likes to split longer inserts/deletes into smaller ones by adding short equals in-between.
    def __init__(self, a, b):
        diff = edlib.align(a, b, task='path', mode='SHW')
        self.opcodes = self._to_difflib(diff['cigar'])

    def _to_difflib(self, cigar):
        # From doc https://pypi.org/project/edlib/
        # Here we are using extended cigar format, which uses following symbols:
        # Match: '=', Insertion to target: 'I', Deletion from target: 'D', Mismatch: 'X'.
        # e.g. cigar of "5=1X1=1I" means "5 matches, 1 mismatch, 1 match, 1 insertion (to target)".

        a_idx = b_idx = 0
        opcodes = []
        for mo in re.finditer(r'\d+', cigar):
            num_c = int(mo.group(0))
            op_c = cigar[mo.span()[1]]
            if op_c == '=':
                opcodes.append(('equal', a_idx, a_idx + num_c, b_idx, b_idx + num_c))
                a_idx += num_c
                b_idx += num_c
            elif op_c == 'D':
                opcodes.append(('insert', a_idx, a_idx, b_idx, b_idx + num_c))
                b_idx += num_c
            elif op_c == 'I':
                opcodes.append(('delete', a_idx, a_idx + num_c, b_idx, b_idx))
                a_idx += num_c
            elif op_c == 'X':
                opcodes.append(('replace', a_idx, a_idx + num_c, b_idx, b_idx + num_c))
                a_idx += num_c
                b_idx += num_c
            else:
                assert False, (op_c, cigar)
        return opcodes


if Levenshtein:
    Diff = Diff_Levenshtein
elif edlib:
    Diff = Diff_edlib
else:
    print('Note: comparing of sequences backed up on standard difflib library. Which is quite slow.\n' +
          'To make it faster, install python-Levenshtein or edlib!')
    Diff = Diff_difflib


def diff_check_memory(a, b):
    if Levenshtein:
        # Note: Levenshtein methods use much more memory!!!
        from common_utils.resources import virtual_memory
        if 10 * len(a) * len(b) > 0.9 * virtual_memory:
            print('Using difflib!')
            return Diff_difflib(a, b)
        return Diff_Levenshtein(a, b)
    return Diff_difflib(a, b)


if __name__ == '__main__':
    # Test by checking IRs for given GenBank file
    import sys
    from Bio import SeqIO
    from zci_bio.chloroplast.utils import find_chloroplast_irs
    seq = SeqIO.read(sys.argv[1], 'genbank')
    if irs := find_chloroplast_irs(seq, check_length=False):
        ira, irb = irs
        ira_s = ira.extract(seq)
        irb_s = irb.extract(seq)
        if ira.strand == irb.strand:
            irb_s = irb_s.reverse_complement()
        # diff = Diff_edlib(str(ira_s.seq), str(irb_s.seq))
        ira = str(ira_s.seq)
        irb = str(irb_s.seq)
        print(f'Lengths {len(ira)} and {len(irb)}')
        diff_cls = Diff
        if len(sys.argv) > 2:
            char = sys.argv[2][0].upper()
            if char == 'L':
                diff_cls = Diff_Levenshtein
            elif char == 'E':
                diff_cls = Diff_edlib
            elif char == 'D':
                diff_cls = Diff_difflib
        print('cls:', diff_cls.__name__)
        diff = diff_cls(ira, irb)
        print(diff.get_opcodes())
    else:
        print('No IRs found!!!')
