import os
import re
from common_utils.file_utils import write_yaml, read_yaml, print_yaml

# Sequencing can be done on different platforms with different read properties.
# More data files are produced and used for later processing (assembling.)
#
# Sequence read description described (lists) raw data files with there properties.
# It is used to store data in clean form. This data can be used as input parameters for tools (assemblers.)
# Data for single reads and paired reads are stored, with there properties.
# - reads
# - - instrument: x
#   - read_length: 150
#   - files: [r.fastq.gz, ...]
# - pair_reads:
# - - instrument: x
#   - read_length: 150
#   - gap_length: 300
#   - files: [[r1.fastq.gz, r2.fastq.gz], ...]

_defaults_reads = dict(instrument=None, read_length=None)
_defaults_paired_reads = dict(instrument=None, read_length=None, gap_length=None)
_omit_files = set(['.adapter.', '.lowqual.'])
_re_12 = re.compile(r'[^0-9][12][^0-9]')


class SequenceReads:
    def __init__(self, reads=None, paired_reads=None, data=None):
        self.reads = []
        self.paired_reads = []
        #
        for r in reads or []:
            self.add_read(r)
        for r in paired_reads or []:
            self.add_paired_read(r)
        if data:
            for r in data.get('reads', []):
                self.add_read(r)
            for r in data.get('paired_reads', []):
                self.add_paired_read(r)

    def _add_read_files(self, in_arr, r, defaults):
        d = dict(defaults)
        d.update(r)
        d['files'] = list(d['files'])  # Copy it!
        in_arr.append(d)

    def add_read(self, r):
        self._add_read_files(self.reads, r, _defaults_reads)

    def add_paired_read(self, r):
        self._add_read_files(self.paired_reads, r, _defaults_paired_reads)

    def write_data(self, filename):
        write_yaml(dict(reads=self.reads, paired_reads=self.paired_reads), filename)

    def print_data(self):
        print_yaml(dict(reads=self.reads, paired_reads=self.paired_reads))

    #
    @staticmethod
    def from_file(filename):
        data = read_yaml(filename)
        return SequenceReads(data=data)

    @staticmethod
    def from_directory(directory, instrument=None, read_length=None, gap_length=None):
        # Put lot of (file naming) heuristics in this!
        reads = []  # Just files
        paired_reads = []

        files_to_proc = set(f for f in os.listdir(directory)
                            if f.endswith('fastq.gz') and
                            all(x not in f for x in _omit_files) and
                            os.path.isfile(os.path.join(directory, f)))

        for f in list(sorted(files_to_proc)):  # Iterate through set copy
            if f not in files_to_proc:
                continue
            files_to_proc.discard(f)  # To be sure :-)

            # Find {1|2} between non numbers
            for m in _re_12.finditer(f):
                s1, s2 = m.span()
                other_f = f"{f[:s1 + 1]}{2 if f[s1 + 1] == '1' else 1}{f[s2 - 1:]}"
                if other_f in files_to_proc:
                    files_to_proc.discard(other_f)
                    paired_reads.append([f, other_f])
                else:
                    reads.append(f)

        #
        return SequenceReads(
            reads=[dict(files=reads, instrument=instrument, read_length=read_length)] if reads else None,
            paired_reads=[
                dict(files=paired_reads, instrument=instrument, read_length=read_length, gap_length=gap_length)]
            if paired_reads else None)
