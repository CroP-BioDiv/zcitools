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
# - - platform: x
#   - read_length: 150
#   - file: r.fastq.gz
# - pair_reads:
# - - platform: x
#   - read_length: 150
#   - insert_length: 300
#   - file_1: r1.fastq.gz
#   - file_2: r2.fastq.gz

_defaults_reads = dict(platform=None, read_length=None)
_defaults_paired_reads = dict(platform=None, read_length=None, insert_length=None)
_omit_files = set(['.adapter.', '.lowqual.'])
_re_12 = re.compile(r'[^0-9][12][^0-9]')
_fasta_gz = re.compile(r'fast[aq]\.gz$')


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

    #
    def __iter__(self):
        for r in self.reads:
            yield 'S', r
        for r in self.paired_reads:
            yield 'PE', r

    #
    def _add_read_files(self, in_arr, r, defaults):
        in_arr.append(dict(defaults))
        in_arr[-1].update(r)

    def add_read(self, r):
        self._add_read_files(self.reads, r, _defaults_reads)

    def add_paired_read(self, r):
        self._add_read_files(self.paired_reads, r, _defaults_paired_reads)

    def write_data(self, filename):
        write_yaml(dict(reads=self.reads, paired_reads=self.paired_reads), filename)

    def print_data(self):
        print_yaml(dict(reads=self.reads, paired_reads=self.paired_reads))

    #
    def add_relative_path(self, _dir):
        for r in self.reads:
            r['file'] = os.path.join(_dir, r['file'])
        for r in self.paired_reads:
            r['file_1'] = os.path.join(_dir, r['file_1'])
            r['file_2'] = os.path.join(_dir, r['file_2'])

    #
    @staticmethod
    def from_file(filename, relative_dir=None):
        sr = SequenceReads(data=read_yaml(filename))
        _dir = os.path.dirname(filename)
        if relative_dir:
            _dir = os.path.join(relative_dir, _dir) if _dir else relative_dir
        if _dir:
            sr.add_relative_path(_dir)
        return sr

    @staticmethod
    def from_directory(directory, platform=None, read_length=None, insert_length=None):
        # Put lot of (file naming) heuristics in this!
        reads = []  # Just files
        paired_reads = []

        files_to_proc = set(f for f in os.listdir(directory)
                            if _fasta_gz.search(f) and
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
            reads=[dict(file=r, platform=platform, read_length=read_length) for r in reads],
            paired_reads=[
                dict(file_1=f1, file_2=f2, platform=platform, read_length=read_length, insert_length=insert_length)
                for f1, f2 in paired_reads])
