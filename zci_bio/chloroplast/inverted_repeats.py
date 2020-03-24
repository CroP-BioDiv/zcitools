import itertools
import os.path
from collections import OrderedDict, namedtuple
from common_utils.file_utils import ensure_directory, copy_file, run_module_script, set_run_instructions, write_yaml
from zci_bio.annotations.steps import AnnotationsStep
from common_utils.terminal_layout import StringColumns
from ..utils.helpers import read_sequence
from ..utils.import_methods import import_bio_seq_io
from . import run_inverted_repeats
from .run_inverted_repeats import run_one


class _RawData(namedtuple('_RawData', 'start_1, end_1, start_2, end_2, match_length, inverted')):
    def in_annotation(self, rpt_type):
        return _Annotation(rpt_type, self.start_1, self.end_1, self.start_2, self.end_2, self.match_length)

    def is_in(self, b):
        return self.start_1 <= b.start_1 <= self.end_1 and \
            self.start_1 <= b.end_1 <= self.end_1 and \
            self.start_2 <= b.start_2 <= self.end_2 and \
            self.start_2 <= b.end_2 <= self.end_2

    def is_after(self, b, max_gap):
        assert self.inverted == b.inverted, (self, b)
        print(abs(b.start_1 - self.end_1) <= max_gap, abs(b.end_2 - self.start_2) <= max_gap, self[:4], b[:4])
        if abs(b.start_1 - self.end_1) <= max_gap:  # Check first repeat
            if self.inverted:                       # Check second repeat
                return abs(b.end_2 - self.start_2) <= max_gap
            return abs(b.start_2 - self.end_2) <= max_gap

    def add_after(self, b):
        assert self.inverted == b.inverted, (self, b)
        match = self.match_length + b.match_length
        if self.inverted:
            match -= min(0, b.start_1 - self.end_1, self.start_2 - b.end_2)
            return _RawData(self.start_1, b.end_1, b.start_2, self.end_2, match, True)
        match -= min(0, b.start_1 - self.end_1, b.start_2 - self.end_2)
        return _RawData(self.start_1, b.end_1, self.start_2, b.end_2, match, False)

    # Cycle merging
    def is_after_end(self, b, max_gap, length):
        assert self.inverted == b.inverted, (self, b)
        # Example of Mummer format in this case:
        # inverted:
        #   83739  152567r  25076  (self)
        #   1      83738r   114    (b)
        # not inverted:
        #   10000  127491   25076  (self)
        #   1       35076   114    (b)
        if length - self.end_2 + b.start_1 <= max_gap:  # Check second repeat
            if self.inverted:                           # Check second repeat
                return length - self.end_2 + b.start_1 <= max_gap
            return abs(self.end_1 - b.start_2) <= max_gap


# Repeat type, one of values: inverted (for IRa and IRb), direct (short direct), other (for other irs)
# Check http://www.insdc.org/controlled-vocabulary-rpttype-qualifier for possible values
_Annotation = namedtuple('_Annotation', 'rpt_type, start_1, end_1, start_2, end_2, matched')

_instructions = """
Steps:
 - copy file calculate.zip onto server
 - unzip it
 - change directory to {step_name}
 - run script: python3 {script_name}
    - to specify number of threads to use run: python3 {script_name} <num_threads>
      default is number of cores.
 - copy file output.zip back into project's step directory {step_name}
 - run zcit command: zcit.py finish {step_name}

Notes:
 - MUMmers repeat match executable (repeat-match) should be on the PATH or
   environment variable MUMMER_REPEAT_MATCH_EXE should point to it.
 - It is good to use command screen for running the script.
   screen -dm "python3 {script_name}"
"""


def create_irs_data(step_data, input_step, params):  # , run):
    # Creates Annotations step from input sequences/annotations
    # Steps subdirectory 'run_dir' contains input and output calculation files
    SeqIO = import_bio_seq_io()
    files_to_zip = []
    calc_seq_idents = []

    step = input_step.project.new_step(AnnotationsStep, step_data)
    # Set sequences
    step.set_sequences(input_step.all_sequences())
    ensure_directory(step.step_file('run_dir'))

    for seq_ident in input_step.all_sequences():
        out_file = step.step_file('run_dir', f'{seq_ident}.out')
        if not os.path.isfile(out_file):
            seq_rec = input_step.get_sequence_record(seq_ident)
            # Set fasta file for calculation
            files_to_zip.append(step.step_file('run_dir', f'{seq_ident}.fa'))
            SeqIO.write([seq_rec], files_to_zip[-1], 'fasta')
            calc_seq_idents.append(seq_ident)
        elif not os.path.isfile(step.step_file(f'{seq_ident}.gb')):
            calc_seq_idents.append(seq_ident)

    if files_to_zip:
        # Store finish.yml
        finish_f = step.step_file('finish.yml')
        write_yaml(dict(fa_files=files_to_zip), finish_f)

        run = True  # ToDo: ...
        step.save(completed=False)
        if run:
            run_module_script(run_inverted_repeats, step)
            finish_irs_data(step, calc_seq_idents=calc_seq_idents)
        else:
            files_to_zip.append(finish_f)
            set_run_instructions(run_inverted_repeats, step, files_to_zip, _instructions)
    #
    elif calc_seq_idents:
        finish_irs_data(step, calc_seq_idents=calc_seq_idents)
        step.save()
    elif params.force_mummer_parse:
        finish_irs_data(step)
        step.save()

    #
    return step


# Finish part
def finish_irs_data(step_obj, calc_seq_idents=None):
    if calc_seq_idents is None:
        calc_seq_idents = [f[:-4] for f in step_obj.step_dir_files('run_dir') if f.endswith('.out')]
    if calc_seq_idents:
        SeqIO = import_bio_seq_io()
        for seq_ident in calc_seq_idents:
            seq_rec = read_sequence(step_obj.step_file('run_dir', f'{seq_ident}.fa'))
            # Read MUMmer output file
            m_res = _MUMmerResult(step_obj.step_file('run_dir', f'{seq_ident}.out'), seq_ident, len(seq_rec))
            m_res.set_annotation(seq_ident, seq_rec, step_obj.step_file(f'{seq_ident}.gb'))

    step_obj.save(completed=True)


class _MUMmerResult:
    def __init__(self, result_filename, name, length):
        self._name = name
        self._length = length
        self._raw_repeats = []  # List of _RawData objects, sorted by start position
        self._annotation = []   # List of _Annotation objects
        self._irs = None        # _Annotation object or None
        self._read_input_data(result_filename)
        #
        rest_inverted_r = self._annotate_irs([r for r in self._raw_repeats if r.inverted])
        self._annotation.extend(r.in_annotation('other') for r in rest_inverted_r)
        self._annotation.extend(r.in_annotation('direct') for r in self._raw_repeats if not r.inverted)

    def _read_input_data(self, result_filename):
        repeats = []
        with open(result_filename, 'r') as output:
            read = False
            for line in output:
                fields = line.strip().split()
                if read:
                    inverted = (fields[1][-1] == 'r')
                    length = int(fields[2])
                    start_1 = int(fields[0])
                    start_2 = int(fields[1][:-1] if inverted else fields[1])
                    if inverted:
                        start_2 = start_2 - length + 1
                    repeats.append(_RawData(
                        start_1, start_1 + length - 1,
                        start_2, start_2 + length - 1,
                        length, inverted))
                else:
                    if fields[0] == 'Start1':
                        read = True
        self._raw_repeats = sorted(repeats, key=lambda r: r.start_1)

    def _annotate_irs(self, inverted_r):
        if not inverted_r:  # Nothing to process!
            print('ERROR: no inverted repeats at all!', self._name, self._length)
            return inverted_r

        # Note: Mummer was run with n=100, which means if two small gaps are on distance <100,
        #       than gap can be of length ~100
        max_gap = 120
        next_i = 0
        concatenated = []  # Tuples (concatenated _RawData, list of used indices)
        while next_i < len(inverted_r):
            _irs = inverted_r[next_i]
            from_i = next_i
            next_i += 1
            while next_i < len(inverted_r):
                if _irs.is_in(inverted_r[next_i]):  # Next one can be substring of current
                    next_i += 1
                    continue
                if not _irs.is_after(inverted_r[next_i], max_gap):
                    break
                _irs = _irs.add_after(inverted_r[next_i])
                next_i += 1
            concatenated.append((_irs, list(range(from_i, next_i))))

        # Check end-start concatenation
        if len(concatenated) > 1:
            l_irs, l_ids = concatenated[-1]
            f_irs, f_ids = concatenated[0]
            if l_irs.is_after_end(f_irs, max_gap, self._length):
                concatenated[0] = (l_irs.add_after(f_irs), l_ids + f_ids)
                concatenated.pop()

        # Find longest match
        irs, used_ids = max(concatenated, key=lambda x: x[0].match_length)

        # Check data
        if irs.match_length < 20000:
            print(f'WARNING: IR for {self._name} is of short length {irs.match_length}!')
        ira = (irs.start_1, irs.end_1)
        irb = (irs.start_2, irs.end_2)
        if ira[0] < self._length // 2:
            print('WARNING: IRA is located to the left.', self._name, self._length, ira, irb)
        if irb[1] < self._length - 30:
            print('WARNING: IRB is located to the left.', self._name, self._length, ira, irb)

        #
        self._annotation.append(irs.in_annotation('inverted'))
        self._irs = irs
        return [x for i, x in enumerate(inverted_r) if i not in used_ids]

    #
    def set_annotation(self, seq_ident, seq_rec, gb_file):
        if self._annotation:
            SeqIO = import_bio_seq_io()
            # If SeqIO exists than these should be installed!
            from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
            from Bio.Alphabet import DNAAlphabet
            _seg_length = lambda s, e: (e - s) if e > s else (s + self._length - e)
            _get_loc = lambda s, e: FeatureLocation(s, e) if e > s else \
                CompoundLocation([FeatureLocation(s, self._length), FeatureLocation(0, e)])
            for a in self._annotation:
                length = max(_seg_length(a.start_1, a.end_1), _seg_length(a.start_1, a.end_1))
                m = 100 if a.matched == length else round(100 * length / a.matched, 2)
                note_1 = f'repeat_hit position {a.start_2} - {a.end_2}, match {m}.%'
                note_2 = f'repeat_hit position {a.start_1} - {a.end_1}, match {m}.%'
                rpt = ('rpt_type', a.rpt_type)
                seq_rec.features.append(
                    SeqFeature(_get_loc(a.start_1, a.end_1), type='repeat_region',
                               qualifiers=OrderedDict([rpt, ('note', note_1)])))
                seq_rec.features.append(
                    SeqFeature(_get_loc(a.start_2, a.end_2), type='repeat_region',
                               qualifiers=OrderedDict([rpt, ('note', note_2)])))
            # Fix other data
            seq_rec.id = seq_ident
            seq_rec.seq.alphabet = DNAAlphabet()
            SeqIO.write([seq_rec], gb_file, 'genbank')
            return True

    #
    @staticmethod
    def _format_row(data):
        x = [f"{f}-{t}" for f, t, _ in data]
        y = [str((-1 if b[2] else 1) * (a[0] - b[1])) for a, b in zip(data[1:], data[:-1])]
        return list(itertools.chain.from_iterable(zip(x, y))) + x[-1:]

    def show_data(self):
        if self._irs:
            i = self._irs
            irs = f'IRS: {i.start_1}-{i.end_1} : {i.start_2}-{i.end_2}'
        else:
            irs = 'No IRS!'

        print(f"{self._name} ({self._length})  {irs}")
        _a = [(r.start_1, r.end_1, False) for r in _raw_repeats]
        _b = [(r.start_2, r.end_2, r.inverted) for r in _raw_repeats]
        print(StringColumns([self._format_row(_a), self._format_row(_b)]))


#
def show_irs_data(step_obj):
    for seq_ident, seq_rec in step_obj._iterate_records():
        m_res = _MUMmerResult(step_obj.step_file('run_dir', f'{seq_ident}.out'), seq_ident, len(seq_rec))
        m_res.show_data()


#
def calculate_and_add_irs_to_seq_rec(step, seq_ident, seq_rec):
    SeqIO = import_bio_seq_io()
    # Store input fasta file
    ensure_directory(step.step_file('run_dir'))
    input_filename = step.step_file('run_dir', f'{seq_ident}.fa')
    SeqIO.write([seq_rec], input_filename, 'fasta')
    # Run MUMmer
    run_one(input_filename)
    m_res = _MUMmerResult(step.step_file('run_dir', f'{seq_ident}.out'), seq_ident, len(seq_rec))
    return m_res.set_annotation(seq_ident, seq_rec, step.step_file(f'{seq_ident}.gb'))
