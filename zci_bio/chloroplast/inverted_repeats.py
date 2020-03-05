import itertools
from collections import OrderedDict, namedtuple
from common_utils.file_utils import ensure_directory, copy_file, run_module_script, set_run_instructions
from zci_bio.annotations.steps import AnnotationsStep
from common_utils.terminal_layout import StringColumns
from ..utils.import_methods import import_bio_seq_io
from . import run_inverted_repeats
from .run_inverted_repeats import run_one

_Repeat = namedtuple('_Repeat', 'first_start, second_start, length, inverted')

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


def create_irs_data(step_data, input_step, run):
    # Creates Annotations step from input sequences/annotations
    # Steps subdirectory 'run_dir' contains input and output calculation files
    SeqIO = import_bio_seq_io()
    files_to_zip = []

    step = AnnotationsStep(input_step.project, step_data, remove_data=True)
    # Set sequences
    step.set_sequences(input_step.all_sequences())
    ensure_directory(step.step_file('run_dir'))

    for seq_ident, seq_rec in input_step._iterate_records():
        # Write initial genbank file.
        # ToDo: if input step is Annotations than copy would work
        SeqIO.write([seq_rec], step.step_file(f'{seq_ident}.gb'), 'genbank')
        # Set fasta file for calculation
        files_to_zip.append(step.step_file('run_dir', f'{seq_ident}.fa'))
        SeqIO.write([seq_rec], step.step_file('run_dir', f'{seq_ident}.fa'), 'fasta')

    # Stores description.yml
    step.save(completed=run)

    if run:
        run_module_script(run_inverted_repeats, step)
        finish_irs_data(step)
    else:
        set_run_instructions(run_inverted_repeats, step, files_to_zip, _instructions)

    #
    return step


# Finish part
class _MUMmerResult:
    def __init__(self, result_filename, name, length):
        self._name = name
        self._length = length
        repeats = []  # List of _Repeat objects
        with open(result_filename, 'r') as output:
            read = False
            for line in output:
                fields = line.strip().split()
                if read:
                    if fields[1][-1] == 'r':
                        s2 = int(fields[1][:-1])
                        inverted = True
                    else:
                        s2 = int(fields[1])
                        inverted = False
                    repeats.append(_Repeat(int(fields[0]), s2, int(fields[2]), inverted=inverted))
                else:
                    if fields[0] == 'Start1':
                        read = True
        #
        self._repeats = repeats = sorted(repeats, key=lambda r: r.first_start)
        self._a = [(r.first_start, r.first_start + r.length - 1) for r in repeats]
        self._b = [(r.second_start, r.second_start - r.length + 1) if r.inverted else
                   (r.second_start, r.second_start + r.length - 1) for r in repeats]

    def get_irs(self):
        repeats = [r for r in self._repeats if r.inverted]
        if not repeats:
            print(f'Error: no IRs for {self._name}!')
            return

        if len(repeats) == 1:
            irs = repeats[0]
        else:
            # Concatenate repeats into one
            # ToDo: check lot of things!!! Gap, change, length of change, ...
            first = repeats[0]
            last = repeats[-1]
            irs = _Repeat(
                first.first_start, first.second_start,
                last.first_start - first.first_start + last.length,
                inverted=True)

        if irs.length < 20000:
            print(f'Error: IR for {self._name} is of short length {irs.length}!')
            return
        ira = (irs.first_start, irs.first_start + irs.length - 1)
        irb = (irs.second_start - irs.length + 1, irs.second_start)
        if ira[0] < self._length // 2:
            print('ERROR: IRA is located to the left.', self._name, self._length, ira, irb)
        if irb[1] < self._length - 30:
            print('ERROR: IRB is located to the left.', self._name, self._length, ira, irb)
        return ira, irb

    @staticmethod
    def _format_row(data, minus=1):
        x = [f"{f}-{t}" for f, t in data]
        y = [str(minus * (a[0] - b[1])) for a, b in zip(data[1:], data[:-1])]
        return list(itertools.chain.from_iterable(zip(x, y))) + x[-1:]

    def show_data(self):
        irs = self.get_irs()
        if irs:
            ira, irb = irs
            irs = f'IRS: {ira[0]}-{ira[1]} : {irb[0]}-{irb[1]}'
        else:
            irs = 'No IRS!'

        print(f"{self._name} ({self._length})  {irs}")
        print(StringColumns([self._format_row(self._a), self._format_row(self._b, minus=-1)]))


def _set_irs_to_seq_rec(seq_rec, irs):
    from Bio.SeqFeature import SeqFeature, FeatureLocation  # If SeqIO exists than these should be installed also
    ira, irb = irs
    qs = OrderedDict([('rpt_type', 'inverted')])
    seq_rec.features.append(SeqFeature(FeatureLocation(*ira), type='repeat_region', qualifiers=qs))
    seq_rec.features.append(SeqFeature(FeatureLocation(*irb), type='repeat_region', qualifiers=qs))


def finish_irs_data(step_obj):
    SeqIO = import_bio_seq_io()

    for seq_ident, seq_rec in step_obj._iterate_records():
        # Read MUMmer output
        m_res = _MUMmerResult(step_obj.step_file('run_dir', f'{seq_ident}.out'), seq_ident, len(seq_rec))
        irs = m_res.get_irs()
        if irs:
            _set_irs_to_seq_rec(seq_rec, irs)
            SeqIO.write([seq_rec], step_obj.step_file(f'{seq_ident}.gb'), 'genbank')

    step_obj.save()


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
    irs = m_res.get_irs()
    if irs:
        _set_irs_to_seq_rec(seq_rec, irs)
        return True
