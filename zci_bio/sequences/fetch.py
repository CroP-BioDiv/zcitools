import os.path
import tempfile
from .steps import SequencesStep
from ..utils.entrez import Entrez
from ..utils.helpers import fetch_our_sequence
from ..utils.import_methods import import_bio_seq_io
from common_utils.file_utils import write_str_in_file, silent_remove_file, remove_directory

_instructions_no_data = """
Not all data is fetched!

Sequences not fetched: {seqs}

Reason can be:
 - not found in sequence database ({sequence_db}).
 - not found as NCBI sequence.
"""


def fetch_sequences(step_data, table_step, common_db, column_name=None):
    step = SequencesStep(table_step.project, step_data, remove_data=True)
    table_step.propagate_step_name_prefix(step)

    seq_idents = table_step.get_column_values_by_type('seq_ident', column_name=column_name)
    to_fetch = do_fetch_sequences(step, seq_idents, common_db)

    # ToDo: remove not referenced sequences

    # Store step data
    # step._check_data()
    step.save(completed=not to_fetch)
    if to_fetch:
        write_str_in_file(
            step.step_file('INSTRUCTIONS.txt'),
            _instructions_no_data.format(sequence_db=sequence_db, seqs=', '.join(sorted(to_fetch))))
    return step


def do_fetch_sequences(step, seq_idents, common_db):
    # Fetch from our sequences
    all_sequences = []
    for ni in sorted(seq_idents):
        if not step.sequence_exists(ni):
            # Fetching our sequences only for base DB
            ext = fetch_our_sequence(ni, step.directory)
            if ext:
                step.add_sequence_file(ni + ext)
            else:
                all_sequences.append(ni)

    # Fetch from cached common_db sequences
    to_fetch = step.get_common_db_records(common_db, all_sequences, info=True)

    if to_fetch:
        entrez = Entrez()

        fetch_nis = to_fetch
        to_fetch = []
        for ni in fetch_nis:
            assert ni.startswith('NC_'), f'For now only NCBI donwload is supported ({ni})!!!'  # ToDo:
            #
            print(f"  fetching '{ni}' in GenBank format.")
            # Fetch genbank format
            filename = step.step_file(ni + '.gb')
            entrez.efetch(filename, db='nuccore',  id=ni, rettype='gb', retmode='text')
            if os.path.isfile(filename):
                if common_db:
                    common_db.set_record(ni, filename)
            else:
                to_fetch.append(ni)

    #
    for ni in set(all_sequences) - set(to_fetch):
        step.add_sequence_file(ni + '.gb')

    return to_fetch


def update_cached_sequences(common_db_base, start_seq_ident, end_seq_ident):
    cdb_seqs = common_db_base.get_relative_db('sequences')
    to_fetch = sorted(x[0] for x in cdb_seqs.get_all_record_ident(startswith='NC_'))
    if start_seq_ident:
        to_fetch = [x for x in to_fetch if x >= start_seq_ident]
    if end_seq_ident:
        to_fetch = [x for x in to_fetch if x <= end_seq_ident]

    if not to_fetch:
        print('Nothing to fetch!')
        return
    print(f'To fetch {len(to_fetch)}')

    entrez = Entrez()
    SeqIO = import_bio_seq_io()
    tmp_fetched_dir = os.path.join(tempfile.gettempdir(), 'zcitools_fatched')
    remove_directory(tmp_fetched_dir, create=True)
    tmp_cached_dir = os.path.join(tempfile.gettempdir(), 'zcitools_cached')
    remove_directory(tmp_cached_dir, create=True)

    for seq_ident in to_fetch:
        # Fetch data
        fetched_filename = os.path.join(tmp_fetched_dir, f'{seq_ident}.gb')
        entrez.efetch(fetched_filename, db='nuccore',  id=seq_ident, rettype='gb', retmode='text')
        if not os.path.isfile(fetched_filename):
            print(f'Warning: no NCBI sequence {seq_ident}!!!')
            continue

        with open(fetched_filename, 'r') as _in:
            fetched_data = _in.read()

        # Cached data
        cdb_seqs.get_record(seq_ident, tmp_cached_dir)
        cached_filename = os.path.join(tmp_cached_dir, f'{seq_ident}.gb')
        with open(cached_filename, 'r') as _in:
            cached_data = _in.read()

        # If data is the same, do nothing
        if fetched_data == cached_data:
            print(f'  {seq_ident}: no change')
            continue

        fetched_seq = SeqIO.read(fetched_filename, 'genbank')
        cached_seq = SeqIO.read(cached_filename, 'genbank')

        # If sequences are not equal, than remove dependent files
        if str(fetched_seq.seq) != str(cached_seq.seq):
            print(f'  {seq_ident}: remove dependent files!')
            for idents in _cdb_dependent_files:
                common_db_base.remove_record(idents)
        #
        print(f'  {seq_ident}: updated cached file!')
        cdb_seqs.set_record(seq_ident, fetched_filename, force=True, remove_directories=True)


def _cdb_dependent_files(seq_ident):
    # sequences: f_, n_, p_
    # GeSeq: NC_, f_, n_, p_
    num = seq_ident[3:]  # Remove NC_
    return [('sequences', f'{x}_{num}') for x in ('f', 'n', 'p')] + \
        [('GeSeq', f'{x}_{num}') for x in ('NC', 'f', 'n', 'p', )]
