import os.path
from .steps import SequencesStep
from ..utils.entrez import Entrez
from ..utils.helpers import fetch_our_sequence
from common_utils.file_utils import write_str_in_file

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

    # Fetch from our sequences
    all_sequences = []
    for ni in table_step.get_column_values_by_type('seq_ident', column_name=column_name):
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

    # ToDo: remove not referenced sequences

    # Store step data
    # step._check_data()
    step.save(completed=not to_fetch)
    if to_fetch:
        write_str_in_file(
            step.step_file('INSTRUCTIONS.txt'),
            _instructions_no_data.format(sequence_db=sequence_db, seqs=', '.join(sorted(to_fetch))))
    return step
