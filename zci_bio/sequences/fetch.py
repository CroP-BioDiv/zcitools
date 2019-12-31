from .steps import SequencesStep
from ..utils.entrez import Entrez


def fetch_sequences(step_data, table_step, cache):
    step = SequencesStep(table_step.project, step_data, update_mode=True)
    all_sequences = [ni for ni in table_step.get_column_values_by_type('seq_ident') if not step.sequence_exists(ni)]

    # Fetch cached sequences
    to_fetch = step.get_cached_records(cache, all_sequences, info=True)

    if to_fetch:
        entrez = Entrez()

        for ni in to_fetch:
            assert ni.startswith('NC_'), 'For now only NCBI donwload is supported!!!'  # ToDo:
            #
            print(f"  fetching '{ni}' in GenBank format.")
            # Fetch genbank format
            filename = step.step_file(ni + '.gb')
            entrez.efetch(filename, db='nuccore',  id=ni, rettype='gb', retmode='text')
            if cache:
                cache.set_record(ni, filename)

    for ni in all_sequences:
        step.add_sequence_file(ni + '.gb')

    # ToDo: remove not referenced sequences

    # Store step data
    step._check_data()
    step.save()
    return step
