from zcitools.steps.sequences import SequencesStep
from zcitools.utils.import_methods import import_bio_entrez


def fetch_sequences(step_data, table_step):
    step = SequencesStep(step_data, update_mode=True)
    cache = step.get_cache_object()
    all_sequences = [ni for ni in table_step.get_column_by_type('seq_ident') if not step.sequence_exists(ni)]

    # Fetch cached sequences
    to_fetch = step.get_cached_records(cache, all_sequences, info=True)

    if to_fetch:
        Entrez = import_bio_entrez()

        for ni in to_fetch:
            assert ni.startswith('NC_'), 'For now only NCBI donwload is supported!!!'  # ToDo:
            #
            print(f"  fetching '{ni}' in GenBank format.")
            # Fetch genbank format
            filename = step.step_file(ni + '.gb')
            with Entrez.efetch(db='nuccore',  id=ni, rettype='gb', retmode='text') as handle:
                with open(filename, 'w') as out:
                    out.write(handle.read())
            if cache:
                cache.set_record(ni, filename)

    for ni in all_sequences:
        step.add_sequence_file(ni + '.gb')

    # ToDo: remove not referenced sequences

    # Store step data
    step._check_data()
    step.save()
    return step
