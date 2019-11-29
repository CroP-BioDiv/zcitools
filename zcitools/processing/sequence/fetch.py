from zcitools.steps.sequences import SequencesStep
from zcitools.utils.import_methods import import_bio_entrez


def fetch_sequences(step_data, table_step, force_download=False):
    Entrez = import_bio_entrez()

    step = SequencesStep(step_data, update_mode=True)
    cache = step.get_cache_object()
    for ni in table_step.get_column_by_type('seq_ident'):
        assert ni.startswith('NC_'), 'For now only NCBI donwload is supported!!!'  # ToDo:

        if force_download or not step.sequence_exists(ni):
            ni_filename = ni + '.gb'
            if cache and cache.has_record(ni):
                cache.get_record(ni, step._step_name, info=True)
                step.add_sequence_file(ni_filename)
                continue

            #
            print(f"  fetching '{ni}' in GenBank format.")
            # ToDo: what with db='nuccore'?
            filename = step.step_file(ni_filename)  # Fetch genbank format
            with Entrez.efetch(db='nuccore',  id=ni, rettype='gb', retmode='text') as handle:
                with open(filename, 'w') as out:
                    out.write(handle.read())
            #
            step.add_sequence_file(ni_filename)
            if cache:
                cache.set_record(ni, step.step_file(ni_filename))
        else:
            print(f"  data for '{ni}' already presented.")

    # ToDo: remove not referenced sequences

    # Store step data
    step._check_data()
    step.save()
    return step
