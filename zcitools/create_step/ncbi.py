from ..steps.sequences import SequencesStep
from ..utils.import_methods import import_bio_entrez


def download_ncbi(step_data, table_step, force_download=False):
    Entrez = import_bio_entrez()

    step = SequencesStep(step_data, update_mode=True)
    for ni in table_step.get_column_by_type('ncbi_ident'):
        if force_download or not step.sequence_exists(ni):
            print(f"  fetching '{ni}' in GenBank format.")
            ni_filename = ni + '.gb'
            # ToDo: what with db='nuccore'?
            filename = step.step_file(ni_filename)  # Fetch genbank format
            with Entrez.efetch(db='nuccore',  id=ni, rettype='gb', retmode='text') as handle:
                with open(filename, 'w') as out:
                    out.write(handle.read())
            step.add_sequence_file(ni_filename)
        else:
            print(f"  data for '{ni}' already presented.")

    # ToDo: remove not referenced sequences

    # Store step data
    step._check_data()
    step.save()
    return step
