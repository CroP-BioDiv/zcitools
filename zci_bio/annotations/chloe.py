import time
from .steps import AnnotationsStep
from common_utils.net_utils import download_url, get_url_json

# Protocol:
#  - Start annotation:
#    https://chloe.plantenergy.edu.au/annotate-ncbi?ncid=<NCBI_acc_num>&force_circular=true
#  - Check annotation status:
#    https://chloe.plantenergy.edu.au/task/<task_id>
#
# Status can be pending, progress and success.
# Returned data can be of these formats:
#  {"status":"PENDING",
#   "task_id":"5e40e05aff01bd341ea7578e4157b508",
#   "task_url":"/task/5e40e05aff01bd341ea7578e4157b508"}
#  {'status': 'PROGRESS',
#   'msg': '[NC_050991.1]± aligned NC_002202 (3836,781) 0.845s',
#   'task_id': '1aba6be29c8e4fe303e90b0817697b2f',
#   'task_url': '/task/1aba6be29c8e4fe303e90b0817697b2f'}
#  {"status":"SUCCESS",
#   "download_gbff":"/download/e8c40bf6a617a4fff79d93b2078b9ae8.gbff",
#   "download_gff":"/download/e8c40bf6a617a4fff79d93b2078b9ae8.gff3",
#   "download_sff":"/download/e8c40bf6a617a4fff79d93b2078b9ae8.sff",
#   "genbank_svg":"/compare-gb/e8c40bf6a617a4fff79d93b2078b9ae8.svg",
#   "ncid":"NC_048512.1",
#   "reference_svg":"/genbank-view/e8c40bf6a617a4fff79d93b2078b9ae8.svg",
#   "task_id":"e8c40bf6a617a4fff79d93b2078b9ae8",
#   "view_svg":"/annotate-view/e8c40bf6a617a4fff79d93b2078b9ae8.svg"}


def chloe_annotate(step_data, input_step, common_db, column_name=None, max_progress=8, max_tries=3):
    step = AnnotationsStep(input_step.project, step_data, remove_data=True)
    input_step.propagate_step_name_prefix(step)

    if input_step.get_data_type() == 'table':
        all_sequences = sorted(input_step.get_column_values_by_type('seq_ident', column_name=column_name))
    else:
        all_sequences = sorted(input_step.all_sequences())

    # Fetch common DB sequences
    to_fetch = step.get_common_db_records(common_db, all_sequences, info=True)

    # Store sequence
    for seq_ident in to_fetch:
        num_tries = 0
        while num_tries < max_tries:
            print(f'  Chloe: {seq_ident}', end='', flush=True)
            if data := get_url_json(f'https://chloe.plantenergy.edu.au/annotate-ncbi?ncid={seq_ident}&force_circular=true'):
                status = data['status'].upper()
                num_progress = 0
                while status in ('PENDING', 'PROGRESS') and num_progress < max_progress:
                    print('.', end='', flush=True)
                    time.sleep(3)
                    if not (data := get_url_json('https://chloe.plantenergy.edu.au' + data['task_url'])):
                        continue
                    status = data['status'].upper()
                    num_progress += 1

                if status != 'SUCCESS':
                    num_tries += 1
                    print()
                    continue

                print('#', end='', flush=True)
                gb_filename = step.step_file(f'{seq_ident}.gb')
                download_url('https://chloe.plantenergy.edu.au' + data['download_gbff'], gb_filename)
                #
                if common_db:
                    common_db.set_record(seq_ident, gb_filename)
                break

        print()

    #
    step.set_sequences(all_sequences)
    step.save()
    return step
