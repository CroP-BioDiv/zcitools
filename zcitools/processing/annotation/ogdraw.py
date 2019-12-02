import os.path
from zcitools.steps.images import ImagesStep
from zcitools.utils.helpers import split_list

_instructions = """
Open web page: https://chlorobox.mpimp-golm.mpg.de/OGDraw.html

For each GenBank file {step_name}/*.gbff do:

FASTA file(s) to annotate
 * Upload file
 * (check) Circular
 * (check) Plastid
 * (check) Tidy up annotation

Inverted Repeat
 * (check) Auto

Output Options
 * check one output format type (PS prefered?)

Actions
 * Submit

When job is finished:
 - Download all results as zip (small disk icon in Results header) into {step_name}

When all files are processed:
 - run zcit command: zcit.py ogdraw {step_name}

-----
FAQ: https://chlorobox.mpimp-golm.mpg.de/OGDraw-FAQ.html
"""


def calculate_ogdraw(step_data, annotations_step):
    step = ImagesStep(step_data, remove_data=True)
    # cache = step.get_cache_object()
    all_images = sorted(annotations_step.all_sequences())

    # Fetch cached sequences
    # to_fetch = step.get_cached_records(cache, all_images, info=True)
    to_fetch = all_images

    # Store sequence
    if to_fetch:
        # Note: it is important that file has extension gbff (multiple sequence data)
        for i, d in enumerate(split_list(to_fetch, 30)):
            annotations_step.concatenate_seqs_genbank(step.step_file(f'sequences_{i + 1}.gbff'), d)

        # Store instructions
        with open(step.step_file('INSTRUCTIONS.txt'), 'w') as out:
            out.write(_instructions.format(step_name=step_data['step_name']))

    #
    step.set_images(all_images)
    step.save(needs_editing=True)
    return step


def finish_ogdraw(step_obj):
    print('FIFFIIFAISFI')
    pass
