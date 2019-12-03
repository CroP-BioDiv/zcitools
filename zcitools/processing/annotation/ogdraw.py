import os.path
import re
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
 * check {image_format}

Actions
 * Submit

When job is finished:
 - Download all results as zip (small disk icon in Results header) into {step_name}

When all files are processed:
 - run zcit command: zcit.py ogdraw {step_name}

-----
FAQ: https://chlorobox.mpimp-golm.mpg.de/OGDraw-FAQ.html
"""


def calculate_ogdraw(step_data, image_format, annotations_step, cache):
    step = ImagesStep(step_data, remove_data=True)
    all_images = sorted(annotations_step.all_sequences())

    # Fetch cached sequences
    to_fetch = step.get_cached_records(cache, all_images, info=True)

    # Store sequence
    if to_fetch:
        # Note: it is important that file has extension gbff (multiple sequence data)
        for i, d in enumerate(split_list(to_fetch, 30)):
            annotations_step.concatenate_seqs_genbank(step.step_file(f'sequences_{i + 1}.gbff'), d)
            # Write sequences for finish command
            with open(step.step_file(f'list_sequences_{i + 1}.txt'), 'w') as r:
                r.write('\n'.join(d))

        # Store instructions
        with open(step.step_file('INSTRUCTIONS.txt'), 'w') as out:
            out.write(_instructions.format(step_name=step_data['step_name'], image_format=image_format))

    #
    step.set_images(all_images)
    step.save(needs_editing=True)
    return step


def finish_ogdraw(step_obj, cache):
    # Note: original files are left in directory
    # Collect sequence idents submited
    seq_ident_map = dict()  # (sequence file idx, line idx) -> seq_ident
    for f in step.step_files(matches=r'^list_sequences_\d+.txt'):
        seq_idx = int(re.findall(r'\d+', f)[0])
        # Note: line idx starts from 1, since files in zip has that numbering
        seq_ident_map.update(((seq_idx, i + 1), seq_ident) for i, seq_ident in enumerate(open(self.step_file(f))))

    # Check files ogdraw-result-<num>-<hash>.zip
    # ogdraw-result-<num>-<hash>/sequences_<num>ff_<num>/ogdraw_job_<hash>-outfile.pdf
    # ogdraw-result-2019122110023-C6bn53/sequences_1ff_4/ogdraw_job_304c33c65de489abed74dfdf80f2b9fd-outfile.ps
    added_images = []
    for filename in step_obj.step_files(matches='^ogdraw-result-[0-9]-.*.zip'):
        with ZipFile(step_obj.step_file(filename), 'r') as zip_f:
            files_in_zip = defaultdict(list)  # (sequence file idx, line idx) -> list of image files
            # Collect what is in the zip
            for z_i in zip_f.infolist():



                m = _re_zip_genbank.search(z_i.filename)
                if m:
                    seq_ident = m.group(1)
                    added_seqs.append(seq_ident)
                    with open(step_obj.step_file(seq_ident + '.gb'), 'wb') as save_f:
                        save_f.write(zip_f.read(z_i.filename))
    # create=False
    pass
