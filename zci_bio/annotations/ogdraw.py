import os.path
import re
from zipfile import ZipFile
from step_project.common.images.steps import ImagesStep
from common_utils.misc import split_list
from common_utils.file_utils import write_str_in_file, write_yaml, read_yaml, extract_from_zip

_re_zip_jpg = re.compile('GeSeqJob-[0-9]*-[0-9]*_(.*)_OGDRAW.jpg')

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


def create_ogdraw(step_data, image_format, annotations_step, common_db, sequences=None):
    step = ImagesStep(annotations_step.project, step_data, remove_data=True)
    all_images = sorted(sequences.split(';') if sequences else annotations_step.all_sequences())

    # Fetch common db sequences
    to_fetch = step.get_common_db_records(common_db, all_images, info=True)

    # If OGDraw is done on GeSeq data, than jpg images are already in
    if image_format == 'jpg':
        # Extract jpg files job-results-<num>/GeSeqJob-<num>-<num>_<seq_ident>_OGDRAW.jpg
        for filename in annotations_step.step_files(matches='^job-results-[0-9]*.zip'):
            with ZipFile(annotations_step.step_file(filename), 'r') as zip_f:
                for z_i in zip_f.infolist():
                    m = _re_zip_jpg.search(z_i.filename)
                    if m:
                        seq_ident = m.group(1)
                        if seq_ident in to_fetch:
                            to_fetch.remove(seq_ident)
                            extract_from_zip(zip_f, z_i.filename, step.step_file(seq_ident + '.jpg'))

    # Store sequence
    if to_fetch:
        # Note: it is important that file has extension gbff (multiple sequence data)
        sequences = dict()
        for i, d in enumerate(split_list(to_fetch, 30)):
            annotations_step.concatenate_seqs_genbank(step.step_file(f'sequences_{i + 1}.gbff'), d)
            sequences[i + 1] = d

        # Store instructions
        write_str_in_file(step.step_file('INSTRUCTIONS.txt'),
                          _instructions.format(step_name=step_data['step_name'], image_format=image_format))
        # Store image format used
        write_yaml(dict(image_format=image_format, sequences=sequences), step.step_file('finish.yml'))

    #
    step.set_images(all_images)
    step.save(completed=not to_fetch)
    return step


def finish_ogdraw(step_obj, common_db):
    # Note: original files are left in directory

    # Check files ogdraw-result-<num>-<hash>.zip
    zip_files = step_obj.step_files(matches='^ogdraw-result-[0-9]+-.*.zip')
    if not zip_files:
        print("Warning: can't find any OGDraw output file (ogdraw-result-*.zip)!")
        return

    # Collect sequence idents submited
    d = read_yaml(step.step_file('finish.yml'))
    image_format = d['image_format']

    seq_ident_map = dict()  # (sequence file idx, line idx) -> seq_ident
    for seq_idx, sequences in d['sequences'].items():
        # Note: line idx starts from 1, since files in zip has that numbering
        seq_ident_map.update(((seq_idx, i + 1), seq_ident) for i, seq_ident in enumerate(sequences))

    # extract ogdraw-result-<num>-<hash>/sequences_<num>ff_<num>/ogdraw_job_<hash>-outfile.<image_format>
    # Zip subdirectory naming depends on naming of OGDraw input files (sequences_<num>.gbff)
    f_end = f'-outfile.{image_format}'
    added_images = []
    for filename in zip_files:
        with ZipFile(step_obj.step_file(filename), 'r') as zip_f:
            for z_i in zip_f.infolist():
                if z_i.filename.endswith(f_end):
                    # Find sequence id of that file
                    rest = z_i.filename.split('sequences_')[1]
                    nums = re.findall(r'\d+', rest)
                    file_idx = int(nums[0])
                    line_idx = int(nums[1])
                    seq_ident = seq_ident_map[(file_idx, line_idx)]
                    #
                    added_images.append(seq_ident)
                    extract_from_zip(zip_f, z_i.filename, step_obj.step_file(f'{seq_ident}.{image_format}'))

    step_obj._check_data()
    step_obj.save(create=False)

    # Set into the common db
    if common_db:
        for image_ident in added_images:
            common_db.set_record(image_ident, step_obj.step_file(f'{image_ident}.{image_format}'))


# Extract into CommonDB
def extract_jpg_to_common_db(annotations_step, common_db):
    # Extract jpg files job-results-<num>/GeSeqJob-<num>-<num>_<seq_ident>_OGDRAW.jpg
    for filename in annotations_step.step_files(matches='^job-results-[0-9]*.zip'):
        with ZipFile(annotations_step.step_file(filename), 'r') as zip_f:
            for z_i in zip_f.infolist():
                m = _re_zip_jpg.search(z_i.filename)
                if m:
                    seq_ident = m.group(1)
                    print(seq_ident)
                    common_db.set_record_from_stream(seq_ident, zip_f.read(z_i.filename), seq_ident + '.jpg')
