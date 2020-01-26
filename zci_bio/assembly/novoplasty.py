import os.path
import subprocess
from .steps import AssemblyStep
from ..utils.sequence_reads import SequenceReads
from common_utils.file_utils import write_str_in_file

_config_data = """Project:
-----------------------
Project name          = {project_name}
Type                  = {organelle_type}
Genome Range          = {genome_range_from}-{genome_range_to}
K-mer                 = {k_mer}
Max memory            =
Extended log          =
Save assembled reads  =
Seed Input            = {seed}
Reference sequence    =
Variance detection    =
Chloroplast sequence  =

{datasets}

Heteroplasmy:
-----------------------
MAF                   =
HP exclude list       =
PCR-free              =

Optional:
-----------------------
Insert size auto      = yes
Insert Range          = 1.9
Insert Range strict   = 1.3
Use Quality Scores    =
"""

_dataset_data = """Dataset {num}:
-----------------------
Read Length           = {read_length}
Insert size           = {insert_size}
Platform              = {platform}
Single/Paired         = {single_paired}
Combined reads        =
Forward reads         = {read_1}
Reverse reads         = {read_2}"""

_organele_types = dict(c='chloro', m='mito')
_organele_range = dict(c=(120000, 180000), m=(12000, 21000))
_config_file = 'novoplasty_config.txt'


def novoplasty_run(project, step_data, params):
    step = AssemblyStep(project, step_data, remove_data=True)

    # Create config file
    sr = SequenceReads.from_file(params.sequence_reads, relative_dir='..')
    datasets = [_dataset_data.format(
        num=num + 1,
        read_length=read.get('read_length', ''),
        insert_size=read.get('insert_length', ''),
        platform=read.get('platform', '').lower() or 'illumina',  # illumina or ion
        read_1=read.get('file' if rt == 'S' else 'file_1'),
        read_2=read.get('file_2', ''),
        single_paired=rt) for num, (rt, read) in enumerate(sr)]

    ot = params.organelle_type[0].lower()
    g_range = _organele_range[ot]
    conf = _config_data.format(
        project_name=params.project_name or step_data['step_name'],
        organelle_type=_organele_types[ot],
        genome_range_from=params.genome_range_from or g_range[0],
        genome_range_to=params.genome_range_to or g_range[1],
        k_mer=params.k_mer,
        seed=os.path.join('..', params.seed),
        datasets='\n\n'.join(datasets),
    )
    write_str_in_file(step.step_file(_config_file), conf)

    # Run NOVOPlasty
    # ToDo: napraviti opcenitije, run ovdje ili server, trazenje exe-a, ...
    print(f"Command: cd {step.directory}; NOVOPlasty -c {_config_file} > /dev/null")
    subprocess.run(['NOVOPlasty', '-c', _config_file], cwd=step.directory, stdout=subprocess.DEVNULL)

    step.save()
    return step
