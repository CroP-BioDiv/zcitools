import os.path
import tarfile
import datetime
import re
from decimal import Decimal
from step_project.common.table.steps import TableStep
from common_utils.file_utils import write_str_in_file

_instructions = """
Steps:
 - open NCBI Assemblies page https://www.ncbi.nlm.nih.gov/assembly
 - search for all assemblies. In search field enter: all[filter]
 - refine search with options on the left
 - click on "Dowload Assemblies" button
 - for Source database choose between GenBank (all assemblies) or RefSeq (referential assemblies)
 - for file type choose under Reports "Assembly statistics report"
 - download data into project's step directory {step_name}. File should be called genome_assemblies.tar
 - run zcit command: zcit.py finish {step_name}

Notes on search refining:
 - It is enough to take Latest. That means latest version of assembly
 - Excludes are good to keep

"""


def fetch_genome_assemblies(project, step_data):
    # Create table step data
    step = TableStep(project, step_data, remove_data=True)
    step.save()  # Takes a care about complete status!

    # Set instructions
    write_str_in_file(step.step_file('INSTRUCTIONS.txt'), _instructions.format(step_name=step.directory))

    return step


def finish_fetch_genome_assemblies(step):
    ga_file = step.step_file('genome_assemblies.tar')
    if not os.path.isfile(ga_file):
        print('Error: No genome_assemblies.tar file!!!')
        return

    if not tarfile.is_tarfile(ga_file):
        print('Error: No genome_assemblies.tar file!!!')
        return

    records = []
    with tarfile.open(ga_file) as t_file:
        for tarinfo in t_file:
            if tarinfo.isreg() and tarinfo.name.endswith('_assembly_stats.txt'):
                file_data = t_file.extractfile(tarinfo)
                records.append(_extract_data(file_data.readlines(), tarinfo.name))

    step.set_table_data(records, _columns)
    step.save()


def _extract_data(lines, file_name):
    # List to dict
    data = dict()
    for line in lines:
        line = line.decode("utf-8").strip()
        if line.startswith('# '):
            if ':' in line:
                line = line[2:]
                i = line.index(':')
                if i < len(line) - 1:
                    data[line[:i].strip()] = line[i + 1:].strip()
        elif line.startswith('all'):
            fields = line.split()
            assert len(fields) == 6, (file_name, line)
            assert all(f == 'all' for f in fields[:4]), (file_name, line)
            data[fields[4]] = int(fields[5])

    not_in = [f for f in _req_fields if f not in data]
    assert not not_in, (file_name, not_in)

    coverage = data.get('Genome coverage')
    if coverage:
        nums = re.findall(r"\d+\.\d+", coverage)
        if len(nums) > 1:
            print(f'Warning: coverage is "{coverage}"')
        coverage = Decimal(nums[0]) if nums else None

    return [
        # datetime.date.fromisoformat(data['Date']),  # In python 3.7 :-/
        datetime.date(*map(int, data['Date'].split('-'))),
        data['Assembly name'],
        data['Organism name'],
        int(data['Taxid']),
        #
        data['BioProject'],
        data.get('BioSample'),
        data['GenBank assembly accession'],
        data.get('RefSeq assembly accession'),
        data.get('WGS project'),
        #
        int(data['Genome representation'] == 'full'),
        data.get('Sequencing technology'),
        coverage,
        data.get('Assembly type', 'haploid'),  # Default
        data['Assembly level'],
        data.get('Assembly method'),
        #
        data['total-length'],
        data.get('contig-count', -1),
        data.get('contig-N50', -1),
        data.get('contig-L50', -1),
        data.get('scaffold-count', -1),
        data.get('scaffold-N50', -1),
        data.get('scaffold-L50', -1),
        data.get('scaffold-N75', -1),
        data.get('scaffold-N90', -1),
    ]


_req_fields = (
    'Assembly name', 'Organism name', 'Taxid', 'BioProject', 'Date',
    'Assembly level', 'GenBank assembly accession', 'Genome representation',
    'total-length',
    # Old assemblies with Assembly level == Complete Genome, dont have contig-*
    # 'contig-count', 'contig-N50', 'contig-L50',
)

_columns = (
    ('date', 'date'),
    ('assembly_name', 'str'),
    ('organism_name', 'str'),
    ('taxid', 'int'),
    #
    ('bio_project', 'str'),
    ('bio_sample', 'str'),
    ('gen_bank', 'str'),
    ('ref_seq', 'str'),
    ('WGS', 'str'),
    #
    ('full_genome', 'int'),
    ('sequencing_technology', 'str'),
    ('coverage', 'decimal'),
    ('assembly_type', 'str'),
    ('assembly_level', 'str'),
    ('assembly_method', 'str'),
    #
    ('total_length', 'int'),
    ('contig_count', 'int'),
    ('contig_N50', 'int'),
    ('contig_L50', 'int'),
    ('scaffold_count', 'int'),
    ('scaffold_N50', 'int'),
    ('scaffold_L50', 'int'),
    ('scaffold_N75', 'int'),
    ('scaffold_N90', 'int'),
)
