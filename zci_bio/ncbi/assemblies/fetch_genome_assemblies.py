import os.path
import tarfile
import re
import csv
from decimal import Decimal
from step_project.common.table.steps import TableStep, Rows2Table
from common_utils.file_utils import write_str_in_file
from common_utils.value_data_types import fromisoformat

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
        fromisoformat(data['Date']),
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


# ---------------------------------------------------------
# Chloroplast
# ---------------------------------------------------------
def _ncbi_ident(x):
    x = x.split('/')[0]         # NC_046760.1/MN087227.1 -> NC_046760.1
    x = x.split('.')[0]         # NC_046760.1            -> NC_046760
    if ':' in x:
        return x.split(':')[1]  # Pltd:NC_007578         -> NC_007578
    return x


def fetch_chloroplast_list(project, step_data, args):
    # Create table step data
    step = TableStep(project, step_data, remove_data=True)
    if args.family:
        step.set_step_name_prefix(args.family)

    if args.csv_filename and os.path.isfile(args.csv_filename):
        column_descs = [
            dict(column='Organism Name', ),
            dict(column='BioSample', optional=True),
            dict(column='BioProject'),
            # ToDo: 1000000 or 1024*1024?
            dict(column='Size(Mb)', output='size', optional=True, tranfer=lambda x: int(float(x) * 1000000)),
            dict(column='GC%', optional=True, type='decimal'),
            dict(column='Type', optional=True, check=lambda x: x == 'chloroplast'),
            dict(column='Replicons', output='ncbi_ident', type='seq_ident', transfer=_ncbi_ident),
            dict(column='Replicons', output='other_ncbi_ident', type='seq_ident', transfer=lambda x: x.split('/')[1]),
            dict(column='CDS', optional=True, type='int'),
            dict(column='Release Date', optional=True, type='date'),
        ]

        with open(args.csv_filename, 'r') as incsv:
            reader = csv.reader(incsv, delimiter=',', quotechar='"')
            header = next(reader)
            header[0] = header[0][1:]
            rows2table = Rows2Table(column_descs, column_names=header)
            rows2table.set_rows(reader)
            rows2table.in_table_step(step)

    elif args.family:
        pass  # ToDo:

    step.save()  # Takes a care about complete status!
    return step
