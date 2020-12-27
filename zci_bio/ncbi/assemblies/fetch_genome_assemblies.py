import os.path
import tarfile
import re
import csv
from decimal import Decimal
from datetime import date
from step_project.common.table.steps import TableStep, Rows2Table
from common_utils.file_utils import write_str_in_file, write_csv
from common_utils.value_data_types import fromisoformat
from ...utils.ncbi_taxonomy import get_ncbi_taxonomy

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


def _to_date(d):
    if d:
        return date(*map(int, d.split('/')))


def fetch_chloroplast_list(project, step_data, args):
    # Create table step data
    step = TableStep(project, step_data, remove_data=True)
    max_taxid = None
    if args.family:
        # step.set_step_name_prefix(args.family)  # No need for this
        max_taxid = get_ncbi_taxonomy().name_2_taxid(args.family)

    if args.csv_filename and os.path.isfile(args.csv_filename):
        column_descs = [
            dict(column='Organism Name', ),
            dict(column='BioSample', optional=True),
            dict(column='BioProject'),
            # ToDo: 1000000 or 1024*1024?
            dict(column='Size(Mb)', output='length', optional=True, tranfer=lambda x: int(float(x) * 1000000)),
            dict(column='GC%', optional=True, type='decimal'),
            dict(column='Type', optional=True, check=lambda x: x == 'chloroplast'),
            dict(column='Replicons', output='ncbi_ident', type='seq_ident', transfer=_ncbi_ident),
            dict(column='Replicons', output='other_ncbi_ident', type='seq_ident', transfer=lambda x: x.split('/')[1]),
            dict(column='CDS', optional=True, type='int'),
            dict(column='Release Date', optional=True, type='date'),
            dict(output='max_taxid', type='int', value=max_taxid),
        ]

        with open(args.csv_filename, 'r') as incsv:
            reader = csv.reader(incsv, delimiter=',', quotechar='"')
            header = next(reader)
            header[0] = header[0][1:]
            rows2table = Rows2Table(column_descs, column_names=header)
            rows2table.set_rows(reader)
            rows2table.in_table_step(step)

    #
    elif args.family:
        from ...utils.entrez import Entrez
        data = Entrez().search_summary(
            'nucleotide',
            term=f'"{args.family}"[Organism] AND ("complete genome"[Title] AND chloroplast[Title]) AND refseq')
        summary_data = dict()

        columns = [('AccessionId', 'int'), ('ncbi_ident', 'seq_ident'),
                   ('tax_id', 'int'), ('length', 'int'),
                   ('create_date', 'date'), ('update_date', 'date'),
                   ('title', 'str'), ('max_taxid', 'int')]
        family_rows = _filter_summary_data(data, max_taxid)
        _set_summary_data_for_genomes(family_rows, 'ncbi_family', summary_data)

        # Filter genomes from same species
        taxid_2_idx = dict((r[2], idx) for idx, r in enumerate(family_rows))
        species, without_sp = get_ncbi_taxonomy().group_taxids_in_species(list(taxid_2_idx.keys()))
        same_sp_rows = []
        remove_idxs = []
        for (taxid, sp_name), l_sps in species.items():
            if len(l_sps) > 1:
                # Sort by (CreateDate, ncbi_ident) desc
                sp_rows = sorted((family_rows[taxid_2_idx[taxid]] for taxid, _, _ in l_sps),
                                 reverse=True, key=lambda r: (r[4], r[0]))
                same_sp_rows.extend(([taxid, sp_name] + r) for r in sp_rows)
                remove_idxs.extend(taxid_2_idx[r[2]] for r in sp_rows[1:])
                #
                print(f'Genomes for species {sp_name} ({taxid}):')
                for row in sp_rows:
                    print(f' - {row[1]} (length: {row[3]}, date: {row[4]}) {row[6]}')

        _set_summary_data_for_genomes([family_rows[i] for i in remove_idxs], 'removed', summary_data)
        if remove_idxs:
            for idx in sorted(remove_idxs, reverse=True):
                family_rows.pop(idx)
            write_csv(step.step_file('same_species.csv'),
                      [('sp_tax_id', 'int'), ('sp_name', 'str')] + columns, same_sp_rows)

        # Add outgroup(s)
        outgroup_rows = []
        if outgroups := args.outgroup:
            data = Entrez().search_summary('nucleotide', term=' OR '.join(f'{o}[Accession]' for o in outgroups))
            if data:
                outgroup_rows = _filter_summary_data(data, 1)
                fetched = [d['Caption'] for d in data]
                if not_in := [o for o in outgroups if o not in fetched]:
                    print('Warning: no data fetched for outgroups:', not_in)
            else:
                print('Warning: no data fetched for specified outgroups!', outgroups)

        #
        all_rows = sorted(family_rows + outgroup_rows, key=lambda r: r[1])
        step.set_table_data(all_rows, columns)
        #
        _set_summary_data_for_genomes(family_rows, 'family', summary_data)
        _set_summary_data_for_genomes(outgroup_rows, 'outgroup', summary_data)
        _set_summary_data_for_genomes(all_rows, 'all', summary_data)
        step.save_summary_data(summary_data)

    step.save()  # Takes a care about complete status!
    return step


def _set_summary_data_for_genomes(rows, desc, summary_data):
    summary_data[f'{desc}_num_genomes'] = len(rows)
    summary_data[f'{desc}_max_length'] = max((r[3] for r in rows), default=None)
    summary_data[f'{desc}_min_length'] = min((r[3] for r in rows), default=None)
    summary_data[f'{desc}_max_create_date'] = max((r[4] for r in rows), default=None)
    summary_data[f'{desc}_min_create_date'] = min((r[4] for r in rows), default=None)


def _filter_summary_data(data, max_taxid):
    # ToDo: ?
    #   'Status': 'live',
    #   'ReplacedBy': '',
    #   'Flags': IntegerElement(768, attributes={}),
    #   'Comment': '  ',
    #   'Extra': 'gi|334702303|ref|NC_015543.1||gnl|NCBI_GENOMES|27514[334702303]',
    #   'AccessionVersion': 'NC_015543.1'},
    return [[int(d['Gi']), d['Caption'],
             int(d['TaxId']), int(d['Length']),
             _to_date(d['CreateDate']), _to_date(d['UpdateDate']),
             d['Title'], max_taxid]
            for d in data]
