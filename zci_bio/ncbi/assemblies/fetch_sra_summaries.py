from collections import defaultdict
from step_project.common.table.steps import TableStep, TableGroupedStep
from common_utils.misc import split_list, YYYYMMDD_2_date, int_2_human, human_2_int, coverage_2_human
from common_utils.xml_dict import XmlDict
from common_utils.step_database import StepDatabase
from common_utils.import_method import import_pandas
from common_utils.data_types.table import HierarchicalTable
from ...utils.entrez import Entrez

_sra_columns = (
    ('bio_project', 'str'),  # Group column
    ('id', 'int'),
    ('created', 'date'),
    ('updated', 'date'),
    # ('bio_sample', 'str'),
    ('study_acc', 'str'),
    ('sample_acc', 'str'),
    ('name', 'str'),
    ('experiment_acc', 'str'),
    ('platform', 'str'),
    ('instrument_model', 'str'),
    ('library_strategy', 'str'),
    ('library_source', 'str'),
    ('library_selection', 'str'),
    ('library_paired', 'str'),
    ('total_spots', 'int'),
    ('total_bases', 'int'),
)


def fetch_sra_summaries(step_data, table_step):
    step = TableGroupedStep(table_step.project, step_data, update_mode=True)
    step.set_columns(_sra_columns)
    step.save()  # To store step as it is
    step.known_groups()
    bp_column = table_step.choose_first_column('bio_project', 'BioProject', error=True)
    #
    to_proc = table_step.get_column_values(bp_column) - set(step.known_groups())
    if to_proc:
        entrez = Entrez()
        num_grouped_sras = 2000
        for proj in to_proc:
            data = entrez.esearch(db='sra', term=f"{proj}[BioProject]", retmax=num_grouped_sras)
            ids = data['IdList']
            sra_data = []

            if not ids:
                print(f"No SRA data for project: {proj}!")
            else:
                records = entrez.esummary(db='sra', id=','.join(ids), retmax=num_grouped_sras)
                print(f"{proj}: {len(ids)} / {len(records)}")
                for record in records:
                    try:
                        sra = extract_data(record)
                        sra_data.append([sra[c] for c, _ in _sra_columns[1:]])
                    except Exception:
                        print(record)
                        raise
                    # ToDo: check
                    #  - is taxid same as in project

            step.set_group_rows(proj, sra_data)

        # Note: it is not possible to retrive data by group of projects since not all SRA summaries have BioProject!!!
        # num_grouped_projects = 8
        # for projs in split_list(sorted(to_proc), num_grouped_projects):
        #     term = ' OR '.join(f"{p}[BioProject]" for p in projs)
        #     data = entrez.esearch(db='sra', term=term, retmax=num_grouped_sras)
        #     #
        #     proj_sras = defaultdict(list)
        #     ids = data['IdList']
        #     if not ids:
        #         print(f"No SRA data for projects: {', '.join(projs)}!")
        #     else:
        #         records = entrez.esummary(db='sra', id=','.join(ids), retmax=num_grouped_sras)
        #         for record in records:
        #             try:
        #                 r_data = extract_data(record)
        #             except:
        #                 print(record)
        #                 raise
        #             proj_sras[r_data['bio_project']].append(r_data)
        #             # ToDo: check
        #             #  - is taxid same as in project

        #         for proj in projs:
        #             substep = step.create_substep(proj)
        #             sra_data = [[sra[c] for c, _ in _sra_columns] for sra in proj_sras.get(proj, [])]
        #             substep.set_table_data(sra_data, _sra_columns)
        #             substep.save()

    return step


def extract_data(record):
    exp_xml = XmlDict.fromstring(record['ExpXml'])
    summary = exp_xml['Summary']
    stat_attrs = summary['Statistics'].attrib
    platform = summary['Platform']
    exp_attrs = exp_xml['Experiment'].attrib
    #
    lib = exp_xml['Library_descriptor']
    library_paired = lib['LIBRARY_LAYOUT'].get('PAIRED')
    #
    runs = XmlDict.fromstring_nodes(record['Runs'])
    #
    return dict(
        id=int(record['Id']),
        created=YYYYMMDD_2_date(record['CreateDate']),
        updated=YYYYMMDD_2_date(record['UpdateDate']),
        #
        taxid=int(exp_xml['Organism'].attrib['taxid']),
        # Note: not all SRA summaries have BioProject!!!
        # bio_project=exp_xml['Bioproject'].text,
        # bio_sample=exp_xml['Biosample'].text,
        study_acc=exp_xml['Study'].attrib['acc'],
        sample_acc=exp_xml['Sample'].attrib['acc'],
        #
        name=exp_attrs['name'],
        experiment_acc=exp_attrs['acc'],
        platform=platform.text,
        instrument_model=platform.attrib['instrument_model'],
        #
        library_strategy=lib['LIBRARY_STRATEGY'].text,
        library_source=lib['LIBRARY_SOURCE'].text,
        library_selection=lib['LIBRARY_SELECTION'].text,
        library_paired=library_paired.attrib.get('NOMINAL_LENGTH') if library_paired else None,
        #
        total_spots=int(stat_attrs.get('total_spots') or 0),
        total_bases=int(stat_attrs.get('total_bases') or 0),
        #
        runs=[dict(accession=r.attrib['acc'],
                   spots=int(r.attrib.get('total_spots') or 0),
                   bases=int(r.attrib.get('total_bases') or 0))
              for r in runs]
    )


def group_sra_data(step_data, sra_step):
    # Creates step from this select statement
    # SELECT bio_project, instrument_model, library_paired, COUNT(*), SUM(total_spots), SUM(total_bases)
    # GROUP BY bio_project, instrument_model, library_paired
    with StepDatabase([sra_step]) as db:
        column_data_types, rows = db.select_all_tables(
            "bio_project, instrument_model, library_paired, " +
            "SUM(total_spots) as total_spots, SUM(total_bases) as total_bases, COUNT(*) as count",
            group_by_part="bio_project, instrument_model, library_paired")

    #
    step = TableStep(sra_step.project, step_data, remove_data=True)
    step.set_table_data(rows, column_data_types)
    step.save()
    return step


def make_report(assemblies, sra, params):
    assem_columns = (
        'bio_project', 'date', 'assembly_method', 'total_length', 'coverage', 'organism_name',
        'contig_count', 'contig_N50')  # , 'contig_L50'
    # 'scaffold_count', 'scaffold_N50', 'scaffold_L50', 'scaffold_N75', 'scaffold_N90')
    sra_columns = ['instrument_model', 'library_paired', 'total_spots', 'total_bases', 'count']
    #
    columns = [['BioProject', 'Date', 'Method', 'Genome length', 'Coverage', 'Organism', 'Contigs', 'C N50'],  # 'C L50'
               ['Instrument', 'Paired', 'Spots', 'Bases', 'Count']]
    data_types = [['str', 'date', 'str', 'int', 'str', 'str', 'int', 'int'],
                  ['str', 'int', 'int', 'int', 'int']]
    table = HierarchicalTable(columns, data_types)
    #
    where_part = []
    if params.method:
        where_part.append(f"UPPER(assembly_method) like '%{params.method.upper()}%'")
    else:
        where_part.append("assembly_method != ''")
    if params.from_genome_size:
        where_part.append(f"total_length >= {human_2_int(params.from_genome_size)}")
    if params.to_genome_size:
        where_part.append(f"total_length <= {human_2_int(params.to_genome_size)}")
    if params.from_date:
        where_part.append(f"date >= '{params.from_date}'")
    if params.to_date:
        where_part.append(f"date >= '{params.to_date}'")

    with StepDatabase([assemblies, sra]) as db:
        projects = db.select_result(
            f"SELECT {','.join(assem_columns)} FROM a WHERE {' AND '.join(where_part)} ORDER BY date")
        if projects:
            db.cursor.execute("CREATE INDEX xx ON b (bio_project)")  # For speed
            for x in projects:
                sras = db.select_result(f"SELECT {','.join(sra_columns)} FROM b WHERE bio_project = '{x[0]}'")
                if sras:
                    gen_l = x[3]
                    cov = coverage_2_human(sum(s[-2] for s in sras) / gen_l) if gen_l else ''
                    x = x[:4] + (f"{coverage_2_human(x[4])}/{cov}",) + x[5:]
                    table.append_row(0, x)
                    table.extend_rows(1, sras)

            # Post-format data in more readable form
            table.update_column(0, 'Genome length', int_2_human, update_data_type='str')
            table.update_column(1, 'Spots', int_2_human, update_data_type='str')
            table.update_column(1, 'Bases', int_2_human, update_data_type='str')

    if params.print:
        table.print(show_limit=7)
    else:
        table.to_excel(params.output)
