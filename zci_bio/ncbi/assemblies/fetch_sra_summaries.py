from collections import defaultdict
from step_project.common.table.steps import TableStep, TableGroupedStep
from common_utils.misc import split_list, YYYYMMDD_2_date
from common_utils.xml_dict import XmlDict
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
    #
    to_proc = table_step.get_column_values('bio_project') - set(step.known_groups())
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
    step = TableStep(sra_step.project, step_data)
    # Assumes _sra_columns
    bp_idx = sra_step._column_index('bio_project')
    platform_idx = sra_step._column_index('platform')
    instr_idx = sra_step._column_index('instrument_model')

    # print(table_step.get_rows())
    columns = (())
    rows = []
    step.set_table_data(rows, columns)
    step.save()
    return step


    # ('bio_project', 'str'),  # Group column
    # ('platform', 'str'),
    # ('instrument_model', 'str'),
    # ('library_strategy', 'str'),
    # ('library_source', 'str'),
    # ('library_selection', 'str'),
    # ('library_paired', 'str'),
    # ('total_spots', 'int'),
    # ('total_bases', 'int'),
