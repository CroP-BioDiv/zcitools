import io
import re
from ..utils.import_methods import import_bio_seq_io
from .utils import find_chloroplast_irs
from common_utils.value_data_types import rows_2_excel


def analyse_genomes_start(table_step, output_file, common_db):
    SeqIO = import_bio_seq_io()
    _length_match = re.compile('.*, length ([0-9]+)')

    # Fetch CommonDB data for variants
    dbs = ('GeSeq', 'IRs_mummer')  # 'sequences',
    c_dbs = [common_db.get_relative_db(db) for db in dbs]

    # Extract IR data
    rows = []
    for seq_ident in sorted(table_step.get_column_values_by_type('seq_ident')):
        row = [seq_ident, None]
        rows.append(row)
        for db, c_db in zip(dbs, c_dbs):
            data = c_db.get_record_data(seq_ident)
            last_num_cols = len(row)
            if data:
                stream = io.StringIO(data.decode("utf-8"))
                seq_rec = SeqIO.read(stream, 'genbank')
                if not row[1]:
                    row[1] = len(seq_rec)
                irs = find_chloroplast_irs(seq_rec)
                if irs:
                    ira, irb = irs
                    loc = ira.location
                    row.append(loc.start)
                    row.append(len(loc))
                    # Match length
                    length = None
                    note = ira.qualifiers.get('note')
                    if note:
                        m = _length_match.search(note[0])
                        if m:
                            length = int(m.group(1))
                    row.append(length)
                    # Quality (IR length in 25+-1kb, IRB ends on end, IRA starts on second half)
                    quality = []
                    if not (24000 <= len(loc) <= 26000):
                        quality.append('L')
                    if len(seq_rec) - irb.location.end > 10:
                        quality.append('B')
                    if loc.start < len(seq_rec) // 2:
                        quality.append('A')
                    row.append(''.join(quality))

            #
            if last_num_cols == len(row):
                row.extend([None] * 4)

    # Create table
    columns = ['Seq', 'Length']
    for db in ('GeSeq', 'Mummer'):  # 'NCBI',
        columns.extend(f'{db} {c}' for c in ('start', 'length', 'match', 'quality'))
    rows_2_excel(output_file, columns, rows)
