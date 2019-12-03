import os.path
from .import_methods import import_bio_seq_io

# Methods for this and that

# Biopyhton helpers
_bio_ext_2_type = dict(
    gb='genbank', gbff='genbank',
    fa='fasta',  fas='fasta',
)
_bio_ext_2_type = dict(('.' + e, t) for e, t in _bio_ext_2_type.items())


def feature_location_desc(location):
    # Returns tuple of indices for describing location as simple as possible
    # Location can be specified in lot of ways, check:
    #   https://biopython.org/DIST/docs/api/Bio.SeqFeature-module.html

    cls_name = location.__class__.__name__  # To prevent importing of classes
    if cls_name == 'FeatureLocation':
        return (int(location.start), int(location.end))

    if cls_name == 'CompoundLocation':
        parts = location.parts
        return (int(parts[0].start), int(parts[-1].end))
        # return (int(location.start), int(location.end))

    assert False, (f'Not supported type {cls_name}!', location)


def feature_qualifiers_to_desc(feature):
    # From doc:
    # qualifiers - A dictionary of qualifiers on the feature.
    #   These are analogous to the qualifiers from a GenBank feature table.
    #   The keys of the dictionary are qualifier names, the values are the qualifier values.
    #   As of Biopython 1.69 this is an ordered dictionary.

    qualifiers = feature.qualifiers
    if feature.type in ('gene', 'CDS'):
        genes = qualifiers['gene']
        assert len(genes) == 1, genes
        return genes[0]

    if feature.type == 'repeat_region':
        r_type = qualifiers['rpt_type']
        assert len(r_type) == 1, r_type
        return r_type[0]

    return str(qualifiers)  # ToDo: For now


# Sequence files
def write_fasta(filename, data):
    with open(filename, 'w') as fa:
        for ident, seq in data:
            fa.write(f">{ident}\n{seq}\n")


def split_sequences(input_filename, output_ext):
    SeqIO = import_bio_seq_io()
    input_type = _bio_ext_2_type[os.path.splitext(input_filename)[1]]
    input_dir = os.path.dirname(input_filename)
    output_type = _bio_ext_2_type[output_ext]
    sequence_ids = []
    with open(input_filename, 'r') as seqs:
        for rec in SeqIO.parse(seqs, input_type):
            out_f = rec.id + output_ext
            if input_dir:
                out_f = os.path.join(input_dir, out_f)
            SeqIO.write([rec], open(out_f, 'w'), output_type)
            sequence_ids.append(rec.id)
    return sequence_ids


def concatenate_sequences(output_filename, input_filenames):
    SeqIO = import_bio_seq_io()
    output_type = _bio_ext_2_type[os.path.splitext(output_filename)[1]]
    with open(output_filename, 'w') as out_seqs:
        for in_f in input_filenames:
            with open(in_f, 'r') as seq:
                SeqIO.write(list(SeqIO.parse(seq, _bio_ext_2_type[os.path.splitext(in_f)[1]])), out_seqs, output_type)


# Other
def split_list(data, num_items):
    assert isinstance(data, list), type(data)
    for i in range((len(data) // num_items) + 1):
        n = i * num_items
        yield data[n:(n + num_items)]
