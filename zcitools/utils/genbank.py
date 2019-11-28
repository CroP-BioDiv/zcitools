"GenBank format helper methods"


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
