
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
