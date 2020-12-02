# Module implements additions on BioPython's location objects (FeatureLocation and CompoundLocation)
# Note: indexing is 0-based, same as in BioPython's location objects.
#       0 <= start < end
from collections import namedtuple, defaultdict
from .helpers import read_sequence, feature_qualifiers_to_desc


class _Interval(namedtuple('_Interval', 'start, end')):
    # Strictly positive oriented interval
    def __new__(cls, start, end):
        assert 0 <= start < end, (start, end)
        return super(_Interval, cls).__new__(cls, start, end)

    def intersects(self, b):
        return (self.start <= b.start < self.end) or \
            (self.start < b.end <= self.end) or \
            (b.start <= self.start < b.end)


class Feature:
    # Location, list of intervals
    # Biopyhton location object ((FeatureLocation and CompoundLocation) or similar data
    def __init__(self, seq_length, name=None, feature=None, intervals=None, interval=None):
        self.seq_length = seq_length
        self.name = name
        if not name and feature:
            self.name = feature_qualifiers_to_desc(feature)
        self.feature = feature

        # Chec
        assert sum(int(bool(x)) for x in (feature, intervals, interval)) == 1, (feature, intervals, interval)
        if feature:
            intervals = [_Interval(p.start, p.end) for p in feature.location.parts]
        elif intervals:
            intervals = [_Interval(s, e) for s, e in intervals]
        else:
            start, end = interval
            start %= seq_length
            if end != seq_length:
                end = (end % seq_length)
            if start < end:
                intervals = [_Interval(start, end)]
            elif end == 0:
                intervals = [_Interval(start, seq_length)]
            else:
                intervals = [_Interval(start, seq_length), _Interval(0, end)]

        assert all(p.end <= seq_length for p in intervals), (seq_length, intervals)
        if len(intervals) > 1:
            # Sorte them in a way that max gap is between last and first
            intervals = sorted(intervals)
            gaps = [(b.start - a.end) % seq_length for a, b in zip(intervals[-1:] + intervals[:-1], intervals)]
            min_index = max(range(len(gaps)), key=gaps.__getitem__)
            if min_index > 0:
                intervals = intervals[min_index:] + intervals[:min_index]
        #
        self._intervals = intervals

        self._wraps = any(p.start == 0 for p in intervals) and any(p.end == seq_length for p in intervals)
        #
        self.real_start = intervals[0].start
        self.real_end = intervals[-1].end
        self.simple = len(self._intervals) == 1 or (len(self._intervals) == 2 and self._wraps)

    strand = property(lambda self: self.feature.location.strand)

    def __lt__(self, b):
        # Used for sorting
        return self.real_start < b.real_start

    def __len__(self):
        return sum((p.end - p.start) for p in self._intervals)

    def intersects(self, b):
        return any(any(p1.intersects(p2) for p2 in b._intervals) for p1 in self._intervals)

    def extract(self, seq_record):
        # Returns SeqRecord object
        assert len(seq_record) == self.seq_length, (len(seq_record), self.seq_length)
        if self.feature:
            return self.feature.extract(seq_record)
        seq = seq_record.seq
        dna = ''.join(str(seq)[p.start:p.end] for p in self._intervals)
        seq = type(seq)(dna)          # Omitting import of Bio.Seq
        return type(seq_record)(seq)  # Omitting import of Bio.SeqRecord

    def ends(self):
        return self.real_start, self.real_end


class Partition:
    def __init__(self, parts, fill=False):
        # Parts is list of Feature objects
        assert all(p.simple for p in parts), [(p.seq_length, p._intervals) for p in parts if not p.simple]
        parts = sorted(parts)

        if fill:
            new_parts = parts[:1]
            for p in parts[1:]:
                if new_parts[-1].real_end != p.real_start:
                    new_parts.append(Feature(p.seq_length, interval=(new_parts[-1].real_end, p.real_start)))
                new_parts.append(p)
            if new_parts[-1].real_end != parts[0].real_start:
                new_parts.append(Feature(p.seq_length, interval=(new_parts[-1].real_end, parts[0].real_start)))
            #
            parts = new_parts
        self._parts = parts

    def get_part_names(self):
        return (p.name for p in self._parts)

    def get_part_by_name(self, name):
        for p in self._parts:
            if p.name == name:
                return p

    __getitem__ = get_part_by_name

    def not_named_parts(self):
        return [p for p in self._parts if not p.name]

    def put_features_in_parts(self, features):
        # Returns dict partition_name -> list of features
        assert all(p.name for p in self._parts), 'All parts have to have a name!'
        ret = defaultdict(list)
        for f in features:
            in_ps = [p for p in self._parts if p.intersects(f)]
            l_ps = len(in_ps)
            assert l_ps, f.name
            if l_ps == 1:
                ret[in_ps[0].name].append(f)
            elif l_ps == 2:
                ret[' '.join(p.name for p in in_ps)].append(f)
            else:
                print(f'  warning: feature {f.name} in more than two parts!', [ff.name for ff in in_ps])
                ret['_more_'].append(f)
        return ret

    def extract_part(self, name, seq_record):
        part = self.get_part_by_name(name)
        if part:
            return part.extract(seq_record)

    def extract(self, seq_record):
        assert all(p.name for p in self._parts), 'All parts have to have a name!'
        return dict((p.name, p.extract(seq_record)) for p in self._parts)


# Helper methods
def find_features_stat(seq_rec, _type):
    # Returns dict with attrs:
    #  - without_location : number of features without location
    #  - annotated   : number of annotated features
    #  - without_name     : number of features without name (from annotated features)
    #  - disjunct    : number of disjunct features
    #  - name_strand : number different combinations of (name, strand) in disjunct features
    #  - name        : number different combinations names in disjunct features
    _fs = [f for f in seq_rec.features if f.type == _type]
    annotated = [f for f in _fs if f.location]
    disjunct = find_disjunct_features(len(seq_rec), annotated)
    name_strand = set((f.name, f.strand) for f in disjunct)
    names = set(f[0] for f in name_strand)
    return dict(
        without_location=len(_fs) - len(annotated),
        annotated=len(annotated),
        without_name=sum(1 for f in annotated if not feature_qualifiers_to_desc(f, do_assert=True)),
        disjunct=len(disjunct),
        name_strand=len(name_strand),
        names=len(names))


def _find_disjunct_features(seq_len, features):
    # Return list of features that do not overlap.
    # Assumes that features should not overlap. They can overlap because:
    #  - genes annotated with more tools (GeSeq)
    #  - wrongly annotated genes (GeSeq, rps12)

    disjunct = []
    grouped = defaultdict(list)
    for g in sorted((Feature(seq_len, feature=f) for f in features), key=lambda x: x.real_start):
        if not disjunct:
            disjunct.append(g)
        elif not disjunct[-1].intersects(g):
            if g.real_start < g.real_end or not disjunct[0].intersects(g):
                disjunct.append(g)
        else:
            # We assume that features do not overlap!
            if g.name != disjunct[-1].name:
                rem_f = disjunct.pop()
                # Find a substitute for feature we removed
                for ss in grouped[g.name]:
                    if ss.intersects(rem_f) and \
                       not ss.intersects(g) and \
                       all(not x.intersects(ss) for x in disjunct):
                        pass
                disjunct.append(g)

        grouped[g.name].append(g)
    return disjunct, grouped


def find_disjunct_features(seq_len, features):
    # Return list of features that do not overlap.
    return _find_disjunct_features(seq_len, features)[0]


def find_disjunct_features_of_type(seq_rec, _type):
    # Return list of features, so that features with same name are disjunct.
    # Note: features with different name can overlap!
    return find_disjunct_features(len(seq_rec), (f for f in seq_rec.features if f.type == _type and f.location))


def find_disjunct_genes(seq_rec):
    return find_disjunct_features_of_type(seq_rec, 'gene')


#
def split_features_in_uniq_dupl(seq_len, features):
    # Split features in two dicts(name -> list of Feature objects).
    # First contains non overlapping features, second features that overlap with features from the first one.
    _disjunct, grouped = _find_disjunct_features(seq_len, features)
    disjunct = defaultdict(list)
    for f in _disjunct:
        disjunct[f.name].append(f)
    #
    uniq = dict()
    duplicated = dict()
    for name, fs in grouped.items():
        if (un := disjunct.get(name, [])):
            uniq[name] = un
        if fs and (dup := [f for f in fs if f not in un]):
            duplicated[name] = dup

    return uniq, duplicated


def check_annotations(files_or_dirs, filter_type, output_filename):
    from .helpers import read_sequence, feature_qualifiers_to_desc
    from common_utils.file_utils import files_from_args
    from common_utils.show import print_table
    from common_utils.value_data_types import rows_2_excel

    checks = []
    for f_name in files_from_args(files_or_dirs, '.gb'):
        print(f'Processing file: {f_name}')
        seq_rec = read_sequence(f_name)

        # Extract features
        if filter_type:
            filter_type = [t.lower() for t in filter_type]
            features = [f for f in seq_rec.features if f.type.lower() in filter_type]
        else:
            features = seq_rec.features

        #
        without_location = [f for f in features if not f.location]
        if without_location:  # Fix if needed
            features = [f for f in features if f.location]
        without_name = [f for f in features if not feature_qualifiers_to_desc(f, do_assert=True)]
        uniq, duplicated = split_features_in_uniq_dupl(len(seq_rec), features)

        checks.append(dict(
            filename=f_name,
            num_without_location=len(without_location),
            without_location=', '.join(sorted(feature_qualifiers_to_desc(f) for f in without_location)),
            num_without_name=len(without_name),
            num_unique_names=len(uniq),
            num_unique_features=sum(len(_l) for _l in uniq.values()),
            num_duplicated_names=len(duplicated),
            num_duplicated_features=sum(len(_l) for _l in duplicated.values()),
        ))

    # Print and save
    print_table(['File', 'No loc', 'No name', 'Uniq names', 'Uniq fs', 'Dup. names', 'Dup. fs'],
                [[c['filename'], c['num_without_location'], c['num_without_name'],
                 c['num_unique_names'], c['num_unique_features'],
                 c['num_duplicated_names'], c['num_duplicated_features']] for c in checks])

    print(f"""
Statistics, num seqeunces with features that:
 - do not have location  : {sum(1 for c in checks if c['num_without_location'])}
 - do not have name      : {sum(1 for c in checks if c['num_without_name'])}
 - have duplicated names : {sum(1 for c in checks if c['num_duplicated_names'])}
""")
    if output_filename:
        # ToDo: in more formats?
        rows_2_excel(output_filename,
                     ['Filename', 'Num without location', 'Without location', 'Num without name',
                      'Num uniq names', 'Num uniq features',
                      'Num duplicated names', 'Num duplicated features'],
                     [[c['filename'], c['num_without_location'], c['without_location'], c['num_without_name'],
                      c['num_unique_names'], c['num_unique_features'],
                      c['num_duplicated_names'], c['num_duplicated_features']] for c in checks])
