# Module implements additions on BioPython's location objects (FeatureLocation and CompoundLocation)
# Note: indexing is 0-based, same as in BioPython's location objects.
#       0 <= start < end, inrades in are range(start, end).
from collections import namedtuple, defaultdict


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
            self.name = feature.qualifiers['gene'][0]  # For now
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
