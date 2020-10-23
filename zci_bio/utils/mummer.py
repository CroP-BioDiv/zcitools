from types import SimpleNamespace


class MummerDelta:
    def __init__(self, delta_filename):
        self._matches = []

        with open(delta_filename, 'r') as _in:
            self.base_filename, self.query_filename = next(_in).strip().split()
            line = next(_in).strip()
            assert line == 'NUCMER', line
            for line in _in:
                line = line.strip()
                if line[0] == '>':
                    sequence, query, sequence_len, query_len = line[1:].split()
                else:
                    fields = line.split()
                    if len(fields) >= 4:
                        qi = (int(fields[2]), int(fields[3]))
                        self._matches.append(SimpleNamespace(
                            sequence=sequence,
                            query=query,
                            sequence_interval=(int(fields[0]), int(fields[1])),
                            query_interval=qi,
                            positive=(qi[1] > qi[0]),
                        ))

    def __bool__(self):
        return bool(self._matches)

    def __len__(self):
        return len(self._matches)

    def aligns(self, s, q):
        return sorted((m for m in self._matches if (s is None or m.sequence == s) and m.query == q),
                      key=lambda x: x.positive,
                      reverse=True)
