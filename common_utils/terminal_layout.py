class Box:
    def __init__(self, direction, boxes=None):
        assert isinstance(direction, str) and direction[0].lower() in 'rc', (direction, type(direction))
        self.direction = direction[0].lower()
        self.boxes = list(boxes or [])
        self._size = None

    children = property(lambda self: self.boxes)

    def add_box(self, box):
        assert isinstance(box, Box), type(box)
        self.boxes.append(box)
        self._size = None

    def get_opposite_direction(self):
        return 'c' if self.direction == 'r' else 'r'

    def set_size_coord(self, c, v):
        size = list(self.get_size())
        size[c] = v
        self._size = tuple(size)

    def get_size(self):
        if self._size is None:
            self._size = self._get_size()
        return self._size

    def _get_size(self):
        # Returns tuple (x, y)
        if not self.boxes:
            return (0, 0)
        else:
            sizes = [b.get_size() for b in self.boxes]
            if self.direction == 'r':
                x = max(s for s, _ in sizes)
                y = sum(s for _, s in sizes)
            else:
                x = sum(s for s, _ in sizes)
                y = max(s for _, s in sizes)
            return x, y

    def get_lines(self, x=None):
        # Returns list of strings that represent box
        if not self.boxes:
            return []

        sx, sy = self.get_size()
        if self.direction == 'r':
            ret = []
            for b in self.boxes:
                x, _ = b.get_size()
                ret.extend(line.ljust(x) for line in b.get_lines(x=sx))
        else:
            ret = ['' for _ in range(sy)]
            for b in self.boxes:
                x, _ = b.get_size()
                i = 0
                for line in b.get_lines(x=sx):
                    ret[i] += line.ljust(x)
                    i += 1
                while i < len(ret):
                    ret[i] += ' ' * x
                    i += 1
        return ret

    def __str__(self):
        return '\n'.join(self.get_lines())

    def print_sizes(self, ident=''):
        print(ident + str(self.get_size()), getattr(self, 'direction', ''))
        for c in self.children:
            c.print_sizes(ident=ident + '  ')


class WrappedBox(Box):
    # Fills box in given direction until given number of items is reached, than 'goes' to next row/column
    def __init__(self, direction, boxes=None, num_items=5):
        super().__init__(direction)
        self.num_items = num_items
        self._num_current_items = 0
        self._opposite_direction = self.get_opposite_direction()
        if boxes:
            for b in boxes:
                self.add_box(b)

    def add_box(self, box):
        if self._num_current_items == 0:
            if self.boxes:
                self.boxes.append(StringBox(''))  # add separator
            self._opposite_box = Box(self._opposite_direction)
            self.boxes.append(self._opposite_box)
        self._num_current_items = (self._num_current_items + 1) % self.num_items
        self._opposite_box.add_box(box)
        self._size = None


class StringListBox(Box):
    # List of strings
    def __init__(self, lines, padding=1, stretch=False):
        self.stretch = stretch
        self._orig_lines = lines
        self.padding = padding
        #
        p = ' ' * padding
        lines = [p + l + p for l in lines]
        self._size = (max(len(l) for l in lines), len(lines)) if lines else (0, 0)
        self._lines = [l.ljust(self._size[0]) for l in lines]

    children = property(lambda self: tuple())

    def get_lines(self, x=None):
        if self.stretch:
            assert x
            p = ' ' * self.padding
            x -= 2 * self.padding
            return [p + l * (x // len(l)) + p for l in self._orig_lines]
        # if x:
        #     return [l.ljust(x) for l in self._lines]
        return self._lines


class StringBox(StringListBox):
    def __init__(self, s, padding=1, stretch=False):
        super().__init__([s], padding=padding, stretch=stretch)


class StringColumns(Box):
    # Simple table from given data
    def __init__(self, rows, header=None, max_data_length=None, padding=1):
        # Note: header element is string or tuple if header is in more columns
        if header:
            num_cols = len(header)
        else:
            num_cols = len(rows[0])
        assert all(len(r) == num_cols for r in rows), (num_cols, [r for r in rows if len(r) != num_cols])

        # Make strings, to be sure
        rows = [['' if cell is None else str(cell) for cell in r] for r in rows]
        if max_data_length:
            l3 = max_data_length - 3

            def _f(c):
                if len(c) <= max_data_length:
                    return c
                return c[:l3] + '...'

            rows = [[_f(cell) for cell in r] for r in rows]

        if header:
            if isinstance(header[0], (tuple, list)):
                columns = [StringListBox(list(h) + ['-' * len(h)] + lines, padding=padding)
                           for h, *lines in zip(header, *rows)]
            else:
                columns = [StringListBox([h, '-' * len(h)] + lines, padding=padding)
                           for h, *lines in zip(header, *rows)]
        else:
            columns = [StringListBox(lines, padding=padding) for lines in zip(*rows)]

        super().__init__('c', boxes=columns)


class TableByColumns(StringColumns):
    def __init__(self, columns, padding=1):
        super().__init__(list(zip(*columns)))


class TreeBox(Box):
    # root has interface attr children (list) and attr label (str)
    def __init__(self, root, direction='column'):
        self.direction = direction[0].lower()  # For get_opposite_direction()
        #
        boxes = [StringBox(root.label)]
        if len(root.children) == 1:
            boxes.append(TreeBox(root.children[0], direction))
        else:
            boxes.append(Box(self.get_opposite_direction(), boxes=[TreeBox(c) for c in root.children]))
        super().__init__(direction, boxes=boxes)


#
def fill_rows(rows):
    max_l = max(len(r) for r in rows)
    for r in rows:
        if len(r) < max_l:
            r.extend([''] * (max_l - len(r)))
