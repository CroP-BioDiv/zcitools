import xml.etree.ElementTree as ET


class XmlDict:
    def __init__(self, node):
        self.node = node
        self._children = dict((child.tag, child) for child in node)
        if len(self._children) != len(list(node)):
            print(node.tag, [c.tag for c in node])
            assert False, 'More children with same tag!!!'

    tag = property(lambda self: self.node.tag)
    attrib = property(lambda self: self.node.attrib)
    text = property(lambda self: self.node.text)
    tail = property(lambda self: self.node.tail)

    def __getitem__(self, attr):
        return XmlDict(self._children[attr])

    def get(self, attr):
        n = self._children.get(attr)
        if n:
            return XmlDict(n)

    def __iter__(self):
        return (XmlDict(child) for child in self.node)

    @classmethod
    def fromstring(cls, s):
        return cls(ET.fromstring(f'<data>{s}</data>'))

    @classmethod
    def fromstring_nodes(cls, s):
        return [cls(n) for n in ET.fromstring(f'<data>{s}</data>')]
