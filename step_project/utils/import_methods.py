from common_utils.import_method import import_method

_missing_pygraphviz = """
PyGraphviz library is missing.

For installation instruction check web page:
https://pygraphviz.github.io/documentation/latest/install.html

Short: pip install pygraphviz
"""


@import_method(_missing_pygraphviz)
def import_pygraphviz():
    import pygraphviz
    return pygraphviz
