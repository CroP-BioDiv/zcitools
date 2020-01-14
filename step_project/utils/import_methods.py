from common_utils.import_method import import_method


@import_method("""
PyGraphviz library is missing.

For installation instruction check web page:
https://pygraphviz.github.io/documentation/latest/install.html

Short: pip install pygraphviz
""")
def import_pygraphviz():
    import pygraphviz
    return pygraphviz
