from functools import wraps

"""
Method for importing non standard python libraries.
Idea of this methods is to inform user how to import needed library if it is not installed.

Methods receive *args as arguments and return:
 - library if args is empty,
 - list of library attributes if args are given.
"""


# Wrapper for import protocol
def import_method(error_msg):
    def wrap(func):
        def wrapped_f(*args):
            try:
                lib = func()
            except ImportError:
                print()
                print('=' * 80)
                print(error_msg)
                print('=' * 80)
                print()
                raise

            if args:
                if len(args) == 1:
                    return getattr(lib, args[0])  # No unpacking!
                return [getattr(lib, a) for a in args]
            return lib

        return wrapped_f
    return wrap


# --------------------------------------------------
@import_method("""
Pandas library is missing.

For installation instruction check web page:
https://pandas.pydata.org/getpandas.html

Short: pip install pandas
""")
def import_pandas():
    import pandas
    return pandas


_matplotlib_desc = """
Matplotlib library is missing.

For installation instruction check web page:
https://matplotlib.org/users/installing.html

Short: pip install matplotlib
"""
@import_method(_matplotlib_desc)
def import_matplotlib_pyplot():
    import matplotlib.pyplot as pyplot
    return pyplot
