from functools import wraps

# Object method cash decorators
_CACHE_PREFIX = '_cache_'


def cache(func):
    # Caches result of method that doesn't receive arguments
    # Cached value is stored in object attr named _cache_<func_name>
    c_attr = _CACHE_PREFIX + func.__name__

    @wraps(func)
    def func_wrapper(self):
        if not hasattr(self, c_attr):
            v = func(self)
            setattr(self, c_attr, v)
            return v
        return getattr(self, c_attr)
    return func_wrapper


def cache_args(func):
    # Caches result of method that receives argument(s)
    # Cached values are stored in object attr named _cache_<func_name> of type dict
    c_attr = _CACHE_PREFIX + func.__name__

    @wraps(func)
    def func_wrapper(self, *args):
        d = getattr(self, c_attr, None)
        if d is None:
            d = dict()
            setattr(self, c_attr, d)
        #
        if args not in d:
            v = func(self, *args)
            d[args] = v
            return v
        return d[args]
    return func_wrapper


def cache_remove_all(self):
    # Remove cache
    for attr in dir(self):
        if attr.startswith(_CACHE_PREFIX):
            delattr(self, attr)


def is_cached(self, method_name):
    c_attr = _CACHE_PREFIX + method_name
    return hasattr(self, c_attr)


def cache_remove(self, method_name):
    c_attr = _CACHE_PREFIX + method_name
    if hasattr(self, c_attr):
        delattr(self, attr)
