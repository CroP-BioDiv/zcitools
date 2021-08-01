try:
    import psutil
    _data = psutil.virtual_memory()
    virtual_memory = _data.total
except ImportError:
    # ToDo: shell command free
    import os
    virtual_memory = int(os.environ.get('VIRTUAL_MEMORY', 10**10))  # 10GB?
