#!/usr/bin/python3

import os
import sys
import re
import subprocess
pattern = re.compile(r'run_.*\.py')

_dirs = sys.argv[1:] if len(sys.argv) > 1 else os.listdir('.')

for d in _dirs:
    if os.path.isdir(d) and not os.path.exists(os.path.join(d, 'output.zip')):
        py_files = [f for f in os.listdir(d) if pattern.search(f)]
        if len(py_files) == 1:
            subprocess.run(['python3', py_files[0]], cwd=d)
