#!/usr/bin/python3

import os
import re
import subprocess
pattern = re.compile(r'run_.*\.py')

for d in os.listdir('.'):
    if os.path.isdir(d) and not os.path.exists(os.path.join(d, 'output.zip')):
        py_files = [f for f in os.listdir(d) if pattern.search(f)]
        if len(py_files) == 1:
            subprocess.run(['python3', py_files[0]], cwd=d)
