#!/usr/bin/python3

import os
import sys
import re
import subprocess
import argparse
pattern = re.compile(r'run_.*\.py')


def _run(_dir, py_module, num_threads):
    subprocess.run(['python3', py_module, num_threads], cwd=_dir)


#
parser = argparse.ArgumentParser(description="Finds bin paths")
parser.add_argument('directories', nargs='*', help='Directories to run. Default all subdirectories.')
parser.add_argument('-t', '--threads', type=int, help="Number of threads per run")
parser.add_argument('-r', '--runs', default=1, type=int, help="Numbner of runs to run in parallel")
params = parser.parse_args()

to_run = []
for d in (params.directories or os.listdir('.')):
    if os.path.isdir(d) and not os.path.exists(os.path.join(d, 'output.zip')):
        py_files = [f for f in os.listdir(d) if pattern.search(f)]
        if len(py_files) == 1:
            to_run.append((d, py_files[0]))

if params.runs <= 1 or not params.threads:
    for d, py_module in to_run:
        subprocess.run(['python3', py_module], cwd=d)
else:
    from concurrent.futures import ThreadPoolExecutor

    num_threads = str(params.threads)
    with ThreadPoolExecutor(max_workers=params.runs) as executor:
        for d, py_module in to_run:
            executor.submit(_run, d, py_module, num_threads)
