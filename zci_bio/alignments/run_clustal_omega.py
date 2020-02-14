#!/usr/bin/python3

import os
import yaml
import shutil
from concurrent.futures import ThreadPoolExecutor
import multiprocessing
from zipfile import ZipFile

_DEFAULT_EXE_NAME = 'clustalo'
_ENV_VAR = 'CLUSTAL_OMEGA_EXE'

# From my observations Clustal Omega usage of more threads is limited to part of calculation.
# Because of that is seems to me that lot of small jobs is good to run parallel and use more
# threads only for long jobs.
# Also long jobs are run with more threads but in parallel
#
# Calculation strategy:
#  - First run short sequences on one thread. Sort them from longer (short) to shorter.
#  - Than run long sequences on more threads. Sort them from shorter to longer.

_install_instructions = """
Clustal Omega is not installed!
Check web page http://www.clustal.org/omega/ for installation instructions.

There are two ways for this script to locate executable to run:
 - environment variable {env_var} points to executable location,
 - or executable is called {exe} and placed on the PATH.
"""


# Note: it would be good that all scripts accept same format envs
def _find_exe(default_exe, env_var):
    exe = os.getenv(env_var, default_exe)
    if not shutil.which(exe):
        print(_install_instructions.format(exe=default_exe, env_var=env_var))
        raise ValueError(f'No Clustal Omega installed! Tried {exe}')
    return exe


def _alignment_file(f):
    return os.path.join(os.path.dirname(f), 'alignment.phy')


def _run_single(clustalo_exe, filename, output_file, threads):
    cmd = f"{clustalo_exe} -i {filename} -o {output_file} --outfmt=phy --threads={threads}"
    print(f"Command: {cmd}")
    os.system(cmd)


def run(locale=True, threads=None, long_parallel=3):  # usually there are max 3 long files: w, gc, cc
    # Note: run from step's directory!!!
    clustalo_exe = _find_exe(_DEFAULT_EXE_NAME, _ENV_VAR)
    threads = threads or multiprocessing.cpu_count()
    outputs = []

    # Files to run
    with open('finish.yml', 'r') as r:
        seq_files = yaml.load(r, Loader=yaml.CLoader)  # dict with attrs: filename, short, max_seq_length
    short_files = sorted((d for d in seq_files if d['short']), key=lambda x: -x['max_seq_length'])
    long_files = sorted((d for d in seq_files if not d['short']), key=lambda x: x['max_seq_length'])

    if short_files:
        with ThreadPoolExecutor(max_workers=threads) as executor:
            for d in short_files:
                outputs.append(_alignment_file(d['filename']))
                executor.submit(_run_single, clustalo_exe, d['filename'], outputs[-1], 1)

    if long_files:
        long_parallel = min(long_parallel, len(long_files))
        threads = threads // long_parallel
        with ThreadPoolExecutor(max_workers=long_parallel) as executor:
            for d in long_files:
                outputs.append(_alignment_file(d['filename']))
                executor.submit(_run_single, clustalo_exe, d['filename'], outputs[-1], threads)

    # Zip files
    if not locale:
        with ZipFile('output.zip', 'w') as output:
            for f in outputs:
                output.write(f)


if __name__ == '__main__':
    import sys
    run(locale=False,
        threads=int(sys.argv[1]) if len(sys.argv) > 1 else None,
        long_parallel=int(sys.argv[2]) if len(sys.argv) > 2 else 3)
