#!/usr/bin/python3

import os
import yaml
import shutil
from concurrent.futures import ThreadPoolExecutor
import multiprocessing
from zipfile import ZipFile

_DEFAULT_EXE_NAME = 'mafft'
_ENV_VAR = 'MAFFT_EXE'

# From my observations ...
#
# Calculation strategy:
#  - First run short sequences on one thread. Sort them from longer (short) to shorter.
#  - Than run long sequences sequential

_install_instructions = """
MAFFT is not installed!
Check web page https://mafft.cbrc.jp/alignment/software/ for installation instructions.

There are two ways for this script to locate executable to run:
 - environment variable {env_var} points to executable location,
 - or executable is called {exe} and placed on the PATH.
"""


# Note: it would be good that all scripts accept same format envs
def _find_exe(default_exe, env_var):
    exe = os.getenv(env_var, default_exe)
    if not shutil.which(exe):
        print(_install_instructions.format(exe=default_exe, env_var=env_var))
        raise ValueError(f'No MAFFT installed! Tried {exe}')
    return exe


def _alignment_file(f):
    return os.path.join(os.path.dirname(f), 'alignment.phy')


def _run_single(mafft_exe, filename, namelength, output_file, threads):
    cmd = f"{mafft_exe} --maxiterate 10 --phylipout --namelength {namelength} " + \
        f"--thread {threads} {filename} > {output_file}"
    print(f"Command: {cmd}")
    os.system(cmd)


def run(locale=True, threads=None):
    # Note: run from step's directory!!!
    mafft_exe = _find_exe(_DEFAULT_EXE_NAME, _ENV_VAR)
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
                executor.submit(_run_single, mafft_exe, d['filename'], d['namelength'], outputs[-1], 1)

    for d in long_files:
        outputs.append(_alignment_file(d['filename']))
        _run_single(mafft_exe, d['filename'], d['namelength'], outputs[-1], threads)

    # Zip files
    if not locale:
        with ZipFile('output.zip', 'w') as output:
            for f in outputs:
                output.write(f)


if __name__ == '__main__':
    import sys
    run(locale=False, threads=int(sys.argv[1]) if len(sys.argv) > 1 else None)
