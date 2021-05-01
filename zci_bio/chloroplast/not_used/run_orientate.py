#!/usr/bin/python3

import os
import yaml
import shutil
from concurrent.futures import ThreadPoolExecutor
import multiprocessing
from zipfile import ZipFile

_DEFAULT_EXE_NAME = 'mafft'
_ENV_VAR = 'MAFFT_EXE'

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


def _run_single(mafft_exe, filename, output_file, threads):
    cmd = f"{mafft_exe} --quiet --thread {threads} {filename} > {output_file}"
    print(f"Command: {cmd}")
    os.system(cmd)


def run(locale=True, threads=None):
    # Note: run from step's directory!!!
    mafft_exe = _find_exe(_DEFAULT_EXE_NAME, _ENV_VAR)
    threads = threads or multiprocessing.cpu_count()
    outputs = []

    # Files to align
    with open('finish.yml', 'r') as r:
        data = yaml.load(r, Loader=yaml.CLoader)  # dict with attrs: filename, short, max_seq_length
        align_files = data['align_files']

    with ThreadPoolExecutor(max_workers=threads) as executor:
        for seq_ident, part_name, direction in align_files:
            input_file = os.path.join(seq_ident, f"{part_name}_{direction}.fa")
            outputs.append(os.path.join(seq_ident, f"align_{part_name}_{direction}.fa"))
            executor.submit(_run_single, mafft_exe, input_file, outputs[-1], 1)

    # Zip files
    if not locale:
        with ZipFile('output.zip', 'w') as output:
            for f in outputs:
                output.write(f)


if __name__ == '__main__':
    import sys
    run(locale=False, threads=int(sys.argv[1]) if len(sys.argv) > 1 else None)
