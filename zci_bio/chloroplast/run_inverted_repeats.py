#!/usr/bin/python3

import os
import shutil
from concurrent.futures import ThreadPoolExecutor
import multiprocessing
from zipfile import ZipFile

_DEFAULT_EXE_NAME = 'repeat-match'
_ENV_VAR = 'MUMMER_REPEAT_MATCH_EXE'
_CALC_DIR = lambda f: os.path.join('run_dir', f)

# Finds chloroplast inverted repeats using mummer's repeat-match program.
# Tolerates few different bases.
# Result is txt file with 4 numbers of format:
# <ira_start> <ira_end>
# <irb_start> <irb_end>
#
# Calculation strategy: jobs are run in parallel.

_install_instructions = """
MUMmer is not installed!
Check web page https://mummer4.github.io/install/install.html for installation instructions.

There are two ways for this script to locate executable to run:
 - environment variable {env_var} points to executable location,
 - or executable is called {exe} and placed on the PATH.
"""


# Note: it would be good that all scripts accept same format envs
def _find_exe(default_exe, env_var):
    exe = os.getenv(env_var, default_exe)
    if not shutil.which(exe):
        print(_install_instructions.format(exe=default_exe, env_var=env_var))
        raise ValueError(f'No MUMmer installed! Tried {exe}')
    return exe


def _run_single(mummer_exe, input_filename, output_filename, result_filename, n):
    cmd = f"{mummer_exe} -n {n} {_CALC_DIR(input_filename)} > {_CALC_DIR(output_filename)}"
    print(f"Command: {cmd}")
    os.system(cmd)


def run(locale=True, threads=None, min_length=15000):
    # Note: run from step's directory!!!
    fa_files = [f for f in os.listdir('run_dir') if f.endswith('.fa')]
    if not fa_files:
        print('No input files!')

    mummer_exe = _find_exe(_DEFAULT_EXE_NAME, _ENV_VAR)
    threads = threads or multiprocessing.cpu_count()
    outputs = []

    with ThreadPoolExecutor(max_workers=threads) as executor:
        for f in fa_files:
            f_base = f[:-3]
            outputs.append(f'{f_base}.out')
            outputs.append(f'{f_base}.irs')
            executor.submit(_run_single, mummer_exe, f, outputs[-2], outputs[-1], min_length // 100)

    # Zip files
    if not locale:
        with ZipFile('output.zip', 'w') as output:
            for f in outputs:
                output.write(_CALC_DIR(f))


if __name__ == '__main__':
    import sys
    run(locale=False, threads=int(sys.argv[1]) if len(sys.argv) > 1 else None)
