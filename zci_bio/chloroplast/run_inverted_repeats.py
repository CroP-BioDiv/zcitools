#!/usr/bin/python3

import os
import shutil
from concurrent.futures import ThreadPoolExecutor
import multiprocessing
from zipfile import ZipFile
from collections import namedtuple

_DEFAULT_EXE_NAME = 'repeat-match'
_ENV_VAR = 'MUMMER_REPEAT_MATCH_EXE'

# Finds chloroplast inverted repeats using mummer's repeat-match program.
# Tolerates few different bases.
# Result is txt file with 4 numbers of format:
# <ira_start> <ira_end>
# <irb_start> <irb_end>
#
# Calculation strategy: jobs are run in parallel.

_install_instructions = """
Mummer is not installed!
Check web page https://mummer4.github.io/install/install.html for installation instructions.

There are two ways for this script to locate executable to run:
 - environment variable {env_var} points to executable location,
 - or executable is called {exe} and placed on the PATH.
"""

_Repeat = namedtuple('_Repeat', 'first_start, second_start, length, inverted')


# Note: it would be good that all scripts accept same format envs
def _find_exe(default_exe, env_var):
    exe = os.getenv(env_var, default_exe)
    if not shutil.which(exe):
        print(_install_instructions.format(exe=default_exe, env_var=env_var))
        raise ValueError(f'No Mummer installed! Tried {exe}')
    return exe


def _run_single(mummer_exe, input_filename, output_filename, result_filename, n):
    cmd = f"{mummer_exe} -n {n} {input_filename} > {output_filename}"
    print(f"Command: {cmd}")
    os.system(cmd)

    # Read mummer output
    repeats = []  # List of _Repeat objects
    with open(output_filename, 'r') as output:
        read = False
        for line in output:
            fields = line.strip().split()
            if read:
                if fields[1][-1] == 'r':
                    s2 = int(fields[1][:-1])
                    inverted = True
                else:
                    s2 = int(fields[1])
                    inverted = False
                repeats.append(_Repeat(int(fields[0]), s2, int(fields[2]), inverted=inverted))
            else:
                if fields[0] == 'Start1':
                    read = True

    # Write result
    with open(result_filename, 'w') as res:  # Note: file is always created
        if repeats:
            irs = None
            if len(repeats) == 1:
                irs = repeats[0]
            else:
                # Concatenate repeats into one
                # ToDo: check lot of things!!! Gap, change, length of change, ...
                repeats.sort(key=lambda r: r.first_start)
                first = repeats[0]
                last = repeats[-1]
                irs = _Repeat(
                    first.first_start, first.second_start,
                    last.first_start - first.first_start + last.length,
                    inverted=first.inverted)
            if irs.inverted:
                res.write(f"{irs.first_start} {irs.first_start + irs.length - 1}\n")
                res.write(f"{irs.second_start - irs.length + 1} {irs.second_start}")
            else:
                print(f"Error: no inverted repeats found for input file {input_filename}!")
                print(repeats)
        #
        else:
            print(f"Error: no repeats at all found for input file {input_filename}!")


def run(locale=True, threads=None, min_length=15000):
    # Note: run from step's directory!!!
    fa_files = [f for f in os.listdir('.') if f.endswith('.fa')]
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
            executor.submit(_run_single, mummer_exe, f, outputs[-2], outputs[-1], min_length // 20)
            # _run_single(mummer_exe, f, outputs[-2], outputs[-1], min_length)

    # Zip files
    if not locale:
        with ZipFile('output.zip', 'w') as output:
            for f in outputs:
                output.write(f)


if __name__ == '__main__':
    import sys
    run(locale=False, threads=int(sys.argv[1]) if len(sys.argv) > 1 else None)
