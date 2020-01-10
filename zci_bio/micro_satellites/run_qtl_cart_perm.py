#!/usr/bin/python3

import os
import shutil
import yaml
from concurrent.futures import ThreadPoolExecutor
import multiprocessing
from zipfile import ZipFile

_DEFAULT_EXE_NAME = 'Zmapqtl'
_ENV_VAR = 'QTL_CART_PERMUTATION_EXE'

_install_instructions = """
QTL Carthographer is not installed or not properly set!
Check web page https://brcwebportal.cos.ncsu.edu/qtlcart/ for installation instructions.

There are two ways for this script to locate executable to run:
 - environment variable {env_var} points to executable location,
 - or executable is called {exe} and placed on the PATH.
"""


# Note: it would be good that all scripts accept same format envs
def _find_exe(default_exe, env_var):
    exe = os.getenv(env_var, default_exe)
    if not shutil.which(exe):
        print(_install_instructions.format(exe=default_exe, env_var=env_var))
        raise ValueError(f'No QTL Cartographer installed! Tried {exe}')
    return exe


def _run_single(qtl_exe, run_dir, num_perm):
    cmd = f"cd {run_dir}; {qtl_exe} -r {num_perm} -A > /dev/null"
    print(f"Command: {cmd}")
    os.system(cmd)


def _find_5_perc(filename):
    store = False
    rows = []  # floats
    with open(filename, 'r') as r:
        for line in r.readlines():
            fields = line.split()
            if store and len(fields) == 2:
                rows.append(float(fields[1]))
            if len(fields) >= 3 and fields[0] == '#' and fields[1] == 'Repetition' and fields[2] == 'GlobalMax':
                store = True
    rows = sorted(rows)
    return rows[len(rows) // 20]


def run(locale=True, threads=None):
    # Note: run from step's directory!!!
    qtl_exe = _find_exe(_DEFAULT_EXE_NAME, _ENV_VAR)
    threads = threads or multiprocessing.cpu_count()
    # outputs = []

    with open('finish.yml', 'r') as r:
        data = yaml.load(r, Loader=yaml.CLoader)

    num_perm = data['permutations']
    trait_dirs = data['trait_dirs']

    with ThreadPoolExecutor(max_workers=threads) as executor:
        for run_dir in [os.path.abspath(d) for d in trait_dirs]:
            executor.submit(_run_single, qtl_exe, run_dir, num_perm)

    #
    results_dir = 'results'
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    elif os.path.isfile(results_dir):
        raise ValueError(f"Can't create directory, file exists with same name ({results_dir})!")

    # Copy qtlcart.z6e files
    summary_data = []
    output_files = [os.path.join(results_dir, "TLR.txt"), os.path.join(results_dir, "TLOD.txt")]
    for run_dir in sorted(trait_dirs):
        input_f = os.path.join(run_dir, 'qtlcart.z6e')
        output_files.append(os.path.join(results_dir, f"{run_dir}.txt"))
        shutil.copyfile(input_f, output_files[-1])
        summary_data.append((run_dir, _find_5_perc(input_f)))

    #
    with open(output_files[0], 'w') as out:
        for r, f in summary_data:
            out.write(f"{r} {f}\n")
    with open(output_files[1], 'w') as out:
        d = 4.60517018599  # 2 * ln(10)
        for r, f in summary_data:
            out.write(f"{r} {f / d}\n")

    # Zip files
    # if not locale:
    with ZipFile('output.zip', 'w') as output:
        for f in output_files:
            output.write(f, arcname=os.path.basename(f))


if __name__ == '__main__':
    import sys
    run(locale=False, threads=int(sys.argv[1]) if len(sys.argv) > 1 else None)
