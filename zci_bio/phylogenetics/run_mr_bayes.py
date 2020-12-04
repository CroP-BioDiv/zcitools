#!/usr/bin/python3

import os
import yaml
import shutil
import multiprocessing
import subprocess
from concurrent.futures import ThreadPoolExecutor
from zipfile import ZipFile


_DEFAULT_EXE_NAME = 'mr_bayes'
_ENV_VAR = 'MR_BAYES_EXE'
_DEFAULT_EXE_NAME_MPI = 'mr_bayes_mpi'
_ENV_VAR_MPI = 'MR_BAYES_MPI_EXE'
# Short files :-)
_OUTPUT_EXTENSIONS = ('.ckp', '.con.tre', '.parts', '.tstat', '.vstat')


# MrBayes (base version) is not multiprocessing!
# To make is multiprocessing compiling with MPI is needed
#   https://github.com/NBISweden/MrBayes/blob/develop/INSTALL
#
# Calculation strategy:
#  - Just run calculations in a pool. Fist long than short.

_install_instructions = """
MrBayes is not installed!
Check web page https://nbisweden.github.io/MrBayes/index.html for installation instructions.

There are two ways for this script to locate executable to run:
 - environment variable {env_var} points to executable location,
 - or executable is called {exe} and placed on the PATH.
"""


# Note: it would be good that all scripts accept same format envs
def _find_exe(default_exe, env_var, to_raise=True):
    exe = os.getenv(env_var, default_exe)
    if not shutil.which(exe):
        print(_install_instructions.format(exe=default_exe, env_var=env_var))
        if to_raise:
            raise ValueError(f'No MrBayes installed! Tried {exe}')
        return
    return exe


def _run_mr_bayes(exe, run_dir, f):
    cmd = [exe, f]
    print(f"Cmd: cd {run_dir}; {' '.join(cmd)}")
    subprocess.run(cmd, cwd=run_dir)  # , stdout=subprocess.DEVNULL)


def _run_mr_bayes_mpi(exe, run_dir, f, threads):
    cmd = ['mpirun', '--use-hwthread-cpus', '-np', str(threads), exe, f]
    print(f"Cmd: cd {run_dir}; {' '.join(cmd)}")
    subprocess.run(cmd, cwd=run_dir)  # , stdout=subprocess.DEVNULL)


def run(locale=True, threads=None):
    threads = threads or multiprocessing.cpu_count()
    mr_bayes_mpi_exe = _find_exe(_DEFAULT_EXE_NAME_MPI, _ENV_VAR_MPI, to_raise=False) if threads > 1 else None
    mr_bayes_exe = _find_exe(_DEFAULT_EXE_NAME, _ENV_VAR)

    # Files to run
    with open('finish.yml', 'r') as r:
        data_files = yaml.load(r, Loader=yaml.CLoader)  # dict with attrs: filename, short
    files = [d['filename'] for d in data_files if not d['short']] + \
        [d['filename'] for d in data_files if d['short']]

    # Find absoulte paths so that threads can set it's own working dir safe
    dir_data = []  # Tuples: filename, relative dir, absolute dir
    for f in files:
        _dir, f = os.path.split(f)
        dir_data.append((f, _dir, os.path.abspath(_dir)))

    step_dir = os.path.abspath(os.getcwd())
    if mr_bayes_mpi_exe:
        for f, _, abs_dir in dir_data:
            _run_mr_bayes_mpi(mr_bayes_mpi_exe, abs_dir, f, threads)
    else:
        with ThreadPoolExecutor(max_workers=threads) as executor:
            for f, _, abs_dir in dir_data:
                executor.submit(_run_mr_bayes, mr_bayes_exe, abs_dir, f)

    # Zip files
    if not locale:
        os.chdir(step_dir)
        with ZipFile('output.zip', 'w') as output:
            for f, loc_dir, _ in dir_data:
                f = f.replace('.nex', '')
                for ext in _OUTPUT_EXTENSIONS:
                    output.write(os.path.join(loc_dir, f + ext))


if __name__ == '__main__':
    import sys
    run(locale=False, threads=int(sys.argv[1]) if len(sys.argv) > 1 else None)
