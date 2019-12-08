#!/usr/bin/python3

import os
import yaml
import shutil
import multiprocessing
from concurrent.futures import ThreadPoolExecutor
from zipfile import ZipFile


_DEFAULT_EXE_NAME = 'mr_bayes'
_ENV_VAR = 'MR_BAYES_EXE'
# Short files :-)
_OUTPUT_EXTENSIONS = ('.ckp', '.con.tre', '.parts', '.tstat', '.vstat')


# MrBayes is not multiprocessing!
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
def _find_exe(default_exe, env_var):
    exe = os.getenv(env_var, default_exe)
    if not shutil.which(exe):
        print(_install_instructions.format(exe=default_exe, env_var=env_var))
        raise ValueError(f'No MrBayes installed! Tried {exe}')
    return exe


def _run_mr_bayes(mr_bayes_exe, run_dir, f):
    cmd = f"cd {run_dir}; {mr_bayes_exe} {f} > /dev/null"
    print(f"Cmd: {cmd}")
    os.system(cmd)


def run(locale=True, threads=None):
    mr_bayes_exe = _find_exe(_DEFAULT_EXE_NAME, _ENV_VAR)
    threads = threads or multiprocessing.cpu_count()

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
