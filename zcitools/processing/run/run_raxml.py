#!/usr/bin/python3

import os
import yaml
import shutil
import multiprocessing
from zipfile import ZipFile


_DEFAULT_EXE_NAME = 'raxmlHPC-PTHREADS'  # probably raxmlHPC-PTHREADS-AVX2
_ENV_VAR = 'RAXML_EXE'
_OUTPUT_FILES = (
    'RAxML_bestTree.raxml_output',
    'RAxML_bipartitionsBranchLabels.raxml_output',
    'RAxML_bipartitions.raxml_output',
    'RAxML_bootstrap.raxml_output',
    'RAxML_info.raxml_output')


# Calculation strategy:
#  - Just run calculations with given number of threads. Fist short, than long.
#
# Notes:
#  - Output files are create in directory where RAxML is run with extension raxml_output (set in command with switch -n)
#  - Files are listed in _OUTPUT_FILES

_install_instructions = """
RAxML is not installed!
Check web page https://cme.h-its.org/exelixis/web/software/raxml/index.html for installation instructions.

There are two ways for this script to locate executable to run:
 - environment variable {env_var} points to executable location,
 - or executable is called {exe} and placed on the PATH.
"""


# Note: it would be good that all scripts accept same format envs
def _find_exe(default_exe, env_var):
    exe = os.getenv(env_var, default_exe)
    if not shutil.which(exe):
        print(_install_instructions.format(exe=default_exe, env_var=env_var))
        raise ValueError(f'No RAxML installed! Tried {exe}')
    return exe


def run(locale=True, threads=None):
    raxml_exe = _find_exe(_DEFAULT_EXE_NAME, _ENV_VAR)
    threads = threads or multiprocessing.cpu_count()

    # Files to run
    with open('finish.yml', 'r') as r:
        data_files = yaml.load(r, Loader=yaml.CLoader)  # dict with attrs: filename, short
    files = [d['filename'] for d in data_files if d['short']] + \
        [d['filename'] for d in data_files if not d['short']]

    output_dirs = []
    step_dir = os.getcwd()
    for f in files:
        _dir, f = os.path.split(f)
        output_dirs.append(_dir)
        os.chdir(_dir)
        # ToDo: other arguments as parameters
        cmd = f"{raxml_exe} -s {f} -n raxml_output -m GTRGAMMA -f a -# 1000 -x 12345 -p 12345 -T {threads} > /dev/null"
        print(f"Cmd: {cmd}")
        os.system(cmd)
        os.chdir(step_dir)

    # Zip files
    if not locale:
        with ZipFile('output.zip', 'w') as output:
            for f in files:
                d = os.path.dirname(f)
                for x in _OUTPUT_FILES:
                    output.write(os.path.join(d, x))


if __name__ == '__main__':
    import sys
    run(locale=False, threads=int(sys.argv[1]) if len(sys.argv) > 1 else None)
