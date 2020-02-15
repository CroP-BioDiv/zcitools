#!/usr/bin/python3

import os
import yaml
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
import multiprocessing
from zipfile import ZipFile


_DEFAULT_EXE_NAME = 'raxml_threads'  # probably raxmlHPC-PTHREADS-AVX2
_ENV_VAR = 'RAXML_EXE'
_OUTPUT_FILES = (
    'RAxML_bestTree.raxml_output',
    'RAxML_bipartitionsBranchLabels.raxml_output',
    'RAxML_bipartitions.raxml_output',
    'RAxML_bootstrap.raxml_output',
    'RAxML_info.raxml_output')


# Calculation strategy:
#  - Short alignment run in parallel
#  - Long ones serial
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


_stat_cmd = "-n raxml_output -m GTRGAMMA -f a -# 1000 -x 12345 -p 12345 > /dev/null"
_stat_args = [
    '-n', 'raxml_output',
    '-m', 'GTRGAMMA',
    '-f', 'a',
    '-#', '1000',
    '-x', '12345',
    '-p', '12345']


def _run(raxml_exe, run_dir, input_file, threads):
    if os.path.isfile(os.path.join(run_dir, 'partitions.ind')):
        print(f"Command: cd {run_dir}; {raxml_exe} -s {input_file} -T {threads} -q partitions.ind {_stat_cmd}")
        subprocess.run([raxml_exe, '-s', input_file, '-T', str(threads), '-q', 'partitions.ind'] + _stat_args,
                       cwd=run_dir, stdout=subprocess.DEVNULL)
    else:
        print(f"Command: cd {run_dir}; {raxml_exe} -s {input_file} -T {threads} {_stat_cmd}")
        subprocess.run([raxml_exe, '-s', input_file, '-T', str(threads)] + _stat_args,
                       cwd=run_dir, stdout=subprocess.DEVNULL)


def run(locale=True, threads=None):
    raxml_exe = _find_exe(_DEFAULT_EXE_NAME, _ENV_VAR)
    threads = threads or multiprocessing.cpu_count()

    # Files to run
    with open('finish.yml', 'r') as r:
        data_files = yaml.load(r, Loader=yaml.CLoader)  # dict with attrs: filename, short
    short_files = sorted((d for d in data_files if d['short']), key=lambda x: -x['length'])
    long_files = sorted((d for d in data_files if not d['short']), key=lambda x: x['length'])

    if short_files:
        with ThreadPoolExecutor(max_workers=threads) as executor:
            for d in short_files:
                _dir, f = os.path.split(d['filename'])
                executor.submit(_run, raxml_exe, os.path.abspath(_dir), f, 1)

    if long_files:
        for d in long_files:
            _dir, f = os.path.split(d['filename'])
            _run(raxml_exe, os.path.abspath(_dir), f, threads)

    # Zip files
    if not locale:
        with ZipFile('output.zip', 'w') as output:
            for f in data_files:
                d = os.path.dirname(f['filename'])
                for x in _OUTPUT_FILES:
                    output.write(os.path.join(d, x))


if __name__ == '__main__':
    import sys
    run(locale=False, threads=int(sys.argv[1]) if len(sys.argv) > 1 else None)
