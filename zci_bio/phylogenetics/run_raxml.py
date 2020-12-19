#!/usr/bin/python3

import os
from concurrent.futures import ThreadPoolExecutor
try:                 # Run locally, with whole project
    import common_utils.exec_utils as exec_utils
except ImportError:  # Run standalone, on server
    import exec_utils


_DEFAULT_EXE_NAME = 'raxml_threads'  # probably raxmlHPC-PTHREADS-AVX2
_ENV_VAR = 'RAXML_EXE'
_OUTPUT_FILES = (
    'RAxML_bestTree.raxml_output',
    'RAxML_bipartitionsBranchLabels.raxml_output',
    'RAxML_bipartitions.raxml_output',
    'RAxML_bootstrap.raxml_output',
    'RAxML_info.raxml_output')
_stat_args = ['-n', 'raxml_output', '-m', 'GTRGAMMA', '-f', 'a', '-N', '1000']
_ps_args = ['-q', 'partitions.ind']


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


def _run(exe, run_dir, input_file, seed, threads):
    _ps = os.path.isfile(os.path.join(run_dir, 'partitions.ind'))
    cmd = [exe, '-s', input_file, '-T', str(threads), '-x', str(seed), '-p', str(seed)] + \
        (_ps_args if _ps else []) + \
        _stat_args
    exec_utils.run_cmd(cmd, cwd=run_dir)


def run(locale=True, threads=None):
    raxml_exe = exec_utils.find_exe(_DEFAULT_EXE_NAME, _ENV_VAR, _install_instructions, 'RAxML')
    threads = threads or exec_utils.get_num_logical_threads()
    log_run = exec_utils.LogRun(threads=threads, use_mpi=use_mpi, raxml_exe=raxml_exe)

    # Files to run
    data_files = exec_utils.load_finish_yml()  # dict with attrs: filename, short
    short_files = sorted((d for d in data_files if d['short']), key=lambda x: -x['length'])
    long_files = sorted((d for d in data_files if not d['short']), key=lambda x: x['length'])

    if short_files:
        with ThreadPoolExecutor(max_workers=threads) as executor:
            for d in short_files:
                _dir, f = os.path.split(d['filename'])
                executor.submit(_run, raxml_exe, os.path.abspath(_dir), f, d['seed'], 1)

    if long_files:
        for d in long_files:
            _dir, f = os.path.split(d['filename'])
            _run(raxml_exe, os.path.abspath(_dir), f, d['seed'], threads)

    log_run.finish()  # Creates run_info.txt file

    # Zip files
    if not locale:
        dirs = [os.path.dirname(f['filename']) for f in data_files]
        exec_utils.zip_output([os.path.join(d, x) for d, x in itertools.product(dirs, _OUTPUT_FILES)])


if __name__ == '__main__':
    import sys
    run(locale=False, threads=int(sys.argv[1]) if len(sys.argv) > 1 else None)
