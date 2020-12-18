#!/usr/bin/python3

import os
import shutil
import itertools
import time
from concurrent.futures import ThreadPoolExecutor
try:                 # Run locally, with whole project
    import common_utils.exec_utils as exec_utils
except ImportError:  # Run standalone, on server
    import exec_utils

_DEFAULT_EXE_NAME = 'mr_bayes'
_ENV_VAR = 'MR_BAYES_EXE'
_DEFAULT_EXE_NAME_MPI = 'mr_bayes_mpi'
_ENV_VAR_MPI = 'MR_BAYES_MPI_EXE'
# Short files :-)
# '.mcmc', '.trprobs'
_OUTPUT_EXTENSIONS = ('.ckp', '.con.tre', '.parts', '.run1.p', '.run1.t', '.run2.p', '.run2.t', '.tstat', '.vstat')

# Calculation strategy for MrBayes multithread executable:
#  - Sort jobs descending by nchains.
#  - Find how large force pool can be.
#  - Run jobs in upper order
#
# Calculation strategy for MrBayes single processor executable:
#  - First execute long jobs, than short

_install_instructions = """
MrBayes is not installed!
Check web page https://nbisweden.github.io/MrBayes/index.html for installation instructions.

There are two ways for this script to locate executable to run:
 - environment variable {env_var} points to executable location,
 - or executable is called {exe} and placed on the PATH.
"""


def _run_mr_bayes(exe, run_dir, f):
    exec_utils.run_cmd([exe, f], cwd=run_dir)


def _run_mr_bayes_mpi(exe, run_dir, f, nchains, threads, job_idx):
    # There must be at least as many chains as MPI processors
    # Chain consists of two chains :-)
    threads = min(threads, nchains * 2)
    cmd = ['mpirun', '-np', threads, exe, f]
    # Use logical cores if needed
    if threads > exec_utils.get_num_threads():
        cmd.insert(1, '--use-hwthread-cpus')
    # if job_idx:
    #     cmd.insert(1, '--tag-output')
    exec_utils.run_cmd(cmd, cwd=run_dir)


def run(locale=True, threads=None, use_mpi=True):
    threads = threads or exec_utils.get_num_logical_threads()
    # find_exe(default_exe, env_var, install_instructions, raise_desc)
    mr_bayes_mpi_exe = exec_utils.find_exe(_DEFAULT_EXE_NAME_MPI, _ENV_VAR_MPI, _install_instructions, None) \
        if (use_mpi and threads > 1) else None
    mr_bayes_exe = exec_utils.find_exe(_DEFAULT_EXE_NAME, _ENV_VAR, _install_instructions, 'MrBayes')
    step_dir = os.path.abspath(os.getcwd())  # Store current directory, for zipping after processing is done

    log_run = exec_utils.LogRun(
        threads=threads, use_mpi=use_mpi, mr_bayes_exe=mr_bayes_exe, mr_bayes_mpi_exe=mr_bayes_mpi_exe)

    # Files to run
    data_files = exec_utils.load_finish_yml()  # dict with attrs: filename, short

    if mr_bayes_mpi_exe:
        if len(data_files) == 1:
            d = data_files[0]
            _dir, f = os.path.split(d['filename'])
            _run_mr_bayes_mpi(mr_bayes_mpi_exe, os.path.abspath(_dir), f, d['nchains'], threads, None)
        else:
            data_files = sorted(data_files, key=lambda d: (-d['nchains'], d['short']))

            # Find how large thread pool can be
            cum_sum = list(itertools.accumulate(2 * d['nchains'] for d in data_files))
            for max_workers, sum_cs in enumerate(cum_sum):
                if sum_cs > threads:
                    break
            if max_workers == 0:
                max_workers = 1
            elif max_workers == len(data_files) - 1 and cum_sum[-1] <= threads:
                max_workers = len(data_files)

            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                for job_idx, d in enumerate(data_files):
                    _dir, f = os.path.split(d['filename'])
                    executor.submit(_run_mr_bayes_mpi, mr_bayes_mpi_exe,
                                    os.path.abspath(_dir), f, d['nchains'], threads, job_idx + 1)

    else:
        data_files = [d for d in data_files if not d['short']] + [d for d in data_files if d['short']]
        with ThreadPoolExecutor(max_workers=threads) as executor:
            for d in data_files:
                # Find absolute paths so that threads can set it's own working dir safe
                _dir, f = os.path.split(d['filename'])
                executor.submit(_run_mr_bayes, mr_bayes_exe, os.path.abspath(_dir), f)

    # Zip files
    if not locale:
        base_names = [d['result_prefix'] for d in data_files]
        exec_utils.zip_files([f + ext for f, ext in itertools.product(base_names, _OUTPUT_EXTENSIONS)],
                             cwd=step_dir, skip_missing=True)

    #
    log_run.finish()


if __name__ == '__main__':
    import sys
    threads = int(sys.argv[1]) if len(sys.argv) > 1 else None
    use_mpi = (len(sys.argv) <= 2)
    run(locale=False, threads=threads, use_mpi=use_mpi)
