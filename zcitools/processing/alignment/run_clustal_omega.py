#!/usr/bin/python3

import os
import yaml
from concurrent.futures import ThreadPoolExecutor
import multiprocessing
from zipfile import ZipFile

# From my observations Clustal Omega usage of more threads is limited to part of calculation.
# Because of that is seems to me that lot of small jobs is good to run parallel and use more
# threads only for long jobs.
#
# Calculation strategy:
#  - First run short sequences on one thread. Sort them from longer (short) to shorter.
#  - Than run long sequences on more threads. Sort them from shorter to longer.


# Note: it would be good that all scripts accept same format envs
def _find_exe(env_var_prefix, default_exe):
    return os.getenv(f'{env_var_prefix}_EXE', default_exe)


def _alignment_file(f):
    return os.path.join(os.path.dirname(f), 'alignment.phy')


def _run_single(clustalo_exe, filename, output_file):
    cmd = f"{clustalo_exe} -i {filename} -o {output_file} --outfmt=phy --threads=1"
    print(f"Command: {cmd}")
    os.system(cmd)


def run(locale=True, threads=None):
    # Note: run from step's directory!!!
    clustalo_exe = _find_exe('CLUSTAL_OMEGA', 'clustalo')
    if not clustalo_exe:
        print('No Clustal Omega exe!!!')
        return

    threads = threads or multiprocessing.cpu_count()
    outputs = []

    # Files to run
    with open('finish.yml', 'r') as r:
        seq_files = yaml.load(r, Loader=yaml.CLoader)  # dict with attrs short, long
    short_files = sorted((d for d in seq_files if d['short']), key=lambda x: -x['max_seq_length'])
    long_files = sorted((d for d in seq_files if not d['short']), key=lambda x: x['max_seq_length'])

    if short_files:
        with ThreadPoolExecutor(max_workers=threads) as executor:
            for d in short_files:
                outputs.append(_alignment_file(d['filename']))
                executor.submit(_run_single, clustalo_exe, d['filename'], outputs[-1])

    for d in long_files:
        for d in long_files:
            outputs.append(_alignment_file(d['filename']))
            cmd = f"{clustalo_exe} -i {d['filename']} -o {outputs[-1]} --outfmt=phy --threads={threads}"
            print(f"Command: {cmd}")
            os.system(cmd)

    # Zip files
    if not locale:
        with ZipFile('output.zip', 'w') as output:
            for f in outputs:
                output.write(f)


if __name__ == '__main__':
    import sys
    run(locale=False, threads=int(sys.argv[1]) if len(sys.argv) > 1 else None)
