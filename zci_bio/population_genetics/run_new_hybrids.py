#!/usr/bin/python3

import os
import shutil
import yaml
import subprocess
from concurrent.futures import ThreadPoolExecutor
import multiprocessing
from zipfile import ZipFile, ZIP_DEFLATED


_DEFAULT_EXE_NAME = 'newhybrids-no-gui-linux.exe'
_ENV_VAR = 'NEW_HYBRIDS_EXE'
results_dir = 'results'

_install_instructions = """
NewHybrids is not installed or not properly set!
Check GitHub repository https://github.com/eriqande/newhybrids for installation instructions.

There are two ways for this script to locate executable to run:
 - environment variable {env_var} points to executable location,
 - or executable is called {exe} and placed on the PATH.
"""


# Note: it would be good that all scripts accept same format envs
def _find_exe(default_exe, env_var):
    exe = os.getenv(env_var, default_exe)
    if not shutil.which(exe):
        print(_install_instructions.format(exe=default_exe, env_var=env_var))
        raise ValueError(f'No NewHybrids installed! Tried {exe}')
    return exe


def _run_single(nh_exe, run_dir, ps, seed):
    # ps has attrs: data_file, gtyp_cat_file, theta_prior, pi_prior, burn_in, num_sweeps
    print(f"""Command: cd {run_dir}; {nh_exe} -d ../{ps['data_file']} -c {ps['gtyp_cat_file']}
--theta-prior {ps['theta_prior']} --pi-prior {ps['pi_prior']}
--burn-in {ps['burn_in']} --num-sweeps {ps['num_sweeps']} -s {seed[0]} {seed[1]} --no-gui > /dev/null""")
    subprocess.run([
        nh_exe,
        '-d', os.path.join('..', ps['data_file']),
        '-c', os.path.join('..', ps['gtyp_cat_file']),
        '--theta-prior', ps['theta_prior'],
        '--pi-prior', ps['pi_prior'],
        '--burn-in', str(ps['burn_in']),
        '--num-sweeps', str(ps['num_sweeps']),
        '-s', str(seed[0]), str(seed[1]),
        '--no-gui'], cwd=run_dir, stdout=subprocess.DEVNULL)


def _fix_prior(data, attr):
    if data[attr] not in ('Jeffreys', 'Uniform'):
        data[attr] = 'Uniform' if data[attr][0].lower() == 'u' else 'Jeffreys'


def run(locale=True, threads=None):
    # Note: run from step's directory!!!
    nh_exe = _find_exe(_DEFAULT_EXE_NAME, _ENV_VAR)
    threads = threads or multiprocessing.cpu_count()
    # outputs = []

    with open('finish.yml', 'r') as r:
        # Attrs: data_file, gtyp_cat_file, theta_prior, pi_prior, burn_in, num_sweeps
        params = yaml.load(r, Loader=yaml.CLoader)

    _fix_prior(params, 'theta_prior')
    _fix_prior(params, 'pi_prior')
    seed_dirs = [f for f in os.listdir('.') if f.startswith('seed_') and os.path.isdir(f)]

    with ThreadPoolExecutor(max_workers=threads) as executor:
        for run_dir in seed_dirs:
            seed = tuple(map(int, run_dir.split('_')[1:]))
            executor.submit(_run_single, nh_exe, os.path.abspath(run_dir), params, seed)

    #
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    elif os.path.isfile(results_dir):
        raise ValueError(f"Can't create directory, file exists with same name ({results_dir})!")

    # Copy result files
    to_copy = (
        'aa-EchoedGtypFreqCats.txt', 'aa-LociAndAlleles.txt', 'aa-Pi.hist', 'aa-PofZ.txt',
        'aa-ThetaAverages.txt', 'aa-Theta.hist', 'EchoedGtypData.txt')

    output_files = []
    for run_dir in seed_dirs:
        for f in to_copy:
            in_f = os.path.join(run_dir, f)
            seed = run_dir.split('_')[1:]
            f_parts = f.split('.')
            out_f = os.path.join(results_dir, f"{f_parts[0]}_{seed[0]}_{seed[1]}.{f_parts[1]}")
            output_files.append(out_f)
            shutil.copyfile(in_f, out_f)

    # Zip files
    if not locale:
        with ZipFile('output.zip', 'w', compression=ZIP_DEFLATED) as output:
            for f in output_files:
                output.write(f)


if __name__ == '__main__':
    import sys
    run(locale=False, threads=int(sys.argv[1]) if len(sys.argv) > 1 else None)
