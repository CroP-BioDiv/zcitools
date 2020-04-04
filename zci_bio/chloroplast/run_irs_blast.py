#!/usr/bin/python3

import os
import yaml
import shutil
from concurrent.futures import ThreadPoolExecutor
import multiprocessing
from zipfile import ZipFile

_DEFAULT_BLASTN_EXE_NAME = 'blastn'
_BLASTN_ENV_VAR = 'BLASTN_EXE'
_CALC_DIR = lambda f: os.path.join('run_dir', f)

# blastn parameters
_QUERY_FA = _CALC_DIR('query.fa')
# # Version 1: Blast referent SSC ends
# _GAP_PENALTY = '-gapopen 1 -gapextend 1 -penalty -2'
# Version 2: Blast referent IRa
_GAP_PENALTY = '-perc_identity 40'
_BLAST_PARAMS = f'-query {_QUERY_FA} {_GAP_PENALTY} -num_alignments 2 -max_hsps 2 -outfmt 5 -num_threads 1'

_install_instructions = """
BLAST+ is not installed!
Check web page https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
for installation instructions.

There are two ways for this script to locate executable to run:
 - environment variable {env_var} points to executable location,
 - or executable is called {exe} and placed on the PATH.
"""


# Note: it would be good that all scripts accept same format envs
def _find_exe(default_exe, env_var):
    exe = os.getenv(env_var, default_exe)
    if not shutil.which(exe):
        print(_install_instructions.format(exe=default_exe, env_var=env_var))
        raise ValueError(f'No BLAST+ installed! Tried {exe}')
    return exe


def _run_single(blastn_exe, input_fa, output_xml):
    cmd = f"{blastn_exe} -subject {input_fa} {_BLAST_PARAMS} > {output_xml}"
    print(f"Command: {cmd}")
    os.system(cmd)


def run(locale=True, threads=None):
    with open('finish.yml', 'r') as r:
        params = yaml.load(r, Loader=yaml.CLoader)
        calc_seq_idents = params['calc_seq_idents']
    if not calc_seq_idents:
        print('No sequences to blast!')
        return

    threads = threads or multiprocessing.cpu_count()
    blastn_exe = _find_exe(_DEFAULT_BLASTN_EXE_NAME, _BLASTN_ENV_VAR)
    outputs = []

    with ThreadPoolExecutor(max_workers=threads) as executor:
        for seq_ident in calc_seq_idents:
            outputs.append(_CALC_DIR(f'{seq_ident}.xml'))
            executor.submit(_run_single, blastn_exe, _CALC_DIR(f'{seq_ident}.fa'), outputs[-1])

    # Zip files
    if not locale:
        with ZipFile('output.zip', 'w') as output:
            for f in outputs:
                output.write('result.xml')


if __name__ == '__main__':
    import sys
    run(locale=False, threads=int(sys.argv[1]) if len(sys.argv) > 1 else None)
