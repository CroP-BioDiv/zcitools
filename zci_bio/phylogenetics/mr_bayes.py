import os
import random
from . import run_mr_bayes
from zci_bio.phylogenetics.steps import MrBayesStep, MrBayesSteps
from common_utils.file_utils import unzip_file, list_zip_files, write_yaml, read_yaml, \
    run_module_script, set_run_instructions, find_executable, silent_remove_file
from common_utils.exceptions import ZCItoolsValueError

_RESULT_PREFIX = 'result'

_SEED_DATA = """begin mrbayes;
    set seed={seed} swapseed={swapseed};
end;

"""

_NEXUS_DATA = """
begin mrbayes;
    set autoclose=yes nowarn=yes autoreplace=no;
    lset nst=6 rates=gamma;
    mcmcp ngen={ngen} printfreq={printfreq} samplefreq={samplefreq} nchains={nchains}
          savebrlens=yes filename={filename_prefix};
    mcmc;
    sumt filename={filename_prefix} {burnin} contype=halfcompat;
end;

"""

_NEXUS_DATA_PARTS = """
begin mrbayes;
    set autoclose=yes nowarn=yes autoreplace=no;
{partitions}
    lset applyto=(all) nst=6 rates=gamma;
    unlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all);
    prset applyto=(all) ratepr=variable;
    prset applyto=(all) statefreqpr=dirichlet(1,1,1,1);
    mcmcp ngen={ngen} printfreq={printfreq} samplefreq={samplefreq} nchains={nchains}
          savebrlens=yes filename={filename_prefix};
    mcmc;
    sumt filename={filename_prefix} {burnin} contype=halfcompat;
end;

"""

_instructions = """
Steps:
 - copy file calculate.zip onto server
 - unzip it
 - change directory to {step_name}
 - run script: {script_name}
    - to specify number of threads to use run: {script_name} <num_threads>
      default is number of cores.
 - copy file output.zip back into project's step directory {step_name}
 - run zcit command: zcit.py finish {step_name}

Notes:
 - MrBayes executable (mr_bayes) should be on the PATH or
   environment variable MR_BAYES_EXE should point to it.
 - It is good to use command screen for running the script.
   screen -dm "python3 {script_name}"
"""


def _copy_alignment_file(align_step, in_step, files_to_proc, args, partitions_obj):
    # ToDo: Nexus i dodati nes na kraj!
    a_f = in_step.step_file('alignment.nex')
    with open(a_f, 'a') as output:
        read_nexus = False
        with open(align_step.get_nexus_file(), 'r') as r:
            # Write header
            line = r.readline()
            assert line.startswith('#NEXUS')
            output.write(line)
            # Seed description
            output.write(_SEED_DATA.format(seed=random.randint(1000, 10000000), swapseed=random.randint(1000, 10000000)))
            # Nexus file content
            output.write(r.read())

        # Add MrBayes data
        # Printing
        ngen = args.ngen
        printfreq = str(max(1, ((ngen // 10000) if ngen > 1000000 else (ngen // 1000))))  # Of type str
        printfreq = int(printfreq[0] + ('0' * (len(printfreq) - 1)))                      # Round it
        # Burnin
        if args.burnin:
            brn = f'relburnin=no burnin={args.burnin}'
        elif (f := args.burninfrac) and 0 <= f < 1:
            brn = f'relburnin=yes burninfrac={f}'
        else:
            brn = ''  # Default is burninfrac=0.25
        # ToDo: Check or set samplefreq?
        params = dict(ngen=ngen, printfreq=printfreq, samplefreq=args.samplefreq, nchains=args.nchains, burnin=brn,
                      filename_prefix=_RESULT_PREFIX)
        if (partitions := partitions_obj.create_mrbayes_partitions(align_step, a_f)):
            output.write(_NEXUS_DATA_PARTS.format(partitions=partitions, **params))
        else:
            output.write(_NEXUS_DATA.format(**params))

    files_to_proc.append(dict(filename=a_f, result_prefix=in_step.step_file(_RESULT_PREFIX),
                              short=align_step.is_short(), nchains=args.nchains))


def create_mr_bayes_data(step_data, alignment_step, args, partitions_obj, run_threads):
    # List of dicts with attrs: filename, short
    # This data is used to optimize calculation
    # ToDo: almost the same as raxml.py. Differs in class types, _copy_alignment_file() and file formats
    files_to_proc = []

    if alignment_step._IS_COLLECTION:
        step = MrBayesSteps(alignment_step.project, step_data, remove_data=True)
        for align_step in alignment_step.step_objects():
            substep = step.create_substep(align_step.get_local_name())
            substep.set_sequences(align_step.all_sequences())
            substep.seq_sequence_type(align_step.get_sequence_type())
            _copy_alignment_file(align_step, substep, files_to_proc, args, partitions_obj)
            #
            substep.save(completed=False)
        if args.num_runs and args.num_runs > 1:
            print('Warning: number of runs for collection of alignments is not supported.')
    else:
        if args.num_runs and args.num_runs > 1:
            step = MrBayesSteps(alignment_step.project, step_data, remove_data=True)
            for run_idx in range(args.num_runs):
                substep = step.create_substep(f'RUN_{run_idx + 1}')
                substep.set_sequences(alignment_step.all_sequences())
                substep.seq_sequence_type(alignment_step.get_sequence_type())
                # ToDo: make symbolic links?
                _copy_alignment_file(alignment_step, substep, files_to_proc, args, partitions_obj)
                #
                substep.save(completed=False)
        else:
            step = MrBayesStep(alignment_step.project, step_data, remove_data=True)
            step.set_sequences(alignment_step.all_sequences())
            step.seq_sequence_type(alignment_step.get_sequence_type())
            _copy_alignment_file(alignment_step, step, files_to_proc, args, partitions_obj)

    # Store files desc
    files_to_zip = [d['filename'] for d in files_to_proc]  # files to zip
    # Remove step directory from files since run script is called from step directory
    for d in files_to_proc:
        d['filename'] = step.strip_step_dir(d['filename'])
        d['result_prefix'] = step.strip_step_dir(d['result_prefix'])
    finish_f = step.step_file('finish.yml')
    write_yaml(files_to_proc, finish_f)

    # Stores description.yml
    step.save(completed=bool(run_threads))

    if run_threads:
        run_module_script(run_mr_bayes, step, threads=run_threads, use_mpi=(not args.no_mpi))
    else:
        files_to_zip.append(finish_f)
        set_run_instructions(run_mr_bayes, step, files_to_zip, _instructions)
    #
    return step


def finish_mr_bayes_data(step_obj):
    output_f = step_obj.step_file('output.zip')
    if not os.path.isfile(output_f):
        raise ZCItoolsValueError('No calculation output file output.zip!')

    allowed_files = set(_RESULT_PREFIX + ext for ext in (
        '.ckp', '.con.tre', '.parts', '.run1.p', '.run1.t', '.run2.p', '.run2.t', '.tstat', '.vstat'))

    # Check are all file MrBayes outputs
    dirs = set(os.path.dirname(d['filename']) for d in read_yaml(step_obj.step_file('finish.yml')))
    for z_file in list_zip_files(output_f):
        parts = z_file.split('/')  # ZipFile uses '/' as separator
        _dir = '' if len(parts) == 1 else os.sep.join(parts[:-1])
        if _dir not in dirs:
            raise ZCItoolsValueError(f'Output contains file(s) in not step directory ({_dir})!')

        if parts[-1] not in allowed_files:
            raise ZCItoolsValueError(f'Not MrBayes output file(s)found in the output ({parts[-1]})!')

    # Unzip data
    unzip_file(output_f, step_obj.directory)

    step_obj._check_data()
    step_obj.save(create=False)


# Helper methods
def run_tracer(step, exe_location):
    exe = find_executable('tracer', exe_location)  # First check is there exe
    os.system(f"{exe} {' '.join(step.get_p_files())}")


def run_consense(step, exe_location):
    from subprocess import Popen, PIPE
    from .utils import read_tree, read_trees

    exe = find_executable('consense', exe_location)  # First check is there exe

    # Find burnin values. Check _copy_alignment_file() method.
    # Note: MrBayes exports also tree on gen 0. That is why index+1 is used.
    args = step.get_command_args()
    if (burnin := args['burnin']):
        _take_m = lambda _l: _l[burnin + 1:]
    else:
        burninfrac = args['burninfrac']
        if burninfrac is None:
            burninfrac = 0.25  # Default is burninfrac=0.25
        _take_m = lambda _l: _l[int(len(_l) * burninfrac) + 1:]

    # Extract trees
    trees = []
    for t_file in step.get_t_files():
        print(f'Reading file {t_file}')
        trees += _take_m(read_trees(t_file, format='nexus'))

    o_f = step.step_file('consense_input.tre')
    print(f'Writting file {o_f}')
    with open(o_f, 'w') as _c_out:
        for t in trees:
            _c_out.write(t.format('newick'))

    # Run consense
    # Note: If presented consense will ask for more things :-/
    silent_remove_file(outfile := step.step_file('outfile'))
    silent_remove_file(outtree := step.step_file('outtree'))
    p = Popen(exe, stdin=PIPE, stdout=PIPE, close_fds=True, cwd=os.path.abspath(step.directory))
    p.stdin.write(b'consense_input.tre\n')
    p.stdin.flush()
    p.stdin.write(b'Y\n')  # sends Enter into process
    p.stdin.flush()
    p.communicate()

    # Move consense output files
    if os.path.isfile(outfile):
        os.rename(outfile, step.step_file('consense_outfile'))
    if os.path.isfile(outtree):
        os.rename(outtree, step.step_file('consense_output.tre'))
