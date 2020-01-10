import os.path
from common_utils.exceptions import ZCItoolsValueError
from common_utils.file_utils import ensure_directory, copy_file, link_file, write_str_in_file, write_yaml, \
    run_module_script, set_run_instructions
from .steps import QTLCartStep
from . import run_qtl_cart_perm

_instructions = """
Steps:
 - copy file calculate.zip onto server
 - unzip it
 - change directory to {step_name}
 - run script: python3 {script_name}
    - to specify number of threads to use run: python3 {script_name} <num_threads>
      default is number of cores.
 - copy file output.zip back into project's step directory {step_name}
 - run zcit command: zcit.py finish {step_name}

Notes:
 - Clustal Omega executable (clustalo) should be on the PATH or
   environment variable CLUSTAL_OMEGA_EXE should point to it.
 - It is good to use command screen for running the script.
   screen -dm "python3 {script_name}"
"""

# Copied from some project!
_qtlcart_rc = """#
#  This is the parameter file for QTLcartographer.
#
#        Some basic parameters
#
-chrom               7  # (Number of chromosomes: Rmap)
-mark              302  # (Average number of markers per chromosome: Rmap)
-vmark         67.8131  # (Standard Deviation in the number of markers/chromosome: Rmap)
-dist           1.1072  # (Average intermarker distance, in cM: Rmap)
-vdist          1.3638  # (Standard Deviation of intermarker distance: Rmap)
-qnum                9  # (Average number of QTLs per trait: Rqtl)
-dom                 1  # (type of dominance: Rqtl)
-beta           2.0000  # (Value of additive beta: Rqtl)
-beta1          2.0000  # (Value of dominance beta1: Rqtl)
-beta2          2.0000  # (Value of dominance beta2: Rqtl)
-traits   {num_traits}  # (Number of quantitative traits: Rqtl)
-whichtrait    {trait}  # (The trait to be analyzed: Zmapqtl)
-Herit          0.5000  # (Heritability: Rcross)
-cross               5  # (Type of cross: Rcross)
-nn                111  # (Sample size: Rcross)
-Model               6  # (Model for analysis: Zmapqtl)
-walk           1.0000  # (Walking speed along chromosomes, in cM: Zmapqtl)
-window        10.0000  # (Window width, in cM: Zmapqtl)
-nbp                5  # (Number of background parameters: Zmapqtl)
-mapfunc             1  # (Map Function)
-gout                1  # (Output flag: Rmap)
-Rmode               0  # (Simulation Mode: Rmap)
-emaplink       0.2500  # (Size of test for linkage: Emap)
-emapseg        0.0100  # (Size of test for segregation: Emap)
-emapmeth           10  #  (Emap method flag)
-emapobj             0  #  (Emap objective function, )
-lodflag             0  # (1 => LOD scores in Preplot, Eqtl)
-srF            0.1000  # (p(Fin): SRmapqtl)
-srB            0.1000  # (p(Fout): SRmapqtl)
-srM                 2  # (Regression type: SRmapqtl)
-srupper           100  # (Maximun number of steps in stepwise regression: SRmapqtl)
-size           0.0500  # (Size, or alpha)
-siglevel       3.8400  # (Size, or alpha)
-LRmax          0.0000  # Maximum LR or LOD score
-Effectmax        0.000000  # Maximum additive effect 
-maxqtl             19  # (Maximum QTL to fit in MImapqtl)
-maxepis            19  # (Maximum Epistatic effect to fit in MImapqtl)
-lodmim         0.0000  # (LOD for adding/dropping parameters in MImapqtl)
-ic                  1  # (Code for IC criterion in MImapqtl)
-phasemim            0  # (Phase of analysis MImapqtl)
-Hypothesis         10  # (Hypothesis test to do/process)
#
# These are the default filenames...
#
-stem                        qtlcart  # (Stem for filenames)
-Thecross                        RI1  # (Type of cross)
-error                   qtlcart.log  # (Log and error file: All)
-map                         tmp.00m  # (Rmap ouput file, linkage map)
-mapin                      temp.00m  # (Rmap input file)
-qtl                     qtlcart.qtl  # (Rqtl output file)
-ifile                       tmp.00c  # (Rcross output file, individual data file)
-iinfile                    temp.00c  # (Rcross input file, individual data file)
-lrfile                   qtlcart.lr  # (Results of Linear Regression analysis)
-srfile                      tmp.00r  # (Results of Stepwise Regression analysis)
-qstat                   qtlcart.qst  # (Results of Qstats)
-zfile       results_T{trait:02}.txt  # (Results of (Composite) Interval Mapping Analysis)
-eqtlfile                   qtlcart.eqt  # (Eqtl output file)
-mimfile                 qtlcart.mim  # (Results of Multiple Interval Mapping Analysis)
-mqtfile                qtlcarti.mqt  # (Multiple Interval Mapping Analysis Model)
-bayesfile               qtlcart.bys  # (Results of Bayesian Analysis)
-mrinput                qtlcart.zr  # (Output of JZmapqtl for input to MultiRegress)
-mroutput                qtlcart.mr  # (Results of MultiRegress)
"""


def create_permutations(project, step_data, raw_file, permutations, num_traits=None, run=False):
    # Check input files
    map_file = raw_file.replace('.raw', '.map')
    data_dir, base_raw_file = os.path.split(raw_file)
    tmp_files = ('tmp.00m', 'tmp.00c', 'tmp.00r')
    for mf in (raw_file, map_file):
        if not os.path.isfile(mf):
            raise ZCItoolsValueError(f"Input MapMaker file {mf} doesn't exist!")
    for qf in tmp_files:
        f = os.path.join(data_dir, qf)
        if not os.path.isfile(f):
            raise ZCItoolsValueError(f"Input Windows QTL Cartographer file {qf} doesn't exist!")

    #
    step = QTLCartStep(project, step_data, remove_data=True)

    # Copy input files
    for qf in tmp_files:
        copy_file(os.path.join(data_dir, qf), step.step_file(qf))

    files_to_zip = []
    # Create trait directories
    # ToDo: find max traits and fix it/set default
    assert num_traits and num_traits > 0, num_traits
    trait_dirs = []
    for t_idx in range(1, num_traits + 1):
        trait_dirs.append(f'T{t_idx:02}')
        t_dir = step.step_file(trait_dirs[-1])
        ensure_directory(t_dir)
        # Create links to input files
        for qf in tmp_files:
            link_file(os.path.join('..', qf), os.path.join(t_dir, qf))
        #
        write_str_in_file(os.path.join(t_dir, 'qtlcart.rc'), _qtlcart_rc.format(trait=t_idx, num_traits=num_traits))

    write_yaml(dict(permutations=permutations, trait_dirs=trait_dirs), step.step_file('finish.yml'))

    # ToDo: files_to_zip = [d['filename'] for d in seq_files]  # files to zip

    # Run or set instructions
    if run:
        run_module_script(run_qtl_cart_perm, step)
    else:
        # files_to_zip.append(finish_f)
        set_run_instructions(run_qtl_cart_perm, step, files_to_zip, _instructions)

    #
    step.save()
    return step


def finish_permutations(step_obj):
    print('ToDo')
    pass
