# Reproduction: Towards the Well-Tempered Chloroplast DNA Sequences

Research is implemented as a [workflow](project.md#workflow) named `chloroplast_normalization`
with required parameters family and outgroup, and optional parameter max_update_date.

**Step 1. main directory (optional)**

Create a directory for a whole research. Make it working directory.


**Step 2. family researches/directories**

Create a research directory for each family. This creates research directories and initialize them with provided data.

```
zcit.py init Apiaceae -w chloroplast_normalization -p 'family=Apiaceae;outgroup=NC_041648;max_update_date=2020-12-31'
zcit.py init Asteraceae -w chloroplast_normalization -p 'family=Asteraceae;outgroup=NC_040857;max_update_date=2020-12-31'
zcit.py init Campanulaceae -w chloroplast_normalization -p 'family=Campanulaceae;outgroup=NC_046571;max_update_date=2020-12-31'
zcit.py init Lamiaceae -w chloroplast_normalization -p 'family=Lamiaceae;outgroup=NC_046852;max_update_date=2020-12-31'
zcit.py init Rosaceae -w chloroplast_normalization -p 'family=Rosaceae;outgroup=NC_046061;max_update_date=2020-12-31'
```


**Step 3. For each family**

Make family research directory working one. Run workflow in it with command `zcit.py workflow run`. This will run 3 steps:

* fetch list of chloroplast accessions from NCBI (step 01_chloroplast_list),
* fetch cpDNA sequences (step 02_seqs),
* create step for GeSeq annotation of fetched sequences (step 03_GeSeq).

[GeSeq (Annotation of Organellar Genomes)](https://chlorobox.mpimp-golm.mpg.de/geseq.html) is web-based tool.
Data has to be posted and results downloaded 'by the hand'.
Follow printed instructions, or same one from the file `02_GeSeq/INSTRUCTIONS.txt`.

After downloading GeSeq results, make `zcit.py workflow run`.
Note: this runs steps so that [demanding calculations](project.md#running_calculations) have to be calculated on a dedicated server.
In case of that, please consult [demanding calculations](project.md#running_calculations). To sun steps locally run command `zcit.py workflow run -r`.


**Notes:**

* In each stage of the calculation it is possible to visually research status with command `zcit.py workflow graph`.
* Workflow uses external programs for alignment and phylogenetic analyses. There are 3 programs that have to be installed and the `PATH` on computer where calculations will be performed (locally or on dedicated server). In case some program is missing, calculation will stop with message what should be installed and set. Programs needed are:
  * [MAFFT](https://mafft.cbrc.jp/alignment/software/). Executable should be named `mafft`, or alternative name set with `MAFFT_EXE` environmnet variable.
  * [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/). Executable should be named `raxml_threads`, or alternative name set with `RAXML_EXE` environmnet variable.
  * [MrBayes](https://github.com/CroP-BioDiv/MrBayes). Executable should be named `mr_bayes` or `mr_bayes_mpi` in case of MPI  version, or alternative names set with `MR_BAYES_EXE` and `MR_BAYES_EXE` environmnet variables.
    * **Note:** MrBayes hangs on some variation of input data ([issue](https://github.com/NBISweden/MrBayes/issues/230)). Please use fixed [version](https://github.com/CroP-BioDiv/MrBayes). Original version is [here](https://nbisweden.github.io/MrBayes/index.html).
