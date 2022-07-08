# How to reproduce results:<br/>Chloroplast genome annotation tools: Challenges and recommendations

Research is implemented as a [workflow](project.md#workflow) named `irs_statistics`
with parameters for annotation methods to check and clades (sequences) to use.

**Step 1. main directory (optional)**

Create a directory for a whole research. Make it working directory.

**Step 2. dataset directories**

Create a dataset directories. Researches were done on datasets IR, IRL and union dataset, with chloroplast sequences published before 2021-12-31.

Create project for dataset IR:
```
zcit.py init dataset_IR -w irs_statistics -p "taxons=asterids,rosids;methods=research;plastids=1;taxa_ranks=order,family,genus;remove_clades=IRL clade,Erodium,Monsonia,Cytinus,Pilostyles,Putranjivaceae,Lophopyxidaceae,Ericaceae;max_update_date=2021-12-31"
```

Create project for dataset IRL:
```
zcit.py init dataset_IRL -w irs_statistics -p "taxons=IRL clade,Erodium,Monsonia,Cytinus,Pilostyles,Putranjivaceae,Lophopyxidaceae,Ericaceae;methods=research;plastids=1;taxa_ranks=order,family,genus;taxa_names=asterids,rosids;remove_clades=;max_update_date=2021-12-31"
```

Create project for union dataset:
```
zcit.py init All -w irs_statistics -p "taxons=asterids,rosids;methods=research;plastids=1;taxa_ranks=order,family,genus;max_update_date=2021-12-31"
```

**Step 3. For each dataset**

Make dataset project directory working one. Run workflow in it with command `zcit.py workflow run`. This will run 4 steps:

* fetch list of chloroplast accessions from NCBI (step 01_chloroplast_list),
* fetch cpDNA sequences (step 02_seqs),
* fetch additional data for methods that need them (02_chloe, 03_chloe, 02_ge_seq, 03_ge_seq),
* run all annotation methods on all sequences and store results in subdirectory 04_stats and excel files (chloroplast_irs_analysis.xlsx, workflow_summary.xlsx).

**Notes:**
* [GeSeq (Annotation of Organellar Genomes)](https://chlorobox.mpimp-golm.mpg.de/geseq.html) is a web-based tool.
  Data has to be posted and results downloaded 'by the hand'.
  Follow printed instructions, or same one from the file `02_ge_seq/INSTRUCTIONS.txt`.
  After downloading GeSeq results, make `zcit.py workflow run`.
* In each stage of the calculation it is possible to visually research status with command `zcit.py workflow graph`.
* Workflow uses external annotation tool executable. These has to be installed and on the `PATH`. Check
  [stand-alone wrapper](https://github.com/CroP-BioDiv/irs_wrappers) project for details.
