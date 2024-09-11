# Case Report: T-Large Granular Lymphocyte-like Leukemia Following CAR T-Cell Therapy

This repository contains code used to produce the results in: Till Braun, Michael Rade, Maximilian Merz et al. T-Large Granular Lymphocyte-like Leukemia Following CAR T-Cell Therapy  ...

## Singularity

All scripts were developed in a Singularity image with Rstudio server. [See README](singularity/) in `./singularity/`. 

## Reproduction

``` sh
$ bash reproduce.sh
```

The script must be run in the base path (./) of this repository. Following the content of `reproduce.sh`:

``` sh
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# scRNA-Seq, VDJ-Seq, ADT-Seq
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Rscript code/multi_omics/01_cellranger.R
Rscript code/multi_omics/02_cellranger_qc.R
Rscript code/multi_omics/03_cellranger_to_seurat
Rscript code/multi_omics/04_qc_prep.R
Rscript code/multi_omics/05_anno_1.R
Rscript code/multi_omics/05_anno_2.R
Rscript code/multi_omics/06_integration.R
Rscript code/multi_omics/07_infercnv.R

# Main Figure 2
Rscript code/publication/main/fig_02.R

# Supplemental figures for Figure 2
Rscript code/publication/supps/fig_02_supps.R

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# WGS
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
...

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# scRNA-seq mutation
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bash code/multi_omics/08_scRNA_mutations/01_indexing.sh
bash code/multi_omics/08_scRNA_mutations/02_pileup.sh
Rscript code/multi_omics/08_scRNA_mutations/03_summarize_data.R
Rscript publication/figure_scripts/supps/fig_03_supp_scRNAseq_mutations.R

```
