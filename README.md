# Case Report: T-Large Granular Lymphocyte-like Leukemia Following CAR T-Cell Therapy

This repository contains code used to produce the results in: Till Braun, Michael Rade, Maximilian Merz et al. T-Large Granular Lymphocyte-like Leukemia Following CAR T-Cell Therapy  ...

## Singularity

All scripts for the single-cell Analysis were developed in a Singularity image with Rstudio server. [See README](singularity/) in `./singularity/`. Viral integration and variant calling analysis was performed using the nf-core pipelines [https://nf-co.re/viralintegration/0.1.1/](https://nf-co.re/viralintegration/0.1.1/) and [https://nf-co.re/sarek/3.4.3/](https://nf-co.re/sarek/3.4.3/). A workflow for the WGS data including variant calling, reverse engineering of the cilta-cel vector sequence, and STAR alignment for integration site analysis is shown in `code/wgs/`

## Reproduction

``` sh
$ bash reproduce.sh
```

The script must be run in the base path (./) of this repository. Following the content of `reproduce.sh`:

``` sh
#!/bin/bash

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# scRNA-Seq, VDJ-Seq, ADT-Seq:
# All analysis performed with the Singularity R studio image
# The script must be executed in the base path (./) of this repo
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Rscript code/multi_omics/01_cellranger.R
Rscript code/multi_omics/02_cellranger_qc.R
Rscript code/multi_omics/03_cellranger_to_seurat.R
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
# WGS:
# STAR and Samtools are required
# https://nf-co.re/viralintegration/0.1.1/
# https://nf-co.re/sarek/3.4.3/
# see also code/wgs/
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# CTAT-VirusIntegrationFinder
# bash code/wgs/ctat_vif/01_prep_db.sh
# bash code/wgs/ctat_vif/02_run_vif.sh

# STAR mappings for integration site analyis
bash code/wgs/STAR/01_STARindex.sh
bash code/wgs/STAR/02_arribaSetting_multiMapper_S1.sh
bash code/wgs/STAR/03_arribaSetting_multiMapper_S2.sh
bash code/wgs/STAR/04_arribaSetting_multiMapper_S3.sh
bash code/wgs/STAR/05_arribaSetting_multiMapper_S4.sh
bash code/wgs/STAR/06_NH1_HI1_AS100_filter_placeholder.sh
Rscript code/wgs/STAR/07_STAR_IS_candidates.R #generates figure S17

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# scRNA-seq mutation:
# For "01_indexing.sh" and "02_pileup.sh samtools" and BCFtools are required
# The R scripts were run in the Singularity R Studio image
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Samtools is required
bash code/multi_omics/08_scRNA_mutations/01_indexing.sh
bash code/multi_omics/08_scRNA_mutations/02_pileup.sh
Rscript code/multi_omics/08_scRNA_mutations/03_summarize_data.R
Rscript publication/figure_scripts/supps/fig_03_supp_scRNAseq_mutations.R

```
