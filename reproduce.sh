#!/bin/bash

# The script must be executed in the base path (./) of this repo

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# scRNA-Seq, VDJ-Seq, ADT-Seq
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Rscript code/multi_omics/01_cellranger.R
# Rscript code/multi_omics/02_cellranger_qc.R
Rscript code/multi_omics/03_cellranger_to_seurat.R
Rscript code/multi_omics/04_qc_prep.R
Rscript code/multi_omics/05_anno_1.R
Rscript code/multi_omics/05_anno_2.R
Rscript code/multi_omics/06_integration.R
Rscript code/multi_omics/07_infercnv.R


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# WGS
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# CTAT-VirusIntegrationFinder
# bash code/wgs/ctat_vif/01_prep_db.sh
# bash code/wgs/ctat_vif/02_run_vif.sh

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# scRNA-seq mutation
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bash code/multi_omics/08_scRNA_mutations/01_indexing.sh
bash code/multi_omics/08_scRNA_mutations/02_pileup.sh
Rscript code/multi_omics/08_scRNA_mutations/03_summarize_data.R
Rscript publication/figure_scripts/supps/fig_03_supp_scRNAseq_mutations.R

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Main Figures
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Rscript code/publication/main/fig_02.R

