library(conflicted)
library(tidyverse)
conflict_prefer("filter", winner = "dplyr")
library(stringr)
library(stringi)
library(readr)
library(vcfR)

# extracted names of annotation columns
CSQ.names = str_split("Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|gnomADe_AF|gnomADe_AFR_AF|gnomADe_AMR_AF|gnomADe_ASJ_AF|gnomADe_EAS_AF|gnomADe_FIN_AF|gnomADe_NFE_AF|gnomADe_OTH_AF|gnomADe_SAS_AF|gnomADg_AF|gnomADg_AFR_AF|gnomADg_AMI_AF|gnomADg_AMR_AF|gnomADg_ASJ_AF|gnomADg_EAS_AF|gnomADg_FIN_AF|gnomADg_MID_AF|gnomADg_NFE_AF|gnomADg_OTH_AF|gnomADg_SAS_AF|MAX_AF|MAX_AF_POPS|FREQS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS", pattern = fixed('|'))[[1]]

# genes of interest for filtering
genes_of_interest = list(
  Huang_2024_PTCL_frequently_mutated = c("APC","APC2","ARID1A","ARID1B","ARID2","ASXL3","ATM","BCOR","BIRC3","BIRC6","CARD11","CCND3","CD58","CDKN2A","CHD8","CHEK2","CIC","CIITA","CREBBP","DNMT3A","EP300","FYN","HLA-A","HLA-B","IDH2","IKBKB","ITPKB","ITPR3","JAK1","JAK2","JAK3","JMY","KDM6B","KMT2A","KMT2C","KMT2D","LRP1B","MEF2A","MGA","MSH2","MSH3","MSH6","MTOR","NCOR2","NF1","NOTCH1","NOTCH2","NOTCH3","PDCD1","PHLPP1","PIK3R1","PLCG1","PLCG2","PMS1","PTEN","PTPN13","PTPRC","PTPRD","PTPRS","REV3L","RHOA","SALL3","SETBP1","SETD1B","SETD2","SMARCA2","SMARCA4","SOCS1","SPEN","STAT3","STAT5B","TAL1","TET1","TET2","TET3","TP53","TRRAP","TSC2","VAV1","YEATS2","YTHDF2","ZEB1"),
  Cheon_2022_TLGLL_putative_driver_mutations = c("STAT3","KMT2D","TNFAIP3","KDM6A","ABCC9","PIK3R1","TET2","PCDHA11","SLC6A15","SULF1","ARHGAP25","DDX59","DNMT3A","FAS","STAT5B"),
  Lone_2018_PTCL_frequently_mutated = c("SETD2","TET3","TET2","IDH2","DNMT3A","PLCG1","CD28","FYN","CARD11","RHOA","STAT5B","STAT3","INO80","PIK3CD","ARID1B","UBR5","MLL3","MLL2","DDX3X","IRF4","CCR4","JAK1","PRDM1","CD58","PHIP","MYD88","PIM1","CD79B","MYC","ID3","TP53","GNA13","EZH2","SMARCA2","CCND1","PTPN1","B2M","TNFAIP3")
)
all_genes_of_interest = unique(unlist(genes_of_interest, use.names = F))

# chromosome order function
example.meta = read.vcfR("../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARneg_CD5pos_merged_vs_Apherese_CARneg_CD5pos.dragmap.mutect2_VEP.filtered.vcf.gz")@meta
chrs = example.meta[startsWith(example.meta, "##contig")]
chrs = chrs[c(1:25,length(chrs))]# all large assemblies + carvykti (although probably not necessary)
trim_chrom <- function(x){substr(x,4,str_length(x)-1)}
chrs = str_extract(chrs, pattern = "ID=.*,") %>% trim_chrom()
CHROM.order = function(x)match(x,chrs)

prep_df <- function(inpath,
                    gtcols = NULL,
                    drop_GT = F){
  df = read.vcfR(inpath, limit = 1e10, cols = gtcols) %>%
    vcfR2tidy(single_frame = T)
  df = df$dat
  
  if (drop_GT) {
    df = df %>%
      select(!(Indiv:gt_GT_alleles)) %>%
      distinct() %>%
      mutate(ID = paste(CHROM, POS, paste0(REF, ">", ALT), sep = "_"))
  } else {
    df = df %>%
      distinct() %>%
      mutate(ID = paste(CHROM, POS, paste0(REF, ">", ALT), sep = "_"))
  }
  
  return(df)
}

# filter functions
filter_vcf_inputs_by_relevance <- function(df_infile,
                                      rarity_threshold=0.001){
  
  csq = df_infile %>%
    select(ID, CSQ) %>%
    dplyr::filter(!is.na(CSQ)) %>% 
    separate_longer_delim(cols = CSQ, delim = ",") %>%
    mutate(CSQ = stri_trim(CSQ, side = "both", pattern = "[^()]")) %>%
    separate_wider_delim(cols = CSQ, delim = "|", names = CSQ.names)
    
  df = df_infile %>% select(!CSQ)
  
  IDs.missing.anno = setdiff(unique(df$ID), unique(csq$ID))
  if (length(IDs.missing.anno)>0) {
    warning(paste("Annotation seemingly missing in", length(IDs.missing.anno), "variants. Skipping them during filtering."))
  }
  
  csq.filtered = csq %>%
    mutate(  MAX_AF = as.double(MAX_AF),
             Rare = (is.na(MAX_AF) | (MAX_AF < rarity_threshold)),
             ImpactHighMed = (IMPACT=="HIGH" | IMPACT=="MODERATE")) %>%
    filter(  str_detect(SIFT, pattern = fixed("deleterious")) |  # predicted deleterious by SIFT
             str_detect(PolyPhen, pattern = fixed("damaging")) | # or polyphen
             ((SYMBOL %in% all_genes_of_interest) & (Rare | ImpactHighMed)) | # or possibly impacting genes of interest
             (str_detect(CLIN_SIG, fixed("pathogenic")) | str_detect(CLIN_SIG, fixed("risk_factor")) | str_detect(CLIN_SIG, fixed("drug_response"))) | # annotated possibly pathogenic or risk factor or drug response associated
             ( (IMPACT=="HIGH" | (IMPACT=="MODERATE" & Rare)) & str_detect(CLIN_SIG, fixed("benign"), negate =T)) # further high-impact or rare moderate-impact variants that are not already known to be benign (note: the line above rescues all lines with conflicting benign&pathogenic annotation)
    )
  
  csq.keepflag = csq.filtered %>%
    select(ID, ) %>%
    mutate(KEEP = T)
  
  filt.df = left_join(df, csq.keepflag, by = "ID") %>%
    filter(KEEP)
  
  return(filt.df)
}

################################################################################

unfiltered.IDs = list()
filtered.IDs = list()

# first: for each sample, read in all calls from both aligners -> note unfiltered and filtered IDs for Venn diagrams

## CARneg_CD5pos_merged ########################################################
vcf.paths.somatic.CARneg_CD5pos_merged = list(
  mutect.filtered.dragmap = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARneg_CD5pos_merged_vs_Apherese_CARneg_CD5pos.dragmap.mutect2_VEP.filtered.vcf.gz",
  strelka.indels.filtered.dragmap = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARneg_CD5pos_merged_vs_Apherese_CARneg_CD5pos.dragmap.strelka.somatic_indels_VEP.filtered.vcf.gz",
  strelka.snvs.filtered.dragmap = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARneg_CD5pos_merged_vs_Apherese_CARneg_CD5pos.dragmap.strelka.somatic_snvs_VEP.filtered.vcf.gz",
  mutect.filtered.bwamem2 = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARneg_CD5pos_merged_vs_Apherese_CARneg_CD5pos.bwamem2.mutect2_VEP.filtered.vcf.gz",
  strelka.indels.filtered.bwamem2 = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARneg_CD5pos_merged_vs_Apherese_CARneg_CD5pos.bwamem2.strelka.somatic_indels_VEP.filtered.vcf.gz",
  strelka.snvs.filtered.bwamem2 = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARneg_CD5pos_merged_vs_Apherese_CARneg_CD5pos.bwamem2.strelka.somatic_snvs_VEP.filtered.vcf.gz"
)
dfs.somatic.CARneg_CD5pos_merged = lapply(vcf.paths.somatic.CARneg_CD5pos_merged, function(inpath)prep_df(inpath, gtcols = 9, drop_GT = T))
unfiltered.IDs[["CARneg_CD5pos_merged"]] = lapply(dfs.somatic.CARneg_CD5pos_merged, function(df)df$ID)
dfs.filtered.somatic.CARneg_CD5pos_merged = lapply(dfs.somatic.CARneg_CD5pos_merged, function(in_df){
    filter_vcf_inputs_by_relevance(in_df) %>%
    select(CHROM:FILTER) %>% 
    mutate(INFO = "DP=50") # dummy
  })
filtered.IDs[["CARneg_CD5pos_merged"]] = lapply(dfs.filtered.somatic.CARneg_CD5pos_merged, function(df)df$ID)
somatic.CARneg_CD5pos_merged.all = list(
  mutect.filtered = inner_join(dfs.filtered.somatic.CARneg_CD5pos_merged$mutect.filtered.dragmap,
                               dfs.filtered.somatic.CARneg_CD5pos_merged$mutect.filtered.bwamem2),
  strelka.indels.filtered = semi_join(dfs.filtered.somatic.CARneg_CD5pos_merged$strelka.indels.filtered.dragmap,
                                       dfs.filtered.somatic.CARneg_CD5pos_merged$strelka.indels.filtered.bwamem2, by = "ID"),
  strelka.snvs.filtered = semi_join(dfs.filtered.somatic.CARneg_CD5pos_merged$strelka.snvs.filtered.dragmap,
                                     dfs.filtered.somatic.CARneg_CD5pos_merged$strelka.snvs.filtered.bwamem2, by = "ID")
) %>% do.call(what=rbind) %>% distinct()


## CARpos_CD5dim_merged ########################################################
vcf.paths.somatic.CARpos_CD5dim_merged = list(
  mutect.filtered.dragmap = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARpos_CD5dim_merged_vs_Apherese_CARneg_CD5pos.dragmap.mutect2_VEP.filtered.vcf.gz",
  strelka.indels.filtered.dragmap = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARpos_CD5dim_merged_vs_Apherese_CARneg_CD5pos.dragmap.strelka.somatic_indels_VEP.filtered.vcf.gz",
  strelka.snvs.filtered.dragmap = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARpos_CD5dim_merged_vs_Apherese_CARneg_CD5pos.dragmap.strelka.somatic_snvs_VEP.filtered.vcf.gz",
  mutect.filtered.bwamem2 = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARpos_CD5dim_merged_vs_Apherese_CARneg_CD5pos.bwamem2.mutect2_VEP.filtered.vcf.gz",
  strelka.indels.filtered.bwamem2 = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARpos_CD5dim_merged_vs_Apherese_CARneg_CD5pos.bwamem2.strelka.somatic_indels_VEP.filtered.vcf.gz",
  strelka.snvs.filtered.bwamem2 = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARpos_CD5dim_merged_vs_Apherese_CARneg_CD5pos.bwamem2.strelka.somatic_snvs_VEP.filtered.vcf.gz"
)
dfs.somatic.CARpos_CD5dim_merged = lapply(vcf.paths.somatic.CARpos_CD5dim_merged, function(inpath)prep_df(inpath, gtcols = 9, drop_GT = T))
unfiltered.IDs[["CARpos_CD5dim_merged"]] = lapply(dfs.somatic.CARpos_CD5dim_merged, function(df)df$ID)
dfs.filtered.somatic.CARpos_CD5dim_merged = lapply(dfs.somatic.CARpos_CD5dim_merged, function(in_df){
  filter_vcf_inputs_by_relevance(in_df) %>%
    select(CHROM:FILTER) %>% 
    mutate(INFO = "DP=50") # dummy
})
filtered.IDs[["CARpos_CD5dim_merged"]] = lapply(dfs.filtered.somatic.CARpos_CD5dim_merged, function(df)df$ID)
somatic.CARpos_CD5dim_merged.all = list(
  mutect.filtered = inner_join(dfs.filtered.somatic.CARpos_CD5dim_merged$mutect.filtered.dragmap,
                               dfs.filtered.somatic.CARpos_CD5dim_merged$mutect.filtered.bwamem2),
  strelka.indels.filtered = semi_join(dfs.filtered.somatic.CARpos_CD5dim_merged$strelka.indels.filtered.dragmap,
                                       dfs.filtered.somatic.CARpos_CD5dim_merged$strelka.indels.filtered.bwamem2, by = "ID"),
  strelka.snvs.filtered = semi_join(dfs.filtered.somatic.CARpos_CD5dim_merged$strelka.snvs.filtered.dragmap,
                                     dfs.filtered.somatic.CARpos_CD5dim_merged$strelka.snvs.filtered.bwamem2, by = "ID")
) %>% do.call(what=rbind) %>% distinct()


## CARpos_CD5pos_merged ########################################################
vcf.paths.somatic.CARpos_CD5pos_merged = list(
  mutect.filtered.dragmap = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARpos_CD5pos_merged_vs_Apherese_CARneg_CD5pos.dragmap.mutect2_VEP.filtered.vcf.gz",
  strelka.indels.filtered.dragmap = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARpos_CD5pos_merged_vs_Apherese_CARneg_CD5pos.dragmap.strelka.somatic_indels_VEP.filtered.vcf.gz",
  strelka.snvs.filtered.dragmap = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARpos_CD5pos_merged_vs_Apherese_CARneg_CD5pos.dragmap.strelka.somatic_snvs_VEP.filtered.vcf.gz",
  mutect.filtered.bwamem2 = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARpos_CD5pos_merged_vs_Apherese_CARneg_CD5pos.bwamem2.mutect2_VEP.filtered.vcf.gz",
  strelka.indels.filtered.bwamem2 = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARpos_CD5pos_merged_vs_Apherese_CARneg_CD5pos.bwamem2.strelka.somatic_indels_VEP.filtered.vcf.gz",
  strelka.snvs.filtered.bwamem2 = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARpos_CD5pos_merged_vs_Apherese_CARneg_CD5pos.bwamem2.strelka.somatic_snvs_VEP.filtered.vcf.gz"
)
dfs.somatic.CARpos_CD5pos_merged = lapply(vcf.paths.somatic.CARpos_CD5pos_merged, function(inpath)prep_df(inpath, gtcols = 9, drop_GT = T))
unfiltered.IDs[["CARpos_CD5pos_merged"]] = lapply(dfs.somatic.CARpos_CD5pos_merged, function(df)df$ID)
dfs.filtered.somatic.CARpos_CD5pos_merged = lapply(dfs.somatic.CARpos_CD5pos_merged, function(in_df){
  filter_vcf_inputs_by_relevance(in_df) %>%
    select(CHROM:FILTER) %>% 
    mutate(INFO = "DP=50") # dummy
})
filtered.IDs[["CARpos_CD5pos_merged"]] = lapply(dfs.filtered.somatic.CARpos_CD5pos_merged, function(df)df$ID)
somatic.CARpos_CD5pos_merged.all = list(
  mutect.filtered = inner_join(dfs.filtered.somatic.CARpos_CD5pos_merged$mutect.filtered.dragmap,
                               dfs.filtered.somatic.CARpos_CD5pos_merged$mutect.filtered.bwamem2),
  strelka.indels.filtered = semi_join(dfs.filtered.somatic.CARpos_CD5pos_merged$strelka.indels.filtered.dragmap,
                                       dfs.filtered.somatic.CARpos_CD5pos_merged$strelka.indels.filtered.bwamem2, by = "ID"),
  strelka.snvs.filtered = semi_join(dfs.filtered.somatic.CARpos_CD5pos_merged$strelka.snvs.filtered.dragmap,
                                     dfs.filtered.somatic.CARpos_CD5pos_merged$strelka.snvs.filtered.bwamem2, by = "ID")
) %>% do.call(what=rbind) %>% distinct()


## CARpos_CD5dim (Run 1) #######################################################
vcf.paths.somatic.CARpos_CD5dim = list(
  mutect.filtered.dragmap = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARpos_CD5dim_vs_Apherese_CARneg_CD5pos.dragmap.mutect2_VEP.filtered.vcf.gz",
  strelka.indels.filtered.dragmap = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARpos_CD5dim_vs_Apherese_CARneg_CD5pos.dragmap.strelka.somatic_indels_VEP.filtered.vcf.gz",
  strelka.snvs.filtered.dragmap = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARpos_CD5dim_vs_Apherese_CARneg_CD5pos.dragmap.strelka.somatic_snvs_VEP.filtered.vcf.gz",
  mutect.filtered.bwamem2 = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARpos_CD5dim_vs_Apherese_CARneg_CD5pos.bwamem2.mutect2_VEP.filtered.vcf.gz",
  strelka.indels.filtered.bwamem2 = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARpos_CD5dim_vs_Apherese_CARneg_CD5pos.bwamem2.strelka.somatic_indels_VEP.filtered.vcf.gz",
  strelka.snvs.filtered.bwamem2 = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARpos_CD5dim_vs_Apherese_CARneg_CD5pos.bwamem2.strelka.somatic_snvs_VEP.filtered.vcf.gz"
)
dfs.somatic.CARpos_CD5dim = lapply(vcf.paths.somatic.CARpos_CD5dim, function(inpath)prep_df(inpath, gtcols = 9, drop_GT = T))
unfiltered.IDs[["CARpos_CD5dim"]] = lapply(dfs.somatic.CARpos_CD5dim, function(df)df$ID)
dfs.filtered.somatic.CARpos_CD5dim = lapply(dfs.somatic.CARpos_CD5dim, function(in_df){
  filter_vcf_inputs_by_relevance(in_df) %>%
    select(CHROM:FILTER) %>% 
    mutate(INFO = "DP=50") # dummy
})
filtered.IDs[["CARpos_CD5dim"]] = lapply(dfs.filtered.somatic.CARpos_CD5dim, function(df)df$ID)
somatic.CARpos_CD5dim.all = list(
  mutect.filtered = inner_join(dfs.filtered.somatic.CARpos_CD5dim$mutect.filtered.dragmap,
                               dfs.filtered.somatic.CARpos_CD5dim$mutect.filtered.bwamem2),
  strelka.indels.filtered = semi_join(dfs.filtered.somatic.CARpos_CD5dim$strelka.indels.filtered.dragmap,
                                       dfs.filtered.somatic.CARpos_CD5dim$strelka.indels.filtered.bwamem2, by = "ID"),
  strelka.snvs.filtered = semi_join(dfs.filtered.somatic.CARpos_CD5dim$strelka.snvs.filtered.dragmap,
                                     dfs.filtered.somatic.CARpos_CD5dim$strelka.snvs.filtered.bwamem2, by = "ID")
) %>% do.call(what=rbind) %>% distinct()


## P2251_CARpos_CD5dim (Run 2) #################################################
vcf.paths.somatic.P2251_CARpos_CD5dim = list(
  mutect.filtered.dragmap = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/P2251_CARpos_CD5dim_vs_Apherese_CARneg_CD5pos.dragmap.mutect2_VEP.filtered.vcf.gz",
  strelka.indels.filtered.dragmap = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/P2251_CARpos_CD5dim_vs_Apherese_CARneg_CD5pos.dragmap.strelka.somatic_indels_VEP.filtered.vcf.gz",
  strelka.snvs.filtered.dragmap = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/P2251_CARpos_CD5dim_vs_Apherese_CARneg_CD5pos.dragmap.strelka.somatic_snvs_VEP.filtered.vcf.gz",
  mutect.filtered.bwamem2 = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/P2251_CARpos_CD5dim_vs_Apherese_CARneg_CD5pos.bwamem2.mutect2_VEP.filtered.vcf.gz",
  strelka.indels.filtered.bwamem2 = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/P2251_CARpos_CD5dim_vs_Apherese_CARneg_CD5pos.bwamem2.strelka.somatic_indels_VEP.filtered.vcf.gz",
  strelka.snvs.filtered.bwamem2 = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/P2251_CARpos_CD5dim_vs_Apherese_CARneg_CD5pos.bwamem2.strelka.somatic_snvs_VEP.filtered.vcf.gz"
)
dfs.somatic.P2251_CARpos_CD5dim = lapply(vcf.paths.somatic.P2251_CARpos_CD5dim, function(inpath)prep_df(inpath, gtcols = 9, drop_GT = T))
unfiltered.IDs[["P2251_CARpos_CD5dim"]] = lapply(dfs.somatic.P2251_CARpos_CD5dim, function(df)df$ID)
dfs.filtered.somatic.P2251_CARpos_CD5dim = lapply(dfs.somatic.P2251_CARpos_CD5dim, function(in_df){
  filter_vcf_inputs_by_relevance(in_df) %>%
    select(CHROM:FILTER) %>% 
    mutate(INFO = "DP=50") # dummy
})
filtered.IDs[["P2251_CARpos_CD5dim"]] = lapply(dfs.filtered.somatic.P2251_CARpos_CD5dim, function(df)df$ID)
somatic.P2251_CARpos_CD5dim.all = list(
  mutect.filtered = inner_join(dfs.filtered.somatic.P2251_CARpos_CD5dim$mutect.filtered.dragmap,
                               dfs.filtered.somatic.P2251_CARpos_CD5dim$mutect.filtered.bwamem2),
  strelka.indels.filtered = semi_join(dfs.filtered.somatic.P2251_CARpos_CD5dim$strelka.indels.filtered.dragmap,
                                       dfs.filtered.somatic.P2251_CARpos_CD5dim$strelka.indels.filtered.bwamem2, by = "ID"),
  strelka.snvs.filtered = semi_join(dfs.filtered.somatic.P2251_CARpos_CD5dim$strelka.snvs.filtered.dragmap,
                                     dfs.filtered.somatic.P2251_CARpos_CD5dim$strelka.snvs.filtered.bwamem2, by = "ID")
) %>% do.call(what=rbind) %>% distinct()

################################################################################

combined.somatic.variants = rbind(somatic.CARneg_CD5pos_merged.all,
                                  somatic.CARpos_CD5dim_merged.all,
                                  somatic.CARpos_CD5pos_merged.all) %>%
  arrange(CHROM.order(CHROM), POS) %>%
  distinct() %>%
  mutate(from_initial_round = T,
         initial_somatic = T)

# Apherese sample has different tools
vcf.paths.preexisting.Apherese_CARneg_CD5pos = list(
  deepvar.filtered.dragmap = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/Apherese_CARneg_CD5pos.dragmap.deepvariant_VEP.filtered.vcf.gz",
  strelka.filtered.dragmap = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/Apherese_CARneg_CD5pos.dragmap.strelka.variants_VEP.filtered.vcf.gz",
  deepvar.filtered.bwamem2 = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/Apherese_CARneg_CD5pos.bwamem2.deepvariant_VEP.filtered.vcf.gz",
  strelka.filtered.bwamem2 = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/Apherese_CARneg_CD5pos.bwamem2.strelka.variants_VEP.filtered.vcf.gz"
)
df.preexisting.Apherese_CARneg_CD5pos = lapply(vcf.paths.preexisting.Apherese_CARneg_CD5pos, function(inpath)prep_df(inpath, gtcols = 9, drop_GT = T))
unfiltered.IDs[["Apherese_CARneg_CD5pos"]] = lapply(df.preexisting.Apherese_CARneg_CD5pos, function(df)df$ID)
dfs.filtered.preexisting.Apherese_CARneg_CD5pos = lapply(df.preexisting.Apherese_CARneg_CD5pos, function(in_df){
  filter_vcf_inputs_by_relevance(in_df) %>%
    select(CHROM:FILTER) %>% 
    mutate(INFO = "DP=50") # dummy
})
filtered.IDs[["Apherese_CARneg_CD5pos"]] = lapply(dfs.filtered.preexisting.Apherese_CARneg_CD5pos, function(df)df$ID)
preexisting.Apherese_CARneg_CD5pos.all = list(
  deepvar.filtered = semi_join(dfs.filtered.preexisting.Apherese_CARneg_CD5pos$deepvar.filtered.dragmap,
                                dfs.filtered.preexisting.Apherese_CARneg_CD5pos$deepvar.filtered.bwamem2, by = "ID"),
  strelka.filtered = semi_join(dfs.filtered.preexisting.Apherese_CARneg_CD5pos$strelka.filtered.dragmap,
                                dfs.filtered.preexisting.Apherese_CARneg_CD5pos$strelka.filtered.bwamem2, by = "ID")
) %>% do.call(what=rbind) %>% 
  distinct() %>%
  mutate(from_initial_round = T,
         initial_somatic = F)

combined.variants = rbind(combined.somatic.variants,
                          preexisting.Apherese_CARneg_CD5pos.all) %>%
  distinct() %>%
  group_by(ID) %>%
  slice_head(n=1) %>% # take first in case of duplicates with different QUALs
  ungroup() %>%
  arrange(CHROM.order(CHROM), POS)

# --------------------------------------------------------------------------------

# Read forced call output file, parse out all GTs and VAFs for all samples,
# filter down to our candidate list, irrespective if FilterMutectCalls assigned PASS or not

final.vcf = read.vcfR(file = "../variant_calls_prefiltered/final/Koeln_1.mutect2.filtered_VEP.normalized.vcf.gz")  
final.obj = vcfR2tidy(final.vcf, single_frame = T)
final.df = final.obj$dat %>%
  mutate(ID = paste(CHROM, POS, paste0(REF, ">", ALT), sep = "_")) %>%
  semi_join(combined.variants, by = c("CHROM", "POS"))

final.df.GT = final.df %>% 
  select(ID, Indiv:gt_GT_alleles)
final.df = final.df %>%
  select(!c(Indiv:gt_GT_alleles)) %>%
  distinct()
final.df.csq = final.df %>%
  select(ID, CSQ) %>%
  dplyr::filter(!is.na(CSQ)) %>% 
  separate_longer_delim(cols = CSQ, delim = ",") %>%
  mutate(CSQ = stri_trim(CSQ, side = "both", pattern = "[^()]")) %>%
  separate_wider_delim(cols = CSQ, delim = "|", names = CSQ.names)
final.df = final.df %>%
  select(!CSQ)

# We need to recompute AF here to get the empirical VAFs we want.
# See: https://github.com/broadinstitute/gatk/issues/8080
final.df.VAF = final.df.GT %>%
  select(ID, Indiv, gt_AF, gt_DP, gt_AD) %>%
  separate_wider_delim(cols = gt_AD, delim = ",", names = c("REF_AD", "ALT_AD")) %>%
  mutate(bayesian_AF = as.double(gt_AF),
         REF_AD = as.integer(REF_AD),
         ALT_AD = as.integer(ALT_AD),
         gt_AF = ALT_AD/(REF_AD+ALT_AD)) %>%
  pivot_wider(names_from = Indiv, values_from = c(gt_AF, bayesian_AF, gt_DP, REF_AD, ALT_AD)) %>%
  rowwise() %>%
  mutate(diff_minmax_bayesian_AFs = max(pick(bayesian_AF_Koeln_1_Apherese_CARneg_CD5pos,
                                        bayesian_AF_Koeln_1_PB_CARneg_CD5pos_merged,
                                        bayesian_AF_Koeln_1_PB_CARpos_CD5dim_merged,
                                        bayesian_AF_Koeln_1_PB_CARpos_CD5pos_merged) - 
                               min(pick(bayesian_AF_Koeln_1_Apherese_CARneg_CD5pos,
                                        bayesian_AF_Koeln_1_PB_CARneg_CD5pos_merged,
                                        bayesian_AF_Koeln_1_PB_CARpos_CD5dim_merged,
                                        bayesian_AF_Koeln_1_PB_CARpos_CD5pos_merged)))) %>%
  ungroup()
  

rarity_threshold = 0.001
csq.clean = final.df.csq %>%
  select(!c(AFR_AF:gnomADg_SAS_AF, MANE_SELECT:CCDS, SWISSPROT:UNIPARC)) %>%
  mutate(  MAX_AF = as.double(MAX_AF),
           Rare = (is.na(MAX_AF) | (MAX_AF < rarity_threshold)),
           ImpactHighMed = (IMPACT=="HIGH" | IMPACT=="MODERATE")) %>%
  filter(  str_detect(SIFT, pattern = fixed("deleterious")) |  # predicted deleterious by SIFT
             str_detect(PolyPhen, pattern = fixed("damaging")) | # or polyphen
             ((SYMBOL %in% all_genes_of_interest) & (Rare | ImpactHighMed)) | # or possibly impacting genes of interest
             (str_detect(CLIN_SIG, fixed("pathogenic")) | str_detect(CLIN_SIG, fixed("risk_factor")) | str_detect(CLIN_SIG, fixed("drug_response"))) | # annotated possibly pathogenic or risk factor or drug response associated
             ( (IMPACT=="HIGH" | (IMPACT=="MODERATE" & Rare)) & str_detect(CLIN_SIG, fixed("benign"), negate =T)) # further high-impact or rare moderate-impact variants that are not already known to be benign (note: the line above rescues all lines with conflicting benign&pathogenic annotation)
  ) %>%
  mutate(Population_maximum_AF = MAX_AF, 
         Population_maximum_AF_which = MAX_AF_POPS) %>%
  relocate(Existing_variation, .after = Gene) %>%
  relocate(c(SIFT, PolyPhen, Population_maximum_AF, Population_maximum_AF_which, CLIN_SIG, MOTIF_NAME:TRANSCRIPTION_FACTORS), .after = HGVSp)
  
force.call.list = combined.variants %>%
  select(ID, from_initial_round, initial_somatic)

cleaned.df = final.df %>%
  mutate(in_Panel_of_Normals = PON,
         Short_tandem_repeat = STR,
         Repeat_Unit = RU) %>%
  select(CHROM:FILTER, in_Panel_of_Normals, Short_tandem_repeat, Repeat_Unit) %>%
  left_join(final.df.VAF %>% select(ID, contains("gt_AF"), diff_minmax_bayesian_AFs), by = "ID") %>%
  left_join(csq.clean, by = "ID") %>%
  left_join(force.call.list, by = "ID") %>%
  relocate(QUAL:Repeat_Unit, .after = CLIN_SIG) %>%
  relocate(Population_maximum_AF, Population_maximum_AF_which, .before = Allele) %>%
  mutate(Gene_of_Interest = SYMBOL %in% all_genes_of_interest) %>%
  relocate(c(Gene_of_Interest, from_initial_round, initial_somatic), .before = "Allele") %>%
  relocate(gt_AF_Koeln_1_PB_CARpos_CD5dim, gt_AF_Koeln_1_P2251_CARpos_CD5dim, .after = "PUBMED") %>%
  rename(VEP_Allele = Allele) %>%
  arrange(IMPACT, desc(diff_minmax_bayesian_AFs))

########## Save everything and write out #######################################

if (!(dir.exists("./Final_Outputs"))) {
  dir.create("Final_Outputs")
}

saveRDS(cleaned.df, file = "Final_Outputs/cleaned.df.rds")
saveRDS(final.obj, file = "Final_Outputs/final.df.unfiltered.rds")
saveRDS(list(combined.somatic.variants.df = combined.somatic.variants,
             combined.variants.df = combined.variants), file = "Final_Outputs/force.calling.variants.rds")
saveRDS(list(filtered.IDs = filtered.IDs,
             unfiltered.IDs = unfiltered.IDs), file = "Final_Outputs/filtering.IDS.rds")
write_excel_csv2(cleaned.df, file = "Final_Outputs/cleaned.variants.all.VAFs.csv")
