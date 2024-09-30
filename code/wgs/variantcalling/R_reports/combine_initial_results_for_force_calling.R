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

# filter functions
filter_vcf_inputs_somatic <- function(vcf_in,
                                      rarity_threshold=0.001){
  
  df = vcfR2tidy(vcf_in, single_frame = T)$dat %>%
    select(!(Indiv:gt_GT_alleles)) %>%
    distinct() %>%
    mutate(ID = paste(CHROM, POS, paste0(REF, ">", ALT), sep = "_"))
  
  csq = df %>%
    select(ID, CSQ) %>%
    dplyr::filter(!is.na(CSQ)) %>% 
    separate_longer_delim(cols = CSQ, delim = ",") %>%
    mutate(CSQ = stri_trim(CSQ, side = "both", pattern = "[^()]")) %>%
    separate_wider_delim(cols = CSQ, delim = "|", names = CSQ.names)
    
  df = df %>% select(!CSQ)
  
  IDs.missing.anno = setdiff(unique(df$ID), unique(csq$ID))
  if (length(IDs.missing.anno)>0) {
    warning(paste("Annotation seemingly missing in", length(IDs.missing.anno), "variants. Skipping them during filtering."))
  }
  
  csq.filtered = csq %>%
    mutate(  MAX_AF = as.double(MAX_AF),
             Rare = (is.na(MAX_AF) | (MAX_AF < rarity_threshold)),
             ImpactHighMed = (IMPACT=="HIGH" | IMPACT=="MODERATE")) %>%
    filter(  str_detect(SIFT, pattern = fixed("deleterious")) |  # predicted deleterious by SIFT
             str_detect(SIFT, pattern = fixed("damaging")) | # or polyphen
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

#########################

# extract all samples separately then combine

## CARneg_CD5pos_merged
vcf.paths.somatic.CARneg_CD5pos_merged = list(
  mutect.filtered = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARneg_CD5pos_merged_vs_Apherese_CARneg_CD5pos.dragmap.mutect2_VEP.filtered.vcf.gz",
  strelka.indels.filtered = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARneg_CD5pos_merged_vs_Apherese_CARneg_CD5pos.dragmap.strelka.somatic_indels_VEP.filtered.vcf.gz",
  strelka.snvs.filtered = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARneg_CD5pos_merged_vs_Apherese_CARneg_CD5pos.dragmap.strelka.somatic_snvs_VEP.filtered.vcf.gz"
)
df.filtered.somatic.CARneg_CD5pos_merged = lapply(vcf.paths.somatic.CARneg_CD5pos_merged, function(inpath){
    read.vcfR(inpath, limit = 1e10, cols = 9) %>% 
    filter_vcf_inputs_somatic() %>% # only format, skip individual genotype cols
    select(CHROM:FILTER) %>% 
    mutate(INFO = "DP=50") # dummy
  }) %>% do.call(what=rbind)

## CARpos_CD5dim_merged
vcf.paths.somatic.CARpos_CD5dim_merged = list(
  mutect.filtered = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARpos_CD5dim_merged_vs_Apherese_CARneg_CD5pos.dragmap.mutect2_VEP.filtered.vcf.gz",
  strelka.indels.filtered = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARpos_CD5dim_merged_vs_Apherese_CARneg_CD5pos.dragmap.strelka.somatic_indels_VEP.filtered.vcf.gz",
  strelka.snvs.filtered = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARpos_CD5dim_merged_vs_Apherese_CARneg_CD5pos.dragmap.strelka.somatic_snvs_VEP.filtered.vcf.gz"
)
df.filtered.somatic.CARpos_CD5dim_merged = lapply(vcf.paths.somatic.CARpos_CD5dim_merged, function(inpath){
  read.vcfR(inpath, limit = 1e10, cols = 9) %>% 
    filter_vcf_inputs_somatic() %>%
    select(CHROM:FILTER) %>% 
    mutate(INFO = "DP=50") # dummy
}) %>% do.call(what=rbind)

## CARpos_CD5pos_merged
vcf.paths.somatic.CARpos_CD5pos_merged = list(
  mutect.filtered = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARpos_CD5pos_merged_vs_Apherese_CARneg_CD5pos.dragmap.mutect2_VEP.filtered.vcf.gz",
  strelka.indels.filtered = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARpos_CD5pos_merged_vs_Apherese_CARneg_CD5pos.dragmap.strelka.somatic_indels_VEP.filtered.vcf.gz",
  strelka.snvs.filtered = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARpos_CD5pos_merged_vs_Apherese_CARneg_CD5pos.dragmap.strelka.somatic_snvs_VEP.filtered.vcf.gz"
)
df.filtered.somatic.CARpos_CD5pos_merged = lapply(vcf.paths.somatic.CARpos_CD5pos_merged, function(inpath){
  read.vcfR(inpath, limit = 1e10, cols = 9) %>% 
    filter_vcf_inputs_somatic() %>%
    select(CHROM:FILTER) %>% 
    mutate(INFO = "DP=50") # dummy
}) %>% do.call(what=rbind)

## CARpos_CD5dim (Run 1)
vcf.paths.somatic.CARpos_CD5dim = list(
  mutect.filtered = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARpos_CD5dim_vs_Apherese_CARneg_CD5pos.dragmap.mutect2_VEP.filtered.vcf.gz",
  strelka.indels.filtered = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARpos_CD5dim_vs_Apherese_CARneg_CD5pos.dragmap.strelka.somatic_indels_VEP.filtered.vcf.gz",
  strelka.snvs.filtered = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/PB_CARpos_CD5dim_vs_Apherese_CARneg_CD5pos.dragmap.strelka.somatic_snvs_VEP.filtered.vcf.gz"
)
df.filtered.somatic.CARpos_CD5dim = lapply(vcf.paths.somatic.CARpos_CD5dim, function(inpath){
  read.vcfR(inpath, limit = 1e10, cols = 9) %>% 
    filter_vcf_inputs_somatic() %>%
    select(CHROM:FILTER) %>% 
    mutate(INFO = "DP=50") # dummy
}) %>% do.call(what=rbind)

## P2251_CARpos_CD5dim (Run 2)
vcf.paths.somatic.P2251_CARpos_CD5dim = list(
  mutect.filtered = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/P2251_CARpos_CD5dim_vs_Apherese_CARneg_CD5pos.dragmap.mutect2_VEP.filtered.vcf.gz",
  strelka.indels.filtered = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/P2251_CARpos_CD5dim_vs_Apherese_CARneg_CD5pos.dragmap.strelka.somatic_indels_VEP.filtered.vcf.gz",
  strelka.snvs.filtered = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/P2251_CARpos_CD5dim_vs_Apherese_CARneg_CD5pos.dragmap.strelka.somatic_snvs_VEP.filtered.vcf.gz"
)
df.filtered.somatic.P2251_CARpos_CD5dim = lapply(vcf.paths.somatic.P2251_CARpos_CD5dim, function(inpath){
  read.vcfR(inpath, limit = 1e10, cols = 9) %>% 
    filter_vcf_inputs_somatic() %>%
    select(CHROM:FILTER) %>% 
    mutate(INFO = "DP=50") # dummy
}) %>% do.call(what=rbind)

combined.somatic.variants = rbind(df.filtered.somatic.CARpos_CD5dim_merged,
                                  df.filtered.somatic.CARpos_CD5pos_merged,
                                  df.filtered.somatic.CARpos_CD5pos_merged) %>%
  arrange(CHROM.order(CHROM), POS) %>%
  distinct()

# Apherese sample has different tools
vcf.paths.preexisting.Apherese_CARneg_CD5pos = list(
  deepvar.filtered = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/Apherese_CARneg_CD5pos.dragmap.deepvariant_VEP.filtered.vcf.gz",
  strelka.filtered = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/Apherese_CARneg_CD5pos.dragmap.strelka.variants_VEP.filtered.vcf.gz"
)
df.filtered.preexisting.Apherese_CARneg_CD5pos = lapply(vcf.paths.preexisting.Apherese_CARneg_CD5pos, function(inpath){
  read.vcfR(inpath, limit = 1e10, cols = 9) %>% # only format 
    filter_vcf_inputs_somatic() %>%
    select(CHROM:FILTER) %>% 
    mutate(INFO = "DP=50") # dummy
}) %>% do.call(what=rbind)

combined.variants = rbind(combined.somatic.variants,
                          df.filtered.preexisting.Apherese_CARneg_CD5pos) %>%
  arrange(CHROM.order(CHROM), POS) %>%
  distinct()

combined.variants = combined.variants %>%
  group_by(ID) %>%
  slice_head(n=1) %>% # take first in case of duplicats with different QUALs
  ungroup() %>%
  arrange(CHROM.order(CHROM), POS)

combined.variants.mt = combined.variants %>%
  mutate(ID = NA, QUAL = NA, POS = paste(POS)) %>%
  as.matrix()

# load one of the files for the header, the rest is practically irrelevant
variants.vcf = read.vcfR(file = "../variant_calls_prefiltered/sarek_GATK_hg38_plus_Carvykti/P2251_CARpos_CD5dim_vs_Apherese_CARneg_CD5pos.dragmap.strelka.somatic_snvs_VEP.filtered.vcf.gz", cols = NA)
# replace content with combined variants and write back
variants.vcf@fix = combined.variants.mt
write.vcf(variants.vcf, file = "../merged_forced_alleles.vcf.gz")
