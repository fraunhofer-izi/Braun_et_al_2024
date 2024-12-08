library(conflicted)
library(tidyverse)
conflict_prefer("filter", winner = "dplyr")
library(stringr)
library(stringi)
library(readr)
library(vcfR)
library(ggplot2)

########## Addendum: all mutations for mutational-burden analysis ##############

final.vcf = read.vcfR(file = "$OUTDIR/variant_calls_prefiltered/final/Koeln_1.mutect2.filtered_VEP.normalized.vcf.gz")
final.obj = vcfR2tidy(final.vcf, single_frame = T)

tmp = final.obj$dat %>%
  mutate(ID = paste(CHROM, POS, paste0(REF, ">", ALT), sep = "_")) %>%
  filter(FILTER == "PASS") %>%
  filter(gt_AF >= 0.05)

tmp.csq = tmp %>%
  select(ID, CSQ) %>%
  dplyr::filter(!is.na(CSQ)) %>%
  separate_longer_delim(cols = CSQ, delim = ",") %>%
  mutate(CSQ = stri_trim(CSQ, side = "both", pattern = "[^()]")) %>%
  separate_wider_delim(cols = CSQ, delim = "|", names = CSQ.names) %>%
  distinct()
tmp = tmp %>%
  select(!CSQ)

rarity_threshold = 0.001
tmp.csq.clean = tmp.csq %>%
  select(!c(AFR_AF:gnomADg_SAS_AF, MANE_SELECT:CCDS, SWISSPROT:UNIPARC)) %>%
  mutate(  MAX_AF = as.double(MAX_AF),
           Rare = (is.na(MAX_AF) | (MAX_AF < rarity_threshold))) %>%
  mutate(Population_maximum_AF = MAX_AF,
         Population_maximum_AF_which = MAX_AF_POPS) %>%
  relocate(Existing_variation, .after = Gene) %>%
  relocate(c(SIFT, PolyPhen, Population_maximum_AF, Population_maximum_AF_which, CLIN_SIG, MOTIF_NAME:TRANSCRIPTION_FACTORS), .after = HGVSp)

tmp.csq.typing = tmp.csq.clean %>%
  filter(Rare) %>%
  group_by(ID) %>%
  summarise(EXON = first(na.omit(na_if(EXON, ""))),
            INTRON = first(na.omit(na_if(INTRON, ""))),
            HGVSp = first(na.omit(na_if(HGVSp, ""))))

tmp.csq.typing = tmp.csq.typing %>%
  mutate(synonymous = str_detect(HGVSp, pattern = fixed("%3D")),
         type = ifelse(is.na(EXON), ifelse(is.na(INTRON), "intergenic", "intronic"), ifelse(synonymous, "exonic_synonymous", "exonic_nonsynonymous")))

tmp.plt = tmp %>%
  filter(Indiv %in% c("Koeln_1_PB_CARpos_CD5pos_merged", "Koeln_1_PB_CARpos_CD5dim_merged", "Koeln_1_PB_CARneg_CD5pos_merged")) %>%
  separate_wider_delim(cols = gt_AD, delim = ",", names = c("REF_AD", "ALT_AD")) %>%
  mutate(bayesian_AF = as.double(gt_AF),
         REF_AD = as.integer(REF_AD),
         ALT_AD = as.integer(ALT_AD),
         gt_AF = ALT_AD/(REF_AD+ALT_AD)) %>%
  filter(ALT_AD>=2 & gt_AF>=0.1) %>%
  inner_join(tmp.csq.typing %>% select(ID, type), by = "ID", relationship = "many-to-one") %>%
  transmute(ID = ID,
            Indiv = Indiv,
            type = type)

############# Plot WGS TMB for all post CAR-T infusion samples #################

hg38_eff_genome_size = 3049 # (in MB)

tmp.plt = tmp.plt %>%
  filter(type!="exonic_synonymous") %>%
  mutate(
    Indiv = factor(Indiv, levels = c("Koeln_1_PB_CARpos_CD5dim_merged", "Koeln_1_PB_CARpos_CD5pos_merged", "Koeln_1_PB_CARneg_CD5pos_merged"))
  )

wgs_freq_df = tmp.plt %>%
  group_by(Indiv) %>%
  summarise( wgs_freq = n()/hg38_eff_genome_size)

tmp.plt %>%
  filter(type %in% c("exonic_nonsynonymous", "exonic_non-coding-region")) %>%
  nrow()

standard_chromosomes = paste0("chr", c(1:22, "X", "Y"))
cds = rtracklayer::import("~/GENOMES/human/gencode_hg38/gencode.v46.primary_assembly.annotation.gtf.gz", format="gtf")
cds_rows = subset(cds, transcript_type == "protein_coding" & type == "exon")
cds_rows = cds_rows[seqnames(cds_rows) %in% standard_chromosomes]
cds_rows = reduce(cds_rows)
sum(cds_rows@ranges@width)

wes_freq_df = tmp.plt %>%
  filter(type %in% c("exonic_nonsynonymous", "exonic_non-coding-region")) %>%
  group_by(Indiv) %>%
  summarise( wes_freq = n()/ (sum(cds_rows@ranges@width)/1000000) )

tmp.plt = tmp.plt %>%
  left_join(wgs_freq_df) %>%
  left_join(wes_freq_df)

# eff genome size in Schrader et al was 64 Mb
range_min = 13/64*hg38_eff_genome_size # lowest Schrader et al
range_max = 150/64*hg38_eff_genome_size # highest Schrader et al

tmp.plt %>%
  ggplot(aes(x=Indiv, fill=type)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = range_min, ymax = range_max),alpha=0.3, fill="grey85") +
  geom_hline(aes(yintercept=1.1*hg38_eff_genome_size), linetype="dashed", color="grey10")+
  geom_bar(position = position_stack(), color="black") +
  scale_y_continuous(limits = c(0,7200), sec.axis = sec_axis(~ . / hg38_eff_genome_size, name = "Tumor mutational burden (Mut/Mb)")) +
  ylab("Somatic mutation frequency") +
  scale_x_discrete(labels = c("PB CAR+CD5-", "PB CAR+CD5+", "PB CAR-CD5+")) +
  scale_fill_manual(name="", values = rev(c("#004488", "#ABD2DFFF" , "#AFD9A5FF", "#DDAA33FF")),labels = c("UTR", "CDS", "Intergenic", "Intronic") ) +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1), axis.title.x = element_blank(), legend.position = "right") +
  annotate("text", x=2.5, y=6500, label="TMB range T-PLL\n(Schrader et al.)", color="grey10") +
  guides(fill = guide_legend(title.position = "top"))

ggsave2(
  "figures/tmb_plot.png",
  height=80, width=90, unit="mm", bg="white", dpi=400, scale=1.4
)
