.cran_packages = c(
  "yaml", "ggplot2","reshape2", "dplyr", "naturalsort", "devtools", "scales",
  "stringr", "Seurat", "tibble", "tidyr", "HGNChelper", "forcats", "cowplot",
  "rlang", "remotes", "scGate", "patchwork", "openxlsx", "scCustomize", "ggpubr",
  "tidyverse", "anndata"
)
.bioc_packages = c("infercnv")

# Install CRAN packages (if not already installed)
.inst = .cran_packages %in% installed.packages()
if (any(!.inst)) {
  install.packages(.cran_packages[!.inst], repos = "http://cran.rstudio.com/")
}

# Install bioconductor packages (if not already installed)
.inst <- .bioc_packages %in% installed.packages()
if (any(!.inst)) {
  library(BiocManager)
  BiocManager::install(.bioc_packages[!.inst], ask = T)
}

list.of.packages = c(.cran_packages, .bioc_packages)

## Loading library
for (pack in list.of.packages) {
  suppressMessages(library(
    pack,
    quietly = TRUE,
    verbose = FALSE,
    character.only = TRUE
  ))
}

# source("code/helper/styles.R")
# source("code/helper/functions_plots.R")
# source("code/helper/functions.R")
#
# theme_set(mytheme(base_size = 8))
# base.size = 8

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Load Data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")

se.meta.t = readRDS(paste0(manifest$multi_omics$work, "integration/seurat_harmony_t.Rds"))
se.meta.t = se.meta.t[, !se.meta.t$celltype %in% "dpT"]
se.meta.t@meta.data = droplevels(se.meta.t@meta.data)

se.meta.t$clonePseudoID = gsub("Clone", "Cl", se.meta.t$clonePseudoID)
se.meta.t$TopClones = ifelse(
  se.meta.t$clonePseudoID %in% c("Cl_1", "Cl_2", "Cl_3"),
  se.meta.t$clonePseudoID,
  "Cl_Other"
)
se.meta.t$TopClones = ifelse(
  is.na(se.meta.t$CTstrict),
  NA, se.meta.t$TopClones
)

clone.col = unique(se.meta.t$TopClones)
clone.col = clone.col[!is.na(clone.col)]
clone.col = setNames(c("#555555", "#004488", "#DDAA33", "#BB5566"), clone.col)

se.meta.t@meta.data = se.meta.t@meta.data %>% mutate(
  TopClones_2 = case_when(
    TopClones == "Cl_1" | TopClones == "Cl_3" ~ "Cl_1_3",
    TRUE ~ TopClones
  )
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Read single cell data and annotations
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Reference Annotation
gene_pos = read.table("https://data.broadinstitute.org/Trinity/CTAT/cnv/gencode_v19_gen_pos.complete.txt",header=F,sep="\t")
gene_pos = cbind(gene_pos,name = sapply(gene_pos[,1],function(x){strsplit(x,"|",fixed=TRUE)[[1]][[1]]}))

if (!dir.exists("data/inferCNV/")) {
  dir.create("data/inferCNV/")
}

se.tumor.cc = subset(se.meta.t, subset = orig.ident != "AphNB")
se.control.cc = subset(se.meta.t, subset = orig.ident == "AphNB")

# Merge
se.comb.cc = merge(
  x=se.tumor.cc,
  y=se.control.cc,collapse = T
)

Type = rep("Normal",length(se.comb.cc$orig.ident))
Type[se.comb.cc$orig.ident != "AphNB"] = "Tumor"

se.comb.cc <- AddMetaData(
  object = se.comb.cc,
  metadata = Type,
  col.name = 'Type'
)


raw_count_matrix = as.matrix(se.comb.cc@assays$RNA@data)
annotations = data.frame(Type = se.comb.cc$Type)

mtc = match(rownames(raw_count_matrix), gene_pos$name)
raw_count_matrix = raw_count_matrix[!is.na(mtc),]
gene_pos = gene_pos[mtc[!is.na(mtc)],]
gene_order = gene_pos[,2:4]
rownames(gene_order) = gene_pos$name

# Buid object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=raw_count_matrix,
                                    annotations_file=annotations,
                                    gene_order_file=gene_order,
                                    ref_group_names=c("Normal"))

infercnv_obj = infercnv::run(infercnv_obj,cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=tempfile(),
                             cluster_by_groups=TRUE,
                             denoise=TRUE,
                             HMM=TRUE,
                             HMM_type = 'i3',
                             num_threads = 20,
                             plot_probabilities=F,
                             plot_steps=F,
                             no_prelim_plot=T
)

plot_cnv(infercnv_obj,output_filename = "data/inferCNV/inferCNV",title = "Tumor=All t-cells")
saveRDS(infercnv_obj,file="data/inferCNV/infercnv_obj_complete.rds")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Analyze subclones
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.control.cc = subset(se.meta.t, subset = orig.ident == "AphNB")

for (i in unique(se.meta.t$TopClones_2)){
  se.tumor.cc = subset(se.meta.t, subset = orig.ident != "AphNB" & se.meta.t$TopClones_2 == i)

  # Merge
  se.comb.cc = merge(
    x=se.tumor.cc,
    y=se.control.cc,collapse = T
  )

  Type = rep("Normal",length(se.comb.cc$orig.ident))
  Type[se.comb.cc$orig.ident != "AphNB"] = paste0("Tumor_",i)

  se.comb.cc <- AddMetaData(
    object = se.comb.cc,
    metadata = Type,
    col.name = 'Type'
  )


  raw_count_matrix = as.matrix(se.comb.cc@assays$RNA@data)
  annotations = data.frame(Type = se.comb.cc$Type)

  mtc = match(rownames(raw_count_matrix), gene_pos$name)
  raw_count_matrix = raw_count_matrix[!is.na(mtc),]
  gene_pos = gene_pos[mtc[!is.na(mtc)],]
  gene_order = gene_pos[,2:4]
  rownames(gene_order) = gene_pos$name

  # Buid object
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=raw_count_matrix,
                                      annotations_file=annotations,
                                      gene_order_file=gene_order,
                                      ref_group_names=c("Normal"))

  infercnv_obj = infercnv::run(infercnv_obj,cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir=tempfile(),
                               cluster_by_groups=TRUE,
                               denoise=TRUE,
                               HMM=TRUE,
                               HMM_type = 'i3',
                               num_threads = 20,
                               plot_probabilities=F,
                               plot_steps=F,
                               no_prelim_plot=T
  )

  plot_cnv(infercnv_obj,output_filename = paste0("data/inferCNV/inferCNV_",i),title = paste0("Tumor = non AphNB & ",i))
  saveRDS(infercnv_obj,file=paste0("data/inferCNV/infercnv_obj_",i,".rds"))

}
