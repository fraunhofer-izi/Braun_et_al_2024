.cran_packages = c(
  "yaml", "ggplot2","reshape2", "dplyr", "naturalsort", "devtools", "scales",
  "stringr", "Seurat", "tibble", "tidyr", "HGNChelper", "forcats", "cowplot",
  "rlang", "remotes", "scGate", "patchwork", "openxlsx", "scCustomize", "ggpubr",
  "tidyverse", "scCustomize", "ggh4x", "ggrepel", "anndata", "scico"
)
.bioc_packages = c("dittoSeq", "scRepertoire","infercnv")

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


if (any(!"ProjecTILs" %in% installed.packages())) {
  remotes::install_github("carmonalab/ProjecTILs")
}
library(ProjecTILs)

# Source: https://www.genenames.org/data/genegroup/#!/group/370
tcell_receptor_genes = read.table("./data/T-cell_receptor_genes.txt",header=T,sep="\t")

# Add column in metadata wether CD4/CD8, CAR are present
cd4cd8_car = function(
    obj
){
  cd8cd4 = FetchData(obj, c("CD8A", "CD8B", "CD4", "CD274"), slot = "counts")
  ct.cd8cd4 = cd8cd4 %>% mutate(
    CD4CD8_BY_EXPRS = case_when(
      CD4 > 0 & CD8A == 0 & CD8B == 0 ~ "CD4+CD8-",
      CD4 == 0 & (CD8A > 0 | CD8B > 0) ~ "CD4-CD8+",
      CD4 == 0 & CD8A == 0 & CD8B == 0 ~ "CD4-CD8-",
      CD4 > 0 & (CD8A > 0 | CD8B > 0) ~ "CD4+CD8+",
      TRUE ~ "unresolved"
    )) %>%
    dplyr::select(CD4CD8_BY_EXPRS)
  rownames(ct.cd8cd4) = rownames(cd8cd4)
  obj = AddMetaData(obj, ct.cd8cd4)

  cd3 = FetchData(obj, c("CD3D", "CD3E", "CD3G"), slot = "counts")
  cd3 = cd3 %>% mutate(
    CD3_BY_EXPRS = case_when(
      CD3D > 0 | CD3E > 0 | CD3G > 0 ~ "CD3",
      TRUE ~ "unresolved"
    )) %>%
    dplyr::select(CD3_BY_EXPRS)
  obj = AddMetaData(obj, cd3)

  if("ciltacel" %in% rownames(GetAssayData(obj, slot = c("counts"), assay = "Spatial"))){
    car.ftr = FetchData(obj, c("ciltacel"), slot = "counts")
    obj$CAR_BY_EXPRS = as.factor(car.ftr$ciltacel > 0)
  } else {
    obj$CAR_BY_EXPRS = FALSE
    obj$CAR_BY_EXPRS = as.factor(obj$CAR_BY_EXPRS)
  }

  obj
}


source("code/helper/styles.R")
source("code/helper/functions_plots.R")
source("code/helper/functions.R")

theme_set(mytheme(base_size = 8))
base.size = 8

theme.custom = mytheme(base_size = 8) &
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.title = element_blank(),
    panel.border = element_blank(),
    plot.title = element_text(hjust = .5, size = rel(1), face = "bold")
  )

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
gc.anno = rtracklayer::import("./data/genes.gtf")
gene_pos = read.table("https://data.broadinstitute.org/Trinity/CTAT/cnv/gencode_v19_gen_pos.complete.txt",header=F,sep="\t")
gene_pos = cbind(gene_pos,name = sapply(gene_pos[,1],function(x){strsplit(x,"|",fixed=TRUE)[[1]][[1]]}))


if (!dir.exists("./analysis/multi_omics/inferCNV/"))
  dir.create("./analysis/multi_omics/inferCNV/")


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

plot_cnv(infercnv_obj,output_filename = "./analysis/multi_omics/inferCNV/inferCNV",title = "Tumor=All t-cells")
saveRDS(infercnv_obj,file="./analysis/multi_omics/inferCNV/infercnv_obj_complete.rds")







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

  plot_cnv(infercnv_obj,output_filename = paste0("./analysis/multi_omics/inferCNV/inferCNV_",i),title = paste0("Tumor = non AphNB & ",i))
  saveRDS(infercnv_obj,file=paste0("./analysis/multi_omics/inferCNV/infercnv_obj_",i,".rds"))

}
