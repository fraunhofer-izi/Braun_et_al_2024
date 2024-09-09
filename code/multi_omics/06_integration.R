# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Libs
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c(
  "yaml", "ggplot2", "reshape2", "dplyr", "foreach", "naturalsort", "ggthemes",
  "cowplot", "devtools", "scales", "stringr", "harmony", "MetBrewer",
  "Seurat", "future", "scCustomize", "scGate"
)
.bioc_packages = c("scRepertoire", "UCell")

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

if (any(!"SeuratWrappers" %in% installed.packages())) {
  Sys.unsetenv("GITHUB_PAT")
  remotes::install_github('satijalab/seurat-wrappers')
}
library(SeuratWrappers)

if (any(!"SignatuR" %in% installed.packages())) {
  Sys.unsetenv("GITHUB_PAT")
  remotes::install_github("carmonalab/SignatuR")
}
library(SignatuR)

if (any(!"STACAS" %in% installed.packages())) {
  Sys.unsetenv("GITHUB_PAT")
  remotes::install_github("carmonalab/STACAS")
}
library(STACAS)

if (any(!"ProjecTILs" %in% installed.packages())) {
  Sys.unsetenv("GITHUB_PAT")
  remotes::install_github("carmonalab/ProjecTILs")
}
library(ProjecTILs)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Functions
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
source("code/helper/styles.R")
source("code/helper/functions.R")
source("code/helper/functions_plots.R")
theme_set(mytheme())

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")

output.file = paste0(manifest$multi_omics$work, "integration/seurat_harmony.Rds")
output.file.t = paste0(manifest$multi_omics$work, "integration/seurat_harmony_t.Rds")
output.file.cd4 = paste0(manifest$multi_omics$work, "integration/seurat_harmony_cd4.Rds")
output.file.cd8 = paste0(manifest$multi_omics$work, "integration/seurat_harmony_cd8.Rds")

se.meta = readRDS(paste0(manifest$multi_omics$work, "03_seurat_anno_2.Rds"))
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Integration
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cell.idents = sort(table(se.meta$celltype_short_3))
print(cell.idents)
se.meta = se.meta[, se.meta$celltype_short_3 %in% names(cell.idents[cell.idents >= 50])]
se.meta@meta.data = droplevels(se.meta@meta.data)

dims.use.rna = 20

se.meta = integration(
  obj = se.meta,
  no.ftrs = 2000,
  threads = 20,
  .nbr.dims = dims.use.rna,
  run.integration = T,
  harmony.group.vars = c("orig.ident")
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Save
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
DefaultAssay(se.meta) = "RNA"
saveRDS(se.meta, output.file)
# se.meta = readRDS(output.file)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Add Clonotypes
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cellranger.dirs = list.dirs(
  path = manifest$multi_omics$cellranger, full.names = T, recursive = F
)
cellranger.samples = basename(cellranger.dirs)
fltrd.vdj = paste0(
  cellranger.dirs, "/outs/per_sample_outs/", cellranger.samples, "/vdj_t",
  "/filtered_contig_annotations.csv"
)
names(fltrd.vdj) = gsub("multi_", "", cellranger.samples)
print(paste0("# of files: ", length(fltrd.vdj)))

contig_list <- lapply(fltrd.vdj, function(x) {
  tryCatch(read.csv(x), error=function(e) NULL)
})
length(contig_list)
print(paste0("Empty file: ", names(lengths(contig_list)[lengths(contig_list) == 0])))

combined <- combineTCR(
  contig_list,
  samples = paste0(names(contig_list))
)

combined = data.table::rbindlist(combined) %>% data.frame()
combined$sample = se.meta$orig.ident[match(combined$barcode, rownames(se.meta@meta.data))]
combined = combined[!is.na(combined$sample), ]

max.clonotypes = max(table(combined$CTstrict))
se.meta <- combineExpression(
  combined, se.meta,
  cloneCall = "strict",
  group.by = NULL,
  # group.by = "sample",
  proportion = F,
  cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=max.clonotypes)
)
se.meta$cloneSize = droplevels(se.meta$cloneSize)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# T-cells
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.meta.t = se.meta[, grepl("^CD4|^CD8|^gd|^dp|Estimable", se.meta$celltype_short_3)]

se.meta.t =subset(se.meta.t, VDJ_T_AVAIL == TRUE)

se.meta.t@meta.data = droplevels(se.meta.t@meta.data)
se.meta.t@meta.data = se.meta.t@meta.data %>% mutate(
  celltype = case_when(
    grepl("^CD4", celltype) & CellCycle == T ~ "CD4.Cycling",
    grepl("^CD8", celltype) & CellCycle == T ~ "CD8.Cycling",
    TRUE ~ celltype
  )
)
se.meta.t$celltype = as.factor(se.meta.t$celltype)

cell.idents = sort(table(se.meta.t$celltype))
se.meta.t = se.meta.t[, se.meta.t$celltype %in% names(cell.idents[cell.idents > 10])]
se.meta@meta.data = droplevels(se.meta@meta.data)
se.meta.t@meta.data = droplevels(se.meta.t@meta.data)

se.meta.t = integration(
  obj = se.meta.t,
  no.ftrs = 1000,
  threads = 20,
  .nbr.dims = 15,
  min.dist = .4,
  harmony.group.vars = c("orig.ident")
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pd = se.meta.t@meta.data
pd = pd[!is.na(pd$CTstrict), ]
pd = pd %>%
  dplyr::group_by(CTstrict) %>%
  dplyr::count(name = "cloneFreqAll") %>%
  dplyr::arrange(desc(cloneFreqAll))
pd$pseudo_id = paste0("Clone_", seq(1:nrow(pd)))
se.meta.t$cloneFreqAll = pd$cloneFreqAll[match(se.meta.t$CTstrict, pd$CTstrict)]
se.meta.t$clonePseudoID = pd$pseudo_id[match(se.meta.t$CTstrict, pd$CTstrict)]

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Save
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
DefaultAssay(se.meta) = "RNA"
saveRDS(se.meta.t, output.file.t)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.cd8 = integration(
  obj = se.meta.t[, grepl("^CD8", se.meta.t$celltype_short_3)],
  no.ftrs = 1500,
  threads = 20,
  .nbr.dims = 15,
  harmony.group.vars = c("orig.ident")
)

se.cd4 = integration(
  obj = se.meta.t[, grepl("^CD4", se.meta.t$celltype_short_3)],
  no.ftrs = 1500,
  threads = 20,
  .nbr.dims = 15,
  harmony.group.vars = c("orig.ident")
)

saveRDS(se.cd8, output.file.cd8)
saveRDS(se.cd4, output.file.cd4)

