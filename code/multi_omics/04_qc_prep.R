# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Libraries and some Functions
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c(
  "Seurat", "yaml", "dplyr", "stringr", "naturalsort", "data.table", "ggplot2",
  "scales", "openxlsx", "cowplot", "scCustomize"
)
.bioc_packages = c("dittoSeq", "clustifyr", "scds", "scDblFinder", "UCell")

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

if (any(!"Azimuth" %in% installed.packages())) {
  remotes::install_github('satijalab/azimuth', ref = 'master')
}
library(Azimuth)

if (any(!"Azimuth" %in% installed.packages())) {
  devtools::install_github('satijalab/seurat-data')
}
library(SeuratData)

if (any(!"SeuratDisk" %in% installed.packages())) {
  remotes::install_github("mojaveazure/seurat-disk")
}
library(SeuratDisk)

if (any(!"ProjecTILs" %in% installed.packages())) {
  Sys.unsetenv("GITHUB_PAT")
  remotes::install_github("carmonalab/ProjecTILs")
}
library(ProjecTILs)

source("code/helper/functions.R")
source("code/helper/functions_plots.R")
source("code/helper/styles.R")
theme_set(mytheme(base_size = 12))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Load objects and phenodata
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")

se.meta = readRDS(paste0(manifest$multi_omics$work, "01_seurat_merged.Rds"))
output.file = paste0(manifest$multi_omics$work, "02_seurat_pre.Rds")

se.meta@meta.data = se.meta@meta.data %>% dplyr::mutate(
  SOURCE = dplyr::case_when(
    orig.ident == "P2248" ~ "BM",
    TRUE ~ "PB"
  )
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Add column in metadata wether CD4/CD8, CAR are present
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.meta = cd4cd8_car_present(obj = se.meta)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Normalize
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.meta = NormalizeData(se.meta, assay = 'RNA', normalization.method = "LogNormalize")
se.meta = NormalizeData(se.meta, assay = "ADT", normalization.method = 'CLR', margin = 2)
slot(object = se.meta[["ADT"]], name = 'data') = as(se.meta@assays$ADT@data, "dgCMatrix")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Add %MT, %Ribosomal and complexity values
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.meta[["Perc_of_mito_genes"]] = Seurat::PercentageFeatureSet(se.meta, pattern = "^MT-")
se.meta[["Perc_of_ribosomal_genes"]] = Seurat::PercentageFeatureSet(se.meta, pattern = "^RPL|^RPS")
se.meta@meta.data$log10GenesPerUMI = log10(se.meta$nFeature_RNA) / log10(se.meta$nCount_RNA)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Cell filtering
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
nFeature_low_cutoff = 500
nFeature_high_cutoff = 6000
nCount_low_cutoff = 1000
nCount_high_cutoff = 50000
mt_cutoff = 15
complx_cutoff = 0.8

p1 = qc_vln_plot_cell(
    se.meta, low_cutoff = nFeature_low_cutoff, high_cutoff = nFeature_high_cutoff,
  )
p2 = qc_vln_plot_cell(
    se.meta, .features = "nCount_RNA", plot_title = "UMIs per Cell",
    y_axis_label = "UMIs", low_cutoff = nCount_low_cutoff, high_cutoff = nCount_high_cutoff,
  )
p3 = qc_vln_plot_cell(
  se.meta, .features = "Perc_of_mito_genes", plot_title = "Mito Gene % per Cell",
  y_axis_label = "% Mito Gene Counts", high_cutoff = mt_cutoff,
)
p4 = qc_vln_plot_cell(
  se.meta,  .features = "log10GenesPerUMI", plot_title = "Cell Complexity",
  y_axis_label = "log10(Genes) / log10(UMIs)", high_cutoff = complx_cutoff
)

ggsave2(
  filename="analysis/multi_omics_mir/figures/qc_pre_processing/stats_tech_per_cell.pdf",
  plot = cowplot::plot_grid(
    p1, p2, p3, p4, scale = .9, nrow = 1,
    labels = "AUTO", label_fontface = "bold", label_size = 14
  ),
  width = 160,
  height = 42,
  dpi = 100,
  bg = "white",
  units = "mm",
  scale = 2
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Cell filtering
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
label_cells_rm = function(obj) {
  obj@meta.data = obj@meta.data %>% mutate(
    KEEP_CELL = case_when(
      (nFeature_RNA < nFeature_low_cutoff) | (nFeature_RNA > nFeature_high_cutoff) |
      (nCount_RNA < nCount_low_cutoff) | (nCount_RNA > nCount_high_cutoff) |
      (Perc_of_mito_genes > mt_cutoff)  | (log10GenesPerUMI < complx_cutoff) ~ FALSE,
      TRUE ~ TRUE
    )
  )
  # print(table(obj$orig.ident, obj$KEEP_CELL))
  obj
}

cell.track = count_cells_per_sample(c(se.meta))
se.meta = label_cells_rm(se.meta)
se.meta = subset(se.meta, subset = KEEP_CELL == TRUE)
se.meta@meta.data = droplevels(se.meta@meta.data)
se.meta$KEEP_CELL = NULL
cell.track = count_cells_per_sample(c(se.meta), cell.track, "n1")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("scds")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
scds_doublets = function(se){
  suppressWarnings({
    suppressMessages({
      sce = as.SingleCellExperiment(se)
      sce = scds::cxds(sce)
      sce = scds::bcds(sce)
      sce = scds::cxds_bcds_hybrid(sce, estNdbl = T)
      CD  = data.frame(sce@colData)
      CD = CD %>% dplyr::select(
        cxds_score, bcds_score, hybrid_score, cxds_call, bcds_call, hybrid_call
      )
      CD$barcode = rownames(CD)
      CD
    })
  })
}

obj.l = Split_Object(se.meta, split.by = "orig.ident", threads = 4)
scds_res = parallel::mclapply(obj.l, function(se){
  scds_doublets(se)
}, mc.cores = 4)

scds_res = do.call("rbind", scds_res)
rownames(scds_res) = scds_res$barcode
scds_res$barcode = NULL
se.meta = AddMetaData(se.meta, scds_res)
se.meta$DOUBLETS_CONSENSUS = se.meta$scDblFinder_class == "doublet" & se.meta$hybrid_call == TRUE

cell.track = count_cells_per_sample(c(subset(se.meta, DOUBLETS_CONSENSUS == FALSE)), cell.track, "n2")

saveRDS(
  cell.track, file = "data/stats_celltrack.Rds"
)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# CellCycleScoring
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
source("data/metadata/signatures/cellCycleMarkers.R")
se.meta = CellCycleScoring(
  se.meta, s.features = s.genes, g2m.features = g2m.genes,
  assay = 'RNA', search = TRUE
)

# se = obj.l$MXMERZ002A_03
estimate_cc = function(se){

  suppressWarnings({
    suppressMessages({

      data(cell.cycle.obj) # ProjecTILs package
      DefaultAssay(se) = "RNA"

      se = se %>%
        FindVariableFeatures(verbose = F) %>%
        ScaleData(verbose = F) %>%
        RunPCA(npcs = 20, verbose = F) %>%
        FindNeighbors(reduction = "pca", dims = 1:20, verbose = F)

      tmp =  tryCatch(
        FindClusters(se, resolution = 2, verbose = F),
        error=function(e) "error"
      )
      if(class(tmp) != "Seurat") {
        se = FindClusters(se, resolution = 1, verbose = F)
      } else {
        se = tmp
      }

      # if(ncol(se) < 31) {
      #   se = RunUMAP(se, reduction = "pca", dims = 1:20, seed.use = 1234, verbose = F)
      # } else {
      #   se = RunUMAP(se, reduction = "pca", dims = 1:20, seed.use = 1234, verbose = F)
      # }

      quiet <- function(x) {
        sink(tempfile())
        on.exit(sink())
        invisible(force(x))
      }

      cc.phase = quiet(
        clustifyr::run_gsea(
          GetAssayData(object = se, assay = "RNA", slot = "data"),
          query_genes = cell.cycle.obj$human$cycling,
          cluster_ids =  se@meta.data[["seurat_clusters"]], n_perm = 1000
        )
      )

      cc.phase$pval_adj = p.adjust(cc.phase$pval, method = "BH")
      cc.cl = rownames(cc.phase[cc.phase$pval < 0.05, ])

      cc.cells = (se@meta.data[["seurat_clusters"]] %in% cc.cl)
      se@meta.data$CellCycle = factor(cc.cells)

      pd.cc = se@meta.data[cc.cells, ]
      pd.cc = pd.cc %>% dplyr::mutate(
        CellCycle_Phase = dplyr::case_when(
          G2M.Score > S.Score ~ "G2M",
          TRUE ~ "S"
        )
      )
      se$CellCycle_Phase = pd.cc$CellCycle_Phase[match(rownames(se@meta.data), rownames(pd.cc))]
      se$CellCycle_Phase[is.na(se$CellCycle_Phase)] = "G1M"
      se$CellCycle_Phase = factor(se$CellCycle_Phase, levels = c("G1M", "S", "G2M"))

      # (DimPlot_scCustom(se, reduction = "umap", group.by = "seurat_clusters", pt.size = 1) & mytheme() & theme(legend.position = "none") |
      #   DimPlot_scCustom(se, reduction = "umap", group.by = "CellCycle", pt.size = 1) & mytheme() |
      #   DimPlot_scCustom(se, reduction = "umap", group.by = "CellCycle_Phase", pt.size = 1) & mytheme() ) /
      #   (
      #   FeaturePlot_scCustom(
      #     se,
      #     reduction = "umap",
      #     features = c("S.Score", "G2M.Score"),
      #     pt.size = .5, na_cutoff = 0.1,
      #     colors_use = rev(MetBrewer::met.brewer("Hokusai1",n=100))
      #   ) & mytheme())

      res = se@meta.data %>% dplyr::select(CellCycle, CellCycle_Phase)
      res$barcode = rownames(res)
      res
    })
  })
}

obj.l = Split_Object(se.meta, split.by = "orig.ident", threads = 4)
cc.res = parallel::mclapply(obj.l, function(se){
  estimate_cc(se)
}, mc.cores = 4)
cc.res = do.call("rbind", cc.res)
rownames(cc.res) = cc.res$barcode
se.meta = AddMetaData(se.meta, cc.res)
se.meta$barcode = NULL

se.meta@meta.data = se.meta@meta.data %>%
  dplyr::mutate_if(is.character, as.factor)

# Annotate cells that were missed by the cluster approach
cutoff = 0.20
se.meta$CellCycle_Phase = dplyr::case_when(
  se.meta$S.Score > cutoff & se.meta$G2M.Score <= cutoff ~ "S",
  se.meta$G2M.Score > cutoff & se.meta$S.Score <= cutoff ~ "G2M",
  (se.meta$S.Score > cutoff & se.meta$G2M.Score > cutoff) & (se.meta$S.Score > se.meta$G2M.Score) ~ "S",
  (se.meta$S.Score > cutoff & se.meta$G2M.Score > cutoff) & (se.meta$G2M.Score > se.meta$S.Score) ~ "G2M",
  TRUE ~ se.meta$CellCycle_Phase
)
se.meta$CellCycle_Phase = factor(se.meta$CellCycle_Phase, levels = c("G1M", "S", "G2M"))
se.meta@meta.data$CellCycle = factor(ifelse(se.meta$CellCycle_Phase == "G1M", "FALSE", "TRUE"))

data(cell.cycle.obj) # ProjecTILs package
se.meta = AddModuleScore_UCell(
  se.meta,
  features = list("CellCycle_SCORE" = cell.cycle.obj$human$cycling),
  assay = "RNA", slot = "counts",
  ncores = 25, force.gc = T
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("Save")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
DefaultAssay(se.meta) = "RNA"
Idents(se.meta) = "orig.ident"

saveRDS(se.meta, output.file)
