# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Libraries
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c("Seurat", "yaml", "dplyr", "doParallel", "parallel", "data.table", "Matrix")
.bioc_packages = c("SingleCellExperiment", "scDblFinder")

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

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Load Rawcounts and create a merged Seurat object
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")
work.path = manifest$multi_omics$cellranger
seurat.path = manifest$multi_omics$work

obj.path = paste0(seurat.path, "01_seurat_merged.Rds")

cellranger.dirs = list.dirs(path = work.path, full.names = T, recursive = F)
cellranger.samples = basename(cellranger.dirs)
fltrd.dirs = paste0(
  cellranger.dirs, "/outs/per_sample_outs/", cellranger.samples,
  "/count/sample_filtered_feature_bc_matrix/"
)
names(fltrd.dirs) = cellranger.samples
print(length(fltrd.dirs))

# i = names(fltrd.dirs)[1]
bpparam = BiocParallel::MulticoreParam(workers = 4)
seurat.l = BiocParallel::bplapply(names(fltrd.dirs), function(i) {

  id = i
  fltrd.counts = Read10X(data.dir = fltrd.dirs[names(fltrd.dirs) == id], gene.column = 2)

  seu.obj = CreateSeuratObject(counts = fltrd.counts[[1]], project = id)
  seu.obj[["ADT"]] = CreateAssayObject(counts = fltrd.counts[[2]])

  seu.obj = RenameCells(seu.obj, new.names = gsub("multi_", "", colnames(seu.obj)))
  seu.obj@meta.data$orig.ident = gsub("multi_", "", id)

  sce = scDblFinder(GetAssayData(seu.obj, slot="counts"))
  df = data.frame(sce@colData) %>% dplyr::select(scDblFinder.score, scDblFinder.class)
  colnames(df) = c("scDblFinder_score", "scDblFinder_class")
  seu.obj = AddMetaData(seu.obj, df)

  seu.obj

}, BPPARAM = bpparam)

se.meta = merge(
  seurat.l[[1]], y = seurat.l[2:length(seurat.l)],
  add.cell.ids = names(seurat.l), project = "car_koeln"
)
se.meta@meta.data$orig.ident = factor(se.meta@meta.data$orig.ident)
se.meta[["RNA"]] <- JoinLayers(se.meta[["RNA"]])

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Assay ADT
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
DefaultAssay(se.meta) = "ADT"
adt.ftrs = gsub("totalseqC-", "", rownames(se.meta))

rownames(se.meta@assays$ADT@counts) = toupper(adt.ftrs)
se.meta@assays$ADT@data = se.meta@assays$ADT@counts
rownames(se.meta@assays$ADT@meta.features) = rownames(se.meta@assays$ADT@counts)
DefaultAssay(se.meta) = "RNA"
Idents(se.meta) = "orig.ident"

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Save
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
saveRDS(se.meta, file = obj.path)
