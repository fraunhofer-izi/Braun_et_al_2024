# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Normalize, Harmony, Clustering
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
integration = function(
    obj,
    obj.l = NULL,
    no.ftrs = 500,
    .assay = "RNA",
    max.cl = 1,
    min.cells.per.sample = 25,
    threads = 5,
    .nbr.dims = 15,
    do.cluster = T,
    do.dimreduc = T,
    cc.regr = F,
    hvg.union = T,
    custom.features = NULL,
    run.integration = T,
    harmony.group.vars = NULL,
    obj.split.by = "orig.ident",
    .perp = 50,
    dmreduc.dims = NULL,
    n.neighbors = 30,
    min.dist = 0.3) {

  library(SignatuR)
  library(parallel)
  library(BiocParallel)
  library(harmony)

  if (any(!"readgmt" %in% installed.packages())) {
    Sys.unsetenv("GITHUB_PAT")
    devtools::install_github("jhrcook/readgmt")
  }
  library(readgmt)


  start.time <- Sys.time()

  if(is.null(dmreduc.dims)) {
    dmreduc.dims = .nbr.dims
  }

  DefaultAssay(obj) = "RNA"
  obj@meta.data = droplevels(obj@meta.data)
  obj = DietSeurat(obj, counts = TRUE, data = TRUE)
  obj = NormalizeData(obj)

  # Gene categories to exclude from variable genes
  bl <- c(
    SignatuR::GetSignature(SignatuR$Hs$Compartments$Mito)[[1]],
    SignatuR::GetSignature(SignatuR$Hs$Compartments$Immunoglobulins)[[1]],
    SignatuR::GetSignature(SignatuR$Hs$Compartments$TCR)[[1]]
    # SignatuR::GetSignature(SignatuR$Hs$Blocklists)[[1]]
  )
  bl = c(bl, c("RPS4Y1", "EIF1AY", "DDX3Y", "KDM5D", "XIST")) # gender genes
  bl <- unique(bl)

  if (hvg.union == T) {

    if(is.null(obj.l)) {
      print("Split object")
      obj.l = Split_Object(obj, split.by = obj.split.by, threads = threads)
    }

    select.bool = unlist(lapply(obj.l, function(x){ncol(x) >= min.cells.per.sample}))
    print(table(select.bool))
    obj.l = obj.l[select.bool]
    length(obj.l)

    print("HVG")
    obj.l = parallel::mclapply(obj.l, function(x) {
      x = x[!rownames(x) %in% bl, ]
      x = FindVariableFeatures(x, selection.method = "vst",  assay = .assay, verbose = FALSE)
      x
    }, mc.cores = threads)

    features = SelectIntegrationFeatures(object.list = obj.l, nfeatures = no.ftrs)

    VariableFeatures(obj) = features
    rm(obj.l); gc()

  } else if (!is.null(custom.features)) {
    VariableFeatures(obj) = custom.features
  } else {
    obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = no.ftrs, assay = .assay)
  }

  if (cc.regr == T) {
    # obj = ScaleData(obj, vars.to.regress = c("S.Score", "G2M.Score"), assay = .assay)
    obj = ScaleData(obj, vars.to.regress = c("Perc_of_mito_genes"), assay = .assay)
  } else {
    obj = ScaleData(obj, assay = .assay)
  }

  obj = RunPCA(obj, assay = .assay)
  plot(ElbowPlot(obj, ndims = 50))

  print(paste0("### PCs used for harmony clustering and dim reduc: ", .nbr.dims, " ###"))

  if(run.integration == T){
    obj = RunHarmony(
      obj, group.by.vars = harmony.group.vars, # theta = c(2,3),
      reduction.use ='pca', dims.use = 1:.nbr.dims, max_iter = 15, ncores = threads
    )
    comp.wrk = 'harmony'
  } else {
    comp.wrk = 'pca'
  }

  if (do.cluster == T) {
    print("Find clusters")
    obj = FindNeighbors(obj, reduction = comp.wrk, dims = 1:.nbr.dims)

    reso = seq(0,max.cl,.1)
    names(reso) = reso
    suppressWarnings({
      suppressMessages({
        findclusters.res = parallel::mclapply(reso, function(x) {
          FindClusters(obj, resolution = x)@meta.data[, "seurat_clusters", drop = F]
        }, mc.cores = length(reso))
      })
    })
    res.names = names(findclusters.res)
    findclusters.res = do.call("cbind", findclusters.res)
    colnames(findclusters.res) = paste0("RNA_snn_res.", res.names)
    stopifnot(identical(rownames(obj@meta.data), rownames(findclusters.res)))
    obj = AddMetaData(obj, findclusters.res)
  }
  if (do.dimreduc == T) {
    # print("tSNE")
    # obj = RunTSNE(
    #   obj, reduction = comp.wrk, dims = 1:dmreduc.dims, seed.use = 1234,
    #   nthreads = threads, tsne.method = "FIt-SNE"
    # )
    print("UMAP")
    obj = RunUMAP(
      obj, reduction = comp.wrk, dims = 1:dmreduc.dims, seed.use = 1234,
      min.dist = min.dist, n.neighbors = n.neighbors
    )
  }

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)

  return(obj)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# FindMarkers: Wilcoxon Test
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
run_wilx = function(
    obj = NULL,
    target = NULL,
    contrast.group = NULL,
    contrast = NULL,
    min.cells = 10,
    slot = "data",
    assay = "RNA",
    rm.var.chains = F,
    padj.thresh = 0.05,
    fc.thresh = 1.5,
    min.pct.thres = 0
){

  DefaultAssay(obj) = assay

  obj = DietSeurat(
    obj,
    counts = TRUE,
    data = T,
    scale.data = FALSE,
    features = NULL,
    assays = assay,
    dimreducs = F,
    graphs = F
  )

  obj@meta.data = droplevels(obj@meta.data)
  obj@meta.data$RM = is.na(obj@meta.data[[target]])
  obj = subset(obj, subset = RM == F)

  obj = obj[, obj@meta.data[[contrast.group]] %in% contrast]
  obj@meta.data = droplevels(obj@meta.data)
  obj.l = Seurat::SplitObject(obj, split.by = target)

  tbl = table(obj@meta.data[[contrast.group]], obj@meta.data[[target]])
  # Filter out celltypes with less than x cells in one group
  obj.l = obj.l[colnames(tbl)[colSums(tbl >= min.cells) == 2]]

  t = names(obj.l)
  options(warn = 1)
  bpparam = BiocParallel::MulticoreParam(workers = length(t))

  wil.res = BiocParallel::bplapply(t, function(x) {

    markers = Seurat::FindMarkers(
      object = obj.l[[x]],
      ident.1 = contrast[1], ident.2 = contrast[2],
      group.by = contrast.group, assay = assay,
      logfc.threshold = 0, min.pct = min.pct.thres
    )

    markers$feature = rownames(markers)
    markers$cluster = x
    markers$group.1 = contrast[1]
    rownames(markers) = paste0(markers$feature, "_", markers$cluster)
    markers

  }, BPPARAM = bpparam)

  wil.res = do.call("rbind", wil.res)
  if (rm.var.chains == T) {
    wil.res = wil.res[!grepl('^IGHV|^IGK|^IGL|^IGL|^TR(B|A)V|JCHAIN', wil.res$feature), ]
  }
  wil.res = wil.res[order(wil.res$p_val_adj, decreasing = F), ]
  wil.res$significant = (wil.res$p_val_adj < padj.thresh) & (abs(wil.res$avg_log2FC) > log2(fc.thresh))
  wil.res
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Presto: Wilcoxon Test: Cell identity marker
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
run_wilx_ct_marker = function(
    obj = NULL,
    target = NULL,
    min.cells = 10,
    slot = "data",
    assay = "RNA",
    downsample = F,
    cells = 1000,
    padj.thresh = 0.05,
    fc.thresh = 1.5
){

  DefaultAssay(obj) = assay

  obj = DietSeurat(
    obj,
    counts = TRUE,
    data = T,
    scale.data = FALSE,
    features = NULL,
    assays = assay,
    dimreducs = F,
    graphs = F
  )

  obj@meta.data = droplevels(obj@meta.data)
  obj@meta.data$RM = is.na(obj@meta.data[[target]])
  obj = subset(obj, subset = RM == F)
  obj@meta.data = droplevels(obj@meta.data)

  if(downsample == T) {
    Idents(obj) = target
    obj = subset(obj, downsample = cells)
  }

  t = unique(as.character(obj@meta.data[[target]]))
  bpparam = BiocParallel::MulticoreParam(workers = length(t))
  wil.res = BiocParallel::bplapply(t, function(x) {

    obj@meta.data$contrast = ifelse(obj@meta.data[[target]] == x, x, "other")
    markers = presto::wilcoxauc(
      obj, group_by = "contrast",
      seurat_assay = assay, assay = slot,
      groups_use = c(x, "other")
    )

    markers = subset(markers, group == x)
    markers$padj = p.adjust(markers$padj, method = "bonferroni")
    rownames(markers) = markers$feature
    markers

  }, BPPARAM = bpparam)

  wil.res = do.call("rbind", wil.res)
  wil.res = wil.res[order(abs(wil.res$logFC), decreasing = T), ]
  wil.res$significant = (wil.res$padj < padj.thresh) & (abs(wil.res$logFC) > log(fc.thresh))
  wil.res
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# FindAllMarker with DOR
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
wlx_test_dor = function(
    obj = NULL,
    se.grp = NULL,
    downsample = TRUE,
    downsample.nbr = 1000,
    only.pos = T,
    min.pct = 0
) {

  # obj@meta.data = droplevels(obj@meta.data)
  # Idents(obj) = se.grp
  # res.wlx = FindAllMarkers(
  #   object = obj, only.pos = only.pos, min.diff.pct = 0, logfc.threshold = 0,
  #   return.thresh = 1, min.pct = 0
  # )
  # # table(res.wlx$cluster)
  # grps = as.character(unique(res.wlx$cluster))

  if (downsample == TRUE) {
    print(paste0("Downsample identity to ", downsample.nbr))
    obj.sub = subset(x = obj, downsample = downsample.nbr)
  } else {
    print("No Downsampling")
    obj.sub = obj
  }

  obj.sub@meta.data = droplevels(obj.sub@meta.data)
  Idents(obj.sub) = se.grp
  res.wlx = FindAllMarkers(
    object = obj.sub, only.pos = only.pos, min.diff.pct = 0, logfc.threshold = 0,
    return.thresh = 1, min.pct = min.pct
  )
  # table(res.wlx$cluster)
  grps = as.character(unique(res.wlx$cluster))

  # i = "0"
  l = list()
  for (i in grps) {

    grps.s = subset(res.wlx, cluster == i)

    grp.cells = rownames(obj.sub@meta.data %>% filter(.data[[se.grp]] == i))

    num.cellsInGroup = length(grp.cells)
    num.cellsOutGroup = nrow(obj.sub@meta.data) - length(grp.cells)

    # get number of cells within & outside the group with these genes
    num.TruePos = rowSums(obj.sub@assays$RNA@data[grps.s$gene , grp.cells] > 0)
    num.FalsePos = rowSums(obj.sub@assays$RNA@data[grps.s$gene, !colnames(obj.sub@assays$RNA@data) %in% grp.cells] > 0)

    # get number of cells in and outside group without genes
    num.FalseNeg = num.cellsInGroup - num.TruePos
    num.TrueNeg <- num.cellsOutGroup - num.FalsePos

    # use these values to calculate log(DOR) w/ a pseudocount of 0.5 to avoid +/- infinity values
    logDOR = log((num.TruePos+0.5)/(num.FalsePos+0.5)/((num.FalseNeg+0.5)/(num.TrueNeg+0.5)))

    stopifnot(identical(names(logDOR), grps.s$gene))

    grps.s$logDOR = logDOR
    l[[i]] = grps.s
  }

  if(nrow(res.wlx) > 0) {
    res = do.call("rbind", l) %>% arrange(cluster, p_val_adj, desc(avg_log2FC))
    return(res)
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Add column in metadata wether CD4/CD8, CAR are present
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cd4cd8_car_present = function(
  obj
){
  cd8cd4 = FetchData(obj, c("CD8A", "CD8B", "CD4"), slot = "counts")
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

  if("ciltacel" %in% rownames(GetAssayData(obj, slot = c("counts"), assay = "RNA"))){
    car.ftr = FetchData(obj, c("ciltacel"), slot = "counts")
    obj$CAR_BY_EXPRS = as.factor(car.ftr$ciltacel > 0)
  } else {
    obj$CAR_BY_EXPRS = FALSE
    obj$CAR_BY_EXPRS = as.factor(obj$CAR_BY_EXPRS)
  }

  obj
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Track number of cell (pre-processing
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
count_cells_per_sample = function(obj = NULL, count.base = NULL, col.name = NULL){

  l = lapply(obj, function(x){
    x@meta.data %>% dplyr::select(STUDY, orig.ident)
  })
  df = do.call("rbind", l)
  df = df %>%  dplyr::count(STUDY, orig.ident)
  if(is.null(count.base)) {
    return(df)
  } else {
    count.base[[col.name]] = df$n[match(count.base$orig.ident, df$orig.ident)]
    return(count.base)
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# get metadata from Seurat
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
get_metadata <- function(obj, ..., embedding = names(obj@reductions), nbr.dim = 2) {

  res = as_tibble(obj@meta.data, rownames = "cell")

  if (!is.null(embedding)) {
    if (any(!embedding %in% names(obj@reductions))) {
      stop(paste0(embedding, " not found in seurat object\n"), call. = FALSE)
    }
    embed_dat = purrr::map(names(obj@reductions), ~obj@reductions[[.x]]@cell.embeddings[, 1:nbr.dim]) %>%
      do.call(cbind, .) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("cell")

    res = dplyr::left_join(res, embed_dat, by = "cell")
  }

  if (length(list(...)) > 0) {
    cols_to_get <- setdiff(..., colnames(obj@meta.data))
    if (length(cols_to_get) > 0) {
      res = Seurat::FetchData(obj, vars = cols_to_get) %>%
        tibble::rownames_to_column("cell") %>%
        dplyr::left_join(res, ., by = "cell")
    }
  }
  res
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Split Seurat object (BiocParallel)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Split_Object = function(object, split.by = "orig.ident", threads = 5) {

  library(parallel)
  library(BiocParallel)

  groupings <- FetchData(object = object, vars = split.by)[, 1]
  groupings <- unique(x = as.character(x = groupings))
  names(groupings) = groupings

  if (is.null(threads)) {
    bpparam = BiocParallel::MulticoreParam(workers = length(groupings))
  } else {
    bpparam = BiocParallel::MulticoreParam(workers = threads)
  }

  obj.list = BiocParallel::bplapply(groupings, function(grp) {
    cells <- which(x = object[[split.by, drop = TRUE]] == grp)
    cells <- colnames(x = object)[cells]
    se = subset(x = object, cells = cells)
    se@meta.data = droplevels(se@meta.data)
    se
  }, BPPARAM = bpparam)

  return(obj.list)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# flag levels of significance
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
add_signif <- function(
    data, p.col = NULL, output.col = NULL,
    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05,  1),
    symbols = c("****", "***", "**", "*",  "ns"),
    pval.relax = F
){

  if(pval.relax == T) {
    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1)
    symbols = c("*****", "****", "***", "**", "*",  "")
  }

  if(is.null(output.col)) {
    output.col <- paste0(p.col, ".signif")
  }
  .p.values <- data %>% pull(!!p.col)
  if(all(is.na(.p.values))) {
    .p.signif <- rep("", length(.p.values))
  }
  else{
    .p.signif <- .p.values %>%
      stats::symnum(cutpoints = cutpoints, symbols = symbols, na = "") %>%
      as.character()
  }
  data %>%
    dplyr::mutate(!!output.col := .p.signif)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Write DGEA results to xlsx
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
write_to_xlsx = function(
    df = NULL,
    sheet = "sheet",
    filename = "Supplementaltable3.xlsx",
    wb = NULL
) {
  library("openxlsx")

  if(is.null(wb)){
    wb <- createWorkbook()
  }

  addWorksheet(wb, sheet)
  writeData(
    wb, sheet,
    df %>% dplyr::select(
      "Gene_symbol" = feature, "LogFC" = logFC,
      "Pvalue" = pval, "FDR" = padj, "Cell_identity" = cluster
    ) %>%
      dplyr::arrange(Cell_identity, desc(LogFC)),
    startRow = 1, startCol = 1)
  saveWorkbook(wb, filename, overwrite = T)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# UCell: Gene set enrichment
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ucell_enrich = function(se.w = NULL, cust.gene.sets = NULL){

  # se.w = integration(
  #   obj = se.w,
  #   no.ftrs = 1000,
  #   threads = 20,
  #   .nbr.dims = 10,
  #   run.integration = F,
  #   do.cluster = F
  # )

  se.w = AddModuleScore_UCell(
    se.w,
    features = cust.gene.sets,
    assay = "RNA", slot = "data",
    ncores = 20, force.gc = T
  )
  ucell.res = se.w@meta.data[, grepl("UCell$", colnames(se.w@meta.data), ignore.case = T), drop = F]

  # ucell.res.smooth = UCell::SmoothKNN(
  #   obj = se.w,
  #   signature.names = colnames(ucell.res),
  #   reduction="pca", k=10
  # )
  # ucell.res.smooth = ucell.res.smooth@meta.data[, grepl("_UCell_kNN", colnames(ucell.res.smooth@meta.data))]

  colnames(ucell.res) = gsub("_UCell", "", colnames(ucell.res))
  ucell.res = ucell.res[, names(cust.gene.sets), drop = F]

  # colnames(ucell.res.smooth) = gsub("_UCell_kNN", "", colnames(ucell.res.smooth))
  # ucell.res.smooth = ucell.res.smooth[, colnames(ucell.res), drop = F]

  stopifnot(identical(rownames(ucell.res), colnames(se.w)))
  se.w[['custom_UCell_score']] = Seurat::CreateAssayObject(t(ucell.res))
  # slot(object = se.w[["custom_UCell_score"]], name = 'counts') <- new(Class = 'matrix')
  # se.w[['custom_UCell_smooth_score']] = Seurat::CreateAssayObject(t(ucell.res.smooth))
  # # slot(object = se.w[["custom_UCell_smooth_score"]], name = 'counts') <- new(Class = 'matrix')

  se.w@meta.data = se.w@meta.data[ , !grepl("UCell|nCount_UCell|nFeature_UCell", colnames(se.w@meta.data), ignore.case = T)]
  se.w

}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Ave Expression and % of present cells
# https://divingintogeneticsandgenomics.com/post/how-to-make-a-multi-group-dotplot-for-single-cell-rnaseq-data/
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
GetMatrixFromSeuratByGroupMulti<- function(
  obj,
  features,
  group1,
  group2,
  .assay = "RNA",
  perc_expr_thres = 0
){
  exp_mat<- obj@assays[[.assay]]@data[features, ,drop=FALSE]
  if(.assay == "ADT") {
    count_mat<- obj@assays[[.assay]]@data[features,,drop=FALSE ]
  } else {
    count_mat<- obj@assays[[.assay]]@counts[features,,drop=FALSE ]
  }


  meta<- obj@meta.data %>%
    tibble::rownames_to_column(var = "cell")

  # get the average expression matrix
  exp_df<- as.matrix(exp_mat) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var="gene") %>%
    tidyr::pivot_longer(!gene, names_to = "cell", values_to = "expression") %>%
    left_join(meta) %>%
    group_by(gene,{{group1}}, {{group2}}) %>%
    summarise(average_expression = mean(expression)) %>%
    # the trick is to make the data wider in columns: cell_type|group
    tidyr::pivot_wider(names_from = c({{group1}}, {{group2}}),
                       values_from= average_expression,
                       names_sep="|")

  # convert to a matrix
  exp_mat<- exp_df[, -1] %>% as.matrix()
  rownames(exp_mat)<- exp_df$gene

  # get percentage of positive cells matrix
  count_df<- as.matrix(count_mat) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var="gene") %>%
    tidyr::pivot_longer(!gene, names_to = "cell", values_to = "count") %>%
    left_join(meta) %>%
    group_by(gene, {{group1}}, {{group2}}) %>%
    summarise(percentage = mean(count > perc_expr_thres)) %>%
    tidyr::pivot_wider(names_from = c({{group1}}, {{group2}}),
                       values_from= percentage,
                       names_sep="|")

  percent_mat<- count_df[, -1] %>% as.matrix()
  rownames(percent_mat)<- count_df$gene

  if (!identical(dim(exp_mat), dim(percent_mat))) {
    stop("the dimension of the two matrice should be the same!")
  }

  if(! all.equal(colnames(exp_mat), colnames(percent_mat))) {
    stop("column names of the two matrice should be the same!")
  }

  if(! all.equal(rownames(exp_mat), rownames(percent_mat))) {
    stop("column names of the two matrice should be the same!")
  }
  return(list(exp_mat = exp_mat, percent_mat = percent_mat))
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# get Spacial object
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
get_special_se = function(work.path, id = NULL){
  samples =  file.path(work.path, "filtered_feature_bc_matrix.h5")
  imgs = file.path(work.path,  "spatial", "tissue_lowres_image.png")
  spotfiles = file.path(work.path,"spatial", "tissue_positions.csv")
  json = file.path(work.path, "spatial", "scalefactors_json.json")
  infoTable <- tibble(samples, imgs, spotfiles, json, sample_id = id)
  se = semla::ReadVisiumData(infoTable)
  se = semla::LoadImages(se) %>% Seurat::NormalizeData()
  se
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# clean Clonotype: Clonotypes must not have CD4 and CD8 cells
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
clean_clonotypes = function(
    obj = NULL,
    celltype = "celltype_short_3",
    clone_id = "CTstrict"
){

  obj$celltype = obj[[celltype]]
  obj$clone_id = obj[[clone_id]]

  pd = obj@meta.data %>%
    dplyr::select(celltype, clone_id)
  pd$barcode = rownames(pd)
  cl.lin.diff =
    pd %>%
    dplyr::filter(!is.na(clone_id)) %>%
    dplyr::filter(grepl("CD4|CD8", celltype)) %>%
    dplyr::group_by(clone_id) %>%
    dplyr::mutate(clone_size = n()) %>%
    dplyr::mutate(cd4 = grepl("CD4", celltype)) %>%
    dplyr::mutate(cd8 = grepl("CD8", celltype)) %>%
    dplyr::mutate(nbr_cd4 = sum(cd4)) %>%
    dplyr::mutate(nbr_cd8 = sum(cd8)) %>%
    dplyr::filter(nbr_cd4 > 0 & nbr_cd8 > 0)

  rm.cl.1 = cl.lin.diff %>%
    dplyr::filter(clone_size <= 5) %>%
    dplyr::pull(barcode) %>%
    unique()

  cl.lin.diff = cl.lin.diff[!cl.lin.diff$barcode %in% rm.cl.1, ]

  if(nrow(cl.lin.diff) == 0){
    return(c(rm.cl.1))
  }

  cl.lin.diff$ratio = NA
  for (i in unique(cl.lin.diff$clone_id)) {
    cl.sub = cl.lin.diff[cl.lin.diff$clone_id == i, ]
    cl.sub = cl.sub[!duplicated(cl.sub$clone_id), ]
    .max = names(which.max(cl.sub[, c("nbr_cd4", "nbr_cd8")]))
    .min = names(which.min(cl.sub[, c("nbr_cd4", "nbr_cd8")]))
    cl.lin.diff[cl.lin.diff$clone_id == i, ]$ratio = cl.sub[[.max]] / cl.sub[[.min]]
  }

  rm.cl.2 = cl.lin.diff %>%
    dplyr::filter(ratio <= 3) %>%
    data.frame() %>%
    dplyr::pull(barcode) %>%
    unique()

  # tmp = cl.lin.diff[cl.lin.diff$barcode %in% rm.cl.2, ]
  # tmp[!duplicated(tmp$clone_id), ] %>% data.frame() %>% select(nbr_cd4, nbr_cd8, ratio)

  cl.lin.diff = cl.lin.diff[!cl.lin.diff$barcode %in% rm.cl.2, ]

  if(nrow(cl.lin.diff) == 0){
    return(c(rm.cl.1, rm.cl.2))
  }

  # cl.lin.diff[!duplicated(cl.lin.diff$clone_id), ] %>% data.frame() %>% select(nbr_cd4, nbr_cd8, ratio)

  rm.cl.3 = cl.lin.diff %>%
    dplyr::group_by(clone_id, celltype) %>%
    dplyr::mutate(lin = n()) %>%
    dplyr::group_by(clone_id) %>%
    dplyr::filter(lin == min(lin)) %>%
    dplyr::pull(barcode) %>%
    unique()

  return(c(rm.cl.1, rm.cl.2, rm.cl.3))
  # cl.lin.diff = cl.lin.diff[!cl.lin.diff$barcode %in% rm.cl.3, ]
  # cl.lin.diff %>%
  #   dplyr::group_by(clone_id) %>%
  #   dplyr::mutate(clone_size = n()) %>%
  #   dplyr::mutate(cd4 = grepl("CD4", celltype)) %>%
  #   dplyr::mutate(cd8 = grepl("CD8", celltype)) %>%
  #   dplyr::mutate(nbr_cd4 = sum(cd4)) %>%
  #   dplyr::mutate(nbr_cd8 = sum(cd8)) %>%
  #   dplyr::filter(nbr_cd4 > 0 & nbr_cd8 > 0) %>%
  #   data.frame()

}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Assign VDJ
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
assign_vdj = function(obj, vdj = "vdj_t", batch = NULL, present.bool = TRUE){

  stopifnot(vdj %in% c("vdj_t", "vdj_b"))

  library(scRepertoire)

  l = list()
  for (i in batch) {
    print(i)
    cellranger.dirs = list.dirs(
      path = i, full.names = T, recursive = F
    )
    cellranger.samples = basename(cellranger.dirs)
    fltrd.vdj = paste0(
      cellranger.dirs, "/outs/per_sample_outs/", cellranger.samples, "/", vdj,
      "/filtered_contig_annotations.csv"
    )
    names(fltrd.vdj) = gsub("multi_", "", cellranger.samples)
    print(length(fltrd.vdj))
    l[[i]] = fltrd.vdj
  }
  names(l) = NULL
  fltrd.vdj = unlist(l, use.names = T)

  ###

  contig_list <- lapply(fltrd.vdj, function(x) {
    tryCatch(read.csv(x), error=function(e) NULL)
  })
  length(contig_list)
  print(paste0("Empty file: ", names(lengths(contig_list)[lengths(contig_list) == 0])))
  contig_list = contig_list[lengths(contig_list) != 0]

  if(vdj == "vdj_t"){
    combined <- scRepertoire::combineTCR(
      contig_list,
      samples = paste0(names(contig_list))
    )
  } else {
    combined <- scRepertoire::combineBCR(
      contig_list,
      samples = paste0(names(contig_list))
    )
  }

  combined = data.table::rbindlist(combined) %>% data.frame()

  if(present.bool == TRUE){
    if(vdj == "vdj_t"){
      obj$VDJ_T_AVAIL = combined$CTnt[match(rownames(obj@meta.data), combined$barcode)]
      obj$VDJ_T_AVAIL = ifelse(is.na(obj$VDJ_T_AVAIL), FALSE, TRUE)
      print(table(obj$VDJ_T_AVAIL))
    } else {
      obj$VDJ_B_AVAIL = combined$CTnt[match(rownames(obj@meta.data), combined$barcode)]
      obj$VDJ_B_AVAIL = ifelse(is.na(obj$VDJ_B_AVAIL), FALSE, TRUE)
      print(table(obj$VDJ_B_AVAIL))
    }
    obj
  } else {

    combined$sample = obj$orig.ident[match(combined$barcode, rownames(obj@meta.data))]
    combined = combined[!is.na(combined$sample), ]
    combined = split(combined, combined$sample)

    max.clonotypes = max(unlist(lapply(combined, function(x){max(unname(table(x$CTstrict)))})))
    obj <- combineExpression(
      combined, obj,
      cloneCall = "strict",
      group.by = "sample",
      proportion = F,
      cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=max.clonotypes)
    )
    obj$cloneSize = droplevels(obj$cloneSize)
    obj
    # df = se.meta@meta.data[, (length(colnames(se.meta@meta.data))-6):length(colnames(se.meta@meta.data))]
    # saveRDS(df, paste0(manifest$meta$work, "integration/clonotypes_b.Rds"))
  }
}

