.cran_packages = c(
  "yaml", "ggplot2","reshape2", "dplyr", "naturalsort", "devtools", "scales",
  "stringr", "Seurat", "tibble", "tidyr", "forcats", "scCustomize",
  "rlang", "remotes", "patchwork", "cowplot", "ggrepel", "scico",
  "ggpubr", "circlize", "gridtext", "data.table", "msigdbr"
)
.bioc_packages = c(
  "simplifyEnrichment", "ComplexHeatmap", "clusterProfiler", "fgsea"
)

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
# Remove similar significantly enriched GO Term
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
simplify_go = function(
    res,
    go.id = "ID",
    cutoff = go.sim.thresh,
    by= "p.adjust",
    ontology = "BP"
) {

  semData = GOSemSim:::godata(ont = "BP")
  sim = GOSemSim:::mgoSim(
    res[[go.id]], res[[go.id]],
    semData = semData,
    measure = "Wang",
    combine=NULL
  )

  go1 <- go2 <- similarity <- NULL

  sim.df <- as.data.frame(sim)
  sim.df$go1 <- row.names(sim.df)
  sim.df <- gather(sim.df, go2, similarity, -go1)

  sim.df <- sim.df[!is.na(sim.df$similarity),]


  res.sub = res %>% dplyr::select(.data[[go.id]], .data[[by]])
  sim.df <- merge(sim.df, res.sub, by.x = "go1", by.y = go.id)
  sim.df$go2 <- as.character(sim.df$go2)

  ID <- res[[go.id]]

  GO_to_remove <- character()
  for (i in seq_along(ID)) {
    ii <- which(sim.df$go2 == ID[i] & sim.df$similarity > cutoff)
    if (length(ii) < 2)
      next

    sim_subset <- sim.df[ii,]

    jj <- which(sim_subset[, by] == base::min(sim_subset[, by]))
    GO_to_remove <- c(GO_to_remove, sim_subset$go1[-jj]) %>% unique
  }

  res[!res[[go.id]] %in% GO_to_remove, ]
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ORA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
run_nmf_ora = function (
    genes,
    universe = NULL,
    category = "C5",
    subcategory = NULL,
    db.sub = NULL,
    padj.thr = 0.05,
    go.sim.thresh = .7
)
{

  if (
    !requireNamespace("fgsea", quietly = TRUE) |
    !requireNamespace("msigdbr", quietly = TRUE) |
    !requireNamespace("simplifyEnrichment", quietly = TRUE)
  ) {
    stop(
      "Function 'runGSEA' requires the 'fgsea' and 'msigdbr' packages. \n Please
      install them.", call. = FALSE
    )
  }

  if(category != "C5"){subcategory = NULL}
  if (any(duplicated(genes))) {
    genes <- genes[!duplicated(genes)]
  }
  msig_df = msigdbr::msigdbr(
    species = "Homo sapiens", category = category, subcategory = subcategory
  )

  if(!is.null(db.sub)) {
    msig_df = msig_df[grepl(db.sub, msig_df$gs_name), ]
  }

  # mir = readRDS("data/metadata/signatures/mir_gene_set_collection.Rds")
  # msig_df = mir[grepl("CD8", mir$TermID), ]
  # msig_df = msig_df[!grepl("Cytopus", msig_df$TermID), ]
  # msig_df$TermID = gsub("CD8-", "", msig_df$TermID)
  # msig_df$TermID = gsub(" .+", "", msig_df$TermID)
  # colnames(msig_df) = c("gs_name", "gene_symbol")

  msig_list = split(x = msig_df$gene_symbol, f = msig_df$gs_name)
  msig_list = msig_list[lengths(msig_list) < 500]
  msig_list = msig_list[lengths(msig_list) > 10]

  fgRes = fgsea::fora(pathways = msig_list, genes = genes,  universe = universe)
  fgRes = fgRes[fgRes$padj < padj.thr, ]
  fgRes$go_id = msig_df$gs_exact_source[match(fgRes$pathway, msig_df$gs_name)]
  if(category == "C5" & !is.null(subcategory)){
    fgRes =  simplify_go(
      fgRes, by = "padj", go.id = "go_id",
      ontology = subcategory, cutoff = go.sim.thresh
    )
  }

  return(fgRes)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Bubble plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ora_bubble = function(
    gsea.res = NULL,
    ftrs.list = NULL,
    nbr.tops = 5,
    min.genes = 3,
    font.size = 8,
    dot.range = c(1.5, 4.5),
    max.value = NULL,
    term.length = 50,
    sort.by.padj = FALSE,
    barwidth = unit(.5, 'lines'),
    barheight = unit(5, 'lines'),
    y.text.size = rel(1)
) {

  gsea.df = data.table::rbindlist(gsea.res, idcol = "MP")
  gsea.df = gsea.df[gsea.df$overlap > min.genes,]
  gsea.df = data.frame(gsea.df)
  gsea.df$richfactor = gsea.df$overlap / gsea.df$size

  gsea.df$pathway = gsub("GOBP_", "", gsea.df$pathway)
  gsea.df$pathway = gsub("BIOLOGICAL_PROCESS", "BP", gsea.df$pathway)
  gsea.df$pathway = gsub("BIOLOGICAL_PROCESS", "BP", gsea.df$pathway)
  # gsea.df$pathway = gsub("REGULATION", "REG.", gsea.df$pathway)
  gsea.df$pathway = gsub("IMMUNOGLOBULIN", "IG", gsea.df$pathway)
  gsea.df$pathway = gsub("RESPONSE", "RESP.", gsea.df$pathway)
  gsea.df$pathway = gsub("ANTIGEN_PROCESSING_AND_PRESENTATION", "APAP", gsea.df$pathway)
  gsea.df$pathway = gsub("TUMOR_NECROSIS_FACTOR", "TNF", gsea.df$pathway)
  gsea.df$pathway = gsub("MHC_CLASS_II", "MHC_II", gsea.df$pathway)
  gsea.df$pathway = gsub("MHC_CLASS_I", "MHC_I", gsea.df$pathway)
  gsea.df$pathway = gsub("REACTOME_", "", gsea.df$pathway)
  gsea.df$pathway = gsub("RESPIRATORY_ELECTRON_TRANSPORT", "ETC", gsea.df$pathway)
  gsea.df$pathway = gsub("INTERFERON", "IFN", gsea.df$pathway)
  gsea.df$pathway = gsub("SIGNALING_PATHWAY", "PATHWAY", gsea.df$pathway)
  gsea.df$pathway = gsub("INTERLEUKIN", "IL", gsea.df$pathway)
  # gsea.df$pathway = gsub("PRODUCTION", "PROD.", gsea.df$pathway)

  if(sort.by.padj == T){
    tops = gsea.df %>%
      dplyr::group_by(MP) %>%
      dplyr::arrange(padj) %>%
      dplyr::slice_head(n = nbr.tops)
  } else {
    tops = gsea.df %>%
      dplyr::group_by(MP) %>%
      dplyr::arrange(desc(richfactor)) %>%
      dplyr::slice_head(n = nbr.tops)

  }
  gsea.df = subset(gsea.df, pathway %in% tops$pathway)

  lvls = gsea.res[[1]]$pathway
  lvls = c(lvls, setdiff(gsea.df$pathway, lvls))
  gsea.df$pathway = factor(gsea.df$pathway, levels = rev(lvls))
  gsea.df$MP = factor(gsea.df$MP, levels = naturalsort(unique(gsea.df$MP)))

  if(!is.null(ftrs.list)) {
    gsea.df$zscore = NA
    for (i in 1:nrow(gsea.df)) {
      ct = as.character(gsea.df[i, ]$MP)
      ftrs = unlist(gsea.df[i, ]$overlapGenes)
      ftrs.lfc = names(ftrs.list[[ct]][(ftrs.list[[ct]]) %in% ftrs])
      gs.ftrs = as.numeric(gsea.df[i, ]$size)
      zscore = (length(ftrs.lfc[ftrs.lfc > 0]) - length(ftrs.lfc[ftrs.lfc < 0])) / sqrt(gs.ftrs)
      gsea.df[i, ]$zscore = zscore
    }
    gsea.df$dot_size = -log10(gsea.df$padj)
    fill.title = "z-score"

    if(is.null(max.value)){
      max.value = max(abs(gsea.df$zscore))
    } else {
      gsea.df$zscore[gsea.df$zscore > max.value] = max.value
      gsea.df$zscore[gsea.df$zscore < -max.value] = -max.value
    }
    gsea.df$fill = gsea.df$zscore
  } else {
    gsea.df$fill = -log10(gsea.df$padj)
    fill.title = "-Log10(FDR)"
  }

  default_labeller <- function(n) {
    function(str){
      str <- gsub("_", " ", str)
      yulab.utils::str_wrap(str, n)
    }
  }


  pl = ggplot(gsea.df, aes(x = MP, y = pathway))
  if(!is.null(ftrs.list)) {
    pl = pl + geom_point(
      aes(fill = fill, size = dot_size), colour="black", pch=21, stroke = .3
    ) +
      scale_size(range = dot.range) +
      scale_fill_scico(
        palette = "vik", midpoint = 0, begin = .1, end = .9,
        limits = c(-max.value, max.value),
        breaks = pretty_breaks(n = 3)
      )
  } else {
    pl = pl + geom_point(
      aes(fill = fill), colour="black", pch=21, stroke = .3, size = 3
    ) +
      scale_fill_gradientn(
        colors = scico(30, palette = 'bilbao', direction = -1, end = .9, begin = .1),
        breaks = pretty_breaks(n = 3),
      )
  }
  pl = pl + mytheme(base_size = font.size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "plain", size = rel(1)),
      axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = rel(1)),
      axis.text.y = element_text(size = y.text.size),
      legend.position = "right"
    ) +
    guides(
      fill = guide_colorbar(
        title =  fill.title,
        barwidth = barwidth,
        barheight = barheight,
        order = 1, ticks.linewidth = 1.5/.pt
      ),
      size = guide_legend(title = "-Log10(FDR)", order = 2)
    ) +
    scale_y_discrete(labels = default_labeller(term.length)) +
    xlab(NULL) + ylab(NULL)

  pl
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Bar plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ora_barpl = function(
    gsea.res = NULL,
    ftrs.list = NULL,
    nbr.tops = 5,
    min.genes = 3,
    font.size = 8,
    max.value = NULL,
    term.length = 50,
    sort.by.padj = TRUE,
    barwidth = unit(.5, 'lines'),
    barheight = unit(5, 'lines'),
    bar.width = 1
) {

  gsea.df = data.table::rbindlist(gsea.res, idcol = "MP")
  gsea.df = gsea.df[gsea.df$overlap > min.genes,]
  gsea.df = data.frame(gsea.df)
  gsea.df$richfactor = gsea.df$overlap / gsea.df$size

  gsea.df$pathway = gsub("GOBP_", "", gsea.df$pathway)
  gsea.df$pathway = gsub("BIOLOGICAL_PROCESS", "BP", gsea.df$pathway)
  gsea.df$pathway = gsub("BIOLOGICAL_PROCESS", "BP", gsea.df$pathway)
  gsea.df$pathway = gsub("REGULATION", "REG.", gsea.df$pathway)
  gsea.df$pathway = gsub("IMMUNOGLOBULIN", "IG", gsea.df$pathway)
  gsea.df$pathway = gsub("RESPONSE", "RESP.", gsea.df$pathway)
  gsea.df$pathway = gsub("ANTIGEN_PROCESSING_AND_PRESENTATION", "APAP", gsea.df$pathway)
  gsea.df$pathway = gsub("TUMOR_NECROSIS_FACTOR", "TNF", gsea.df$pathway)
  gsea.df$pathway = gsub("MHC_CLASS_II", "MHC_II", gsea.df$pathway)
  gsea.df$pathway = gsub("MHC_CLASS_I", "MHC_I", gsea.df$pathway)
  gsea.df$pathway = gsub("REACTOME_", "", gsea.df$pathway)
  gsea.df$pathway = gsub("RESPIRATORY_ELECTRON_TRANSPORT", "ETC", gsea.df$pathway)
  gsea.df$pathway = gsub("INTERFERON", "IFN", gsea.df$pathway)
  gsea.df$pathway = gsub("SIGNALING_PATHWAY", "PATHWAY", gsea.df$pathway)
  gsea.df$pathway = gsub("INTERLEUKIN", "IL", gsea.df$pathway)
  # gsea.df$pathway = gsub("PRODUCTION", "PROD.", gsea.df$pathway)

  if(sort.by.padj == T){
    tops = gsea.df %>%
      dplyr::arrange(padj) %>%
      dplyr::slice_head(n = nbr.tops)
  } else {
    tops = gsea.df %>%
      dplyr::arrange(desc(richfactor)) %>%
      dplyr::slice_head(n = nbr.tops)

  }
  gsea.df = subset(gsea.df, pathway %in% tops$pathway)

  lvls = tops$pathway
  gsea.df$pathway = factor(gsea.df$pathway, levels = rev(lvls))

  if(!is.null(ftrs.list)) {
    gsea.df$zscore = NA
    for (i in 1:nrow(gsea.df)) {
      ct = as.character(gsea.df[i, ]$MP)
      ftrs = unlist(gsea.df[i, ]$overlapGenes)
      ftrs.lfc = names(ftrs.list[[ct]][(ftrs.list[[ct]]) %in% ftrs])
      gs.ftrs = as.numeric(gsea.df[i, ]$size)
      zscore = (length(ftrs.lfc[ftrs.lfc > 0]) - length(ftrs.lfc[ftrs.lfc < 0])) / sqrt(gs.ftrs)
      gsea.df[i, ]$zscore = zscore
    }
    gsea.df$dot_size = -log10(gsea.df$padj)
    fill.title = "z-score"

    if(is.null(max.value)){
      max.value = max(abs(gsea.df$zscore))
    } else {
      gsea.df$zscore[gsea.df$zscore > max.value] = max.value
      gsea.df$zscore[gsea.df$zscore < -max.value] = -max.value
    }
    gsea.df$fill = gsea.df$zscore
  }

  default_labeller <- function(n) {
    function(str){
      str <- gsub("_", " ", str)
      yulab.utils::str_wrap(str, n)
    }
  }

  if(sort.by.padj == TRUE){
    gsea.df$X_AXIS = -log10(gsea.df$padj)
    x_title = "-Log10(FDR)"
    if(is.null(ftrs.list)) {
      gsea.df$fill = gsea.df$richfactor
      fill.title = "Richfactor"
    }
  } else {
    gsea.df$X_AXIS = gsea.df$richfactor
    x_title = "Richfactor"
    if(is.null(ftrs.list)) {
      gsea.df$fill = -log10(gsea.df$padj)
      fill.title = "-Log10(FDR)"
    }
  }

  pl = ggplot(data = gsea.df, aes(x = pathway, y = X_AXIS)) +
    coord_flip()
  if(!is.null(ftrs.list)) {
    pl = pl + geom_bar(aes(fill = fill), stat = "identity", colour="black", linewidth = .2, width = bar.width) +
      scale_fill_scico(
        palette = "vik", midpoint = 0, begin = .1, end = .9,
        limits = c(-max.value, max.value),
        breaks = pretty_breaks(n = 4)
      )
  } else {
    pl = pl + geom_bar(aes(fill = fill), stat = "identity", colour="black", linewidth = .2, width = bar.width) +
      scale_fill_gradientn(
        colors = scico(30, palette = 'bilbao', direction = -1, end = .9, begin = .1),
        breaks = pretty_breaks(n = 4),
      )
  }
  pl =
    pl + mytheme(base_size = font.size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "plain", size = rel(1)),
      axis.text.x = element_text(size = rel(1)),
      legend.position = "right",
      legend.ticks.length = unit(0.05, 'cm')
    ) +
    guides(
      fill = guide_colorbar(
        title =  fill.title,
        barwidth = barwidth,
        barheight = barheight,
        order = 1, ticks.linewidth = 1.5/.pt
      )
    ) +
    scale_x_discrete(labels = default_labeller(term.length)) +
    xlab(NULL) + ylab(x_title)

  pl
}
