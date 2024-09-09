# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Violin plot (pre filtering) with cutoffs
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
qc_vln_plot_cell = function(
    obj = se.meta,
    .features = "nFeature_RNA",
    .group.by = "orig.ident",
    plot_title = "Genes Per Cell",
    x_axis_label = NULL,
    y_axis_label = "Features",
    low_cutoff = NULL,
    high_cutoff = NULL
){
  library(ggplot2)
  library(ggthemes)

  if(!is.list(obj)) {
    df = obj@meta.data %>% dplyr::select(.data[[.group.by]], .data[[.features]])
    df
  }
  if(class(obj) == "list") {
    l = lapply(obj, function(x){
      df = x@meta.data %>% dplyr::select(.data[[.group.by]], .data[[.features]])
    })
    df = do.call("rbind", l)
    df
  }
  if(class(obj) == "data.frame") {
    df = obj %>% dplyr::select(.data[[.group.by]], .data[[.features]])
    df
  }

  ggplot(data = df, mapping = aes(x = .data[[.group.by]], y = .data[[.features]])) +
    geom_violin(
      size = .1,
      width = 1,
      scale = "area",
      na.rm = TRUE
    ) +
    stat_summary(
      fun.min = function(z) { quantile(z,0.25) },
      fun.max = function(z) { quantile(z,0.75) },
      fun = median, colour = "#0077BB", size = .2) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle=45, vjust=1, hjust=1),
      axis.ticks.x = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(size = rel(1.2))
    ) +
    geom_hline(yintercept = c(low_cutoff, high_cutoff), linetype = "dashed", color = "#BB5566", size = .3) +
    xlab(x_axis_label) +
    ylab(y_axis_label) +
    ggtitle(plot_title)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
count_cells_per_sample = function(obj = NULL, count.base = NULL, col.name = NULL){

  l = lapply(obj, function(x){
    x@meta.data %>% dplyr::select(orig.ident)
  })
  df = do.call("rbind", l)
  df = df %>%  dplyr::count(orig.ident)
  if(is.null(count.base)) {
    return(df)
  } else {
    count.base[[col.name]] = df$n[match(count.base$orig.ident, df$orig.ident)]
    return(count.base)
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Dimension reduction: for features (assay)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dimreduc_features = function(
  .obj = NULL,
  features = NULL,
  .reduc = "tsne",
  pl.points = F,
  pt.size = .3,
  plot.grid = T,
  base.size = 8,
  .x.title = "UMAP 1",
  .y.title = "UMAP 2",
  .leg.title = NULL,
  .title.size = 1,
  .assay = "RNA",
  log.scale.counts = F,
  .quantile.fltr = T,
  .na_cutoff = NULL,
  min.max = NULL,
  legend.wh = c(.3, 4),
  order = T,
  .raster.scattermore = FALSE,
  .raster.scattermore.pixel = c(512,512),
  .raster.scattermore.pointsize = 0,
  .raster = FALSE,
  .raster.dpi = 150,
  .colors = rev(MetBrewer::met.brewer("Hokusai1",n=100))
) {

  DefaultAssay(.obj) = .assay

  library(scales)
  library(cowplot)
  library(scico)

  if(!is.null(.x.title) & !is.null(.y.title)){
    if(grepl("sne", .reduc, ignore.case = T)) {.x.title = "tSNE 1"; .y.title = "tSNE 2"}
    if(grepl("umap", .reduc, ignore.case = T)) {.x.title = "UMAP 1"; .y.title = "UMAP 2"}
  }

  features = features[features %in% rownames(.obj)]
  if(length(features) == 0) { stop("None of the genes are present in the data set")}

  if(log.scale.counts == T){
    exprs.sub = .obj@assays[[.assay]]@counts[features, , drop = F]
    exprs.sub = log10(exprs.sub + 1)
  } else {
    exprs.sub = .obj@assays[[.assay]]@data[features, , drop = F]
  }

  ftr.pl = list()
  for (i in 1:length(features)) {
    ftr = features[i]
    reduc = data.frame(.obj@reductions[[.reduc]]@cell.embeddings)
    colnames(reduc) = c("DIM_1", "DIM_2")
    reduc$EXPRS = exprs.sub[ftr, ]

    if(order == T) {
      reduc = reduc[order(reduc$EXPRS, decreasing = F), ] # For ggplot: highest values on the top
    } else {
      set.seed(1234)
      reduc = reduc %>% dplyr::sample_frac(1L, replace = FALSE) # permute rows randomly
    }

    ftr.exprs = reduc$EXPRS

    if(.quantile.fltr) {
      qu.max = quantile(ftr.exprs[ftr.exprs > 0], .999)
      ftr.exprs[ftr.exprs > qu.max] = qu.max
      reduc$EXPRS = ftr.exprs
    }

    if(!is.null(min.max)) {
      e.min = min.max[1]
      e.max = min.max[2]
      ftr.exprs[ftr.exprs > e.max] = e.max
    } else if(!is.null(.na_cutoff) & is.null(min.max)) {
      e.min = min(ftr.exprs[ftr.exprs > .na_cutoff])
      e.max = max(ftr.exprs)
    } else if (is.null(.na_cutoff) & is.null(min.max)) {
      e.min = min(ftr.exprs)
      e.max = max(ftr.exprs)
    }

    if(is.null(names(ftr))) {names(ftr) = ""}

    pl =
      ggplot(reduc, aes(x = DIM_1, y = DIM_2, color = EXPRS))
    if (.raster.scattermore == T) {
      pl = pl + scattermore::geom_scattermore(
        pointsize = .raster.scattermore.pointsize, pixels = .raster.scattermore.pixel
      )
    } else {
      if (pl.points == F) {
        pl = pl + geom_point(shape = ".", alpha = 1)
      } else {
        pl = pl + geom_point(size = pt.size)
      }
    }
    pl = pl +  scale_color_gradientn(
      colors = .colors,
      na.value = "#DDDDDD",
      limits = c(e.min, e.max),
      breaks = pretty_breaks(4)
    )
    pl = pl +  guides(
      color = guide_colorbar(
        title = .leg.title, title.hjust = 0, barwidth = unit(legend.wh[1],'lines'),
        barheight = unit(legend.wh[2], 'lines'), ticks.linewidth = 1.5/.pt
      )
    ) +
    {if(nchar(names(ftr)) == 0)ggtitle(ftr)} +
    {if(!nchar(names(ftr)) == 0)ggtitle(names(ftr))} +
    mytheme(base_size = base.size) +
    theme(
      aspect.ratio = 1,
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.spacing = unit(1, "lines"),
      plot.title = element_text(hjust = 0.5, colour = "black", size = rel(.title.size))
    ) +
    xlab(.x.title) + ylab(.y.title)

    if (.raster == T) {
      pl = ggrastr::rasterize(pl, layers='Point', dpi=.raster.dpi)
    }

    ftr.pl[[i]] = pl
  }

  if (plot.grid == T) {
    if(length(ftr.pl) <= 4) {ftr.pl = c(ftr.pl, vector(mode = "list", length = (4 - length(ftr.pl))))}
    if(length(ftr.pl) > 4 && length(ftr.pl) <= 8) {ftr.pl = c(ftr.pl, vector(mode = "list", length = (8 - length(ftr.pl))))}
    plot_grid(plotlist = ftr.pl, ncol = 4, scale = .95, align = "vh")
  } else {
    ftr.pl
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Spatial Featureplot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
spatial_features = function(
    .obj = NULL,
    features = NULL,
    .assay = "Spatial",
    slot = "data",
    log.counts = F,
    pl.points = T,
    pt.size = .3,
    plot.grid = T,
    base.size = 8,
    .leg.title = NULL,
    .title.size = 1,
    .quantile.fltr = T,
    .na_cutoff = NULL,
    min.max = NULL,
    legend.wh = c(.3, 4),
    .raster.scattermore = FALSE,
    .raster.scattermore.pixel = c(512,512),
    .raster.scattermore.pointsize = 0,
    .raster = FALSE,
    .raster.dpi = 150,
    aspect.ratio = 1,
    .colors = rev(MetBrewer::met.brewer("Hokusai1",n=100))
) {

  DefaultAssay(.obj) = .assay

  library(scales)
  library(cowplot)
  library(scico)

  features = features[features %in% rownames(.obj)]
  if(length(features) == 0) { stop("None of the genes are present in the data set")}

  ftr.pl = list()
  for (i in 1:length(features)) {
    ftr = features[i]
    data_use <- semla::GetStaffli(.obj)@meta_data %>%
      dplyr::bind_cols(Seurat::FetchData(.obj, vars = ftr, slot = slot, clean = FALSE))
    data_use$EXPRS = data_use[[ftr]]

    if(log.counts == T){
      data_use$EXPRS = log10(data_use$EXPRS + 1)
    }

    ftr.exprs = data_use$EXPRS



    if(.quantile.fltr) {
      qu.max = quantile(ftr.exprs[ftr.exprs > 0], .999)
      ftr.exprs[ftr.exprs > qu.max] = qu.max
      data_use$EXPRS = ftr.exprs
    }

    if(!is.null(min.max)) {
      e.min = min.max[1]
      e.max = min.max[2]
      ftr.exprs[ftr.exprs > e.max] = e.max
    } else if(!is.null(.na_cutoff) & is.null(min.max)) {
      e.min = min(ftr.exprs[ftr.exprs > .na_cutoff])
      e.max = max(ftr.exprs)
    } else if (is.null(.na_cutoff) & is.null(min.max)) {
      e.min = min(ftr.exprs)
      e.max = max(ftr.exprs)
    }

    if(is.null(names(ftr))) {names(ftr) = ""}

    pl =
      ggplot(data_use, aes(x = x, y = y, color = EXPRS)) +
      scale_y_reverse()
    if (.raster.scattermore == T) {
      pl = pl + scattermore::geom_scattermore(
        pointsize = .raster.scattermore.pointsize, pixels = .raster.scattermore.pixel
      )
    } else {
      if (pl.points == F) {
        pl = pl + geom_point(shape = ".", alpha = 1)
      } else {
        pl = pl + geom_point(size = pt.size)
      }
    }
    pl = pl +  scale_color_gradientn(
      colors = .colors,
      na.value = "#DDDDDD",
      limits = c(e.min, e.max),
      breaks = pretty_breaks(3)
    )
    pl = pl +  guides(
      color = guide_colorbar(
        title = .leg.title, title.hjust = 0, barwidth = unit(legend.wh[1],'lines'),
        barheight = unit(legend.wh[2], 'lines'), ticks.linewidth = 1.5/.pt
      )
    ) +
      {if(nchar(names(ftr)) == 0)ggtitle(ftr)} +
      {if(!nchar(names(ftr)) == 0)ggtitle(names(ftr))} +
      mytheme(base_size = base.size) +
      theme(
        aspect.ratio = aspect.ratio,
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.spacing = unit(1, "lines"),
        plot.title = element_text(hjust = 0.5, colour = "black", size = rel(.title.size))
      ) +
      xlab(NULL) + ylab(NULL)

    if (.raster == T) {
      pl = ggrastr::rasterize(pl, layers='Point', dpi=.raster.dpi)
    }

    ftr.pl[[i]] = pl
  }

  if (plot.grid == T) {
    if(length(ftr.pl) <= 4) {ftr.pl = c(ftr.pl, vector(mode = "list", length = (4 - length(ftr.pl))))}
    if(length(ftr.pl) > 4 && length(ftr.pl) <= 8) {ftr.pl = c(ftr.pl, vector(mode = "list", length = (8 - length(ftr.pl))))}
    plot_grid(plotlist = ftr.pl, ncol = 4, scale = .95, align = "vh")
  } else {
    ftr.pl
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Dimension reduction: for phenodata (meta.data)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dimreduc_pheno = function(
  .obj,
  .target = NULL,
  .reduc = "umap",
  .raster.ggrastr = FALSE,
  .raster.scattermore = FALSE,
  .raster.scattermore.pixel = c(512,512),
  .raster.scattermore.pointsize = 0,
  .raster.dpi = 150,
  .sort = F,
  na_cutoff = NULL,
  .col.palette = "2",
  .col.pal.dicrete = ggthemes::tableau_color_pal()(10),
  .col.scico = "roma",
  .col.scico.d = -1,
  .col.scico.b = 0,
  .col.scico.e = 1,
  quantile.fltr = F,
  pl.points = F,
  pt.size = 1,
  leg.title = NULL
) {

  reduc = data.frame(.obj@reductions[[.reduc]]@cell.embeddings)
  colnames(reduc) = c("DIM_1", "DIM_2")
  reduc = cbind(.obj@meta.data, reduc)

  if(is.null(.target)) {
    reduc$tmp = "1"
    .target = "tmp"
    cont = F
  } else {
    cont = is.numeric(reduc[[.target]])
  }

  if(grepl("sne", .reduc, ignore.case = T)) {
    .x.title = "tSNE 1"; .y.title = "tSNE 2"
  } else {
    .x.title = "DIM 1"; .y.title = "DIM 2"
  }
  if(grepl("umap", .reduc, ignore.case = T)) {
    .x.title = "UMAP 1"; .y.title = "UMAP 2"
  } else {
    .x.title = "DIM 1"; .y.title = "DIM 2"
  }

  if(quantile.fltr) {
    qu.max = quantile(reduc[[.target]][ reduc[[.target]] > 0 ], .999)
    reduc[[.target]][reduc[[.target]] > qu.max] = qu.max
  }

  if(.sort == T) {
    reduc = reduc[order(reduc[[.target]], decreasing = F), ]
  } else {
    set.seed(1234)
    reduc = reduc %>% dplyr::sample_frac(1L, replace = FALSE) # permute rows randomly
  }

  if(!is.null(na_cutoff)) {
    e.min = min(reduc[[.target]][reduc[[.target]] > na_cutoff])
    e.max = max(reduc[[.target]])
  } else if (is.null(na_cutoff) & cont) {
    e.min = min(reduc[[.target]])
    e.max = max(reduc[[.target]])
  }

  col.cont = list(
    "1" = rev(MetBrewer::met.brewer("Hokusai1",n=100)),
    "2" = scico::scico(30, palette = .col.scico, direction = .col.scico.d, begin = .col.scico.b, end = .col.scico.e)
  )

  pl =
    ggplot(data = reduc, aes(x = DIM_1, y = DIM_2, col = .data[[.target]]))
    # scattermore::geom_scattermore(pointsize = 5, color="black")+
    # scattermore::geom_scattermore(pointsize = 4, color="white")
    if (.raster.scattermore == T) {
      pl = pl + scattermore::geom_scattermore(
        pointsize = .raster.scattermore.pointsize, pixels = .raster.scattermore.pixel
      )
    } else {
      if (pl.points == F) {
        pl = pl + geom_point(shape = ".", alpha = 1)
      } else {
        pl = pl +
          geom_point(data = reduc[is.na(reduc[[.target]]), ], size = pt.size) +
          geom_point(data = reduc[!is.na(reduc[[.target]]), ], size = pt.size)
          # geom_point(size = pt.size)
      }
    }
    pl = pl + theme(
      aspect.ratio = 1,
      panel.spacing = unit(1, "lines"),
      legend.justification = c(0,.5),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0, face = "bold", colour = "black", size = rel(1))
    ) +
    xlab(.x.title) + ylab(.y.title)
  if(cont) {
    pl = pl +
      scale_color_gradientn(
        colors = col.cont[[.col.palette]],
        na.value = "#DDDDDD",
        limits = c(e.min, e.max),
        breaks = scales::pretty_breaks(3)
      ) +
      guides(
        color = guide_colorbar(
          title.hjust = 0, barwidth = unit(.4, 'lines'), barheight = unit(6, 'lines')
        )
      )
  } else {
    pl = pl +
      scale_color_manual(values = .col.pal.dicrete, na.value = "#DDDDDD") +
      guides(alpha = 'none') +
      guides(colour = guide_legend(
        title = leg.title, ncol = 1,
        override.aes = list(size=3, shape = 16, alpha = 1)
      ))

  }

  if (.raster.ggrastr == T) {
    ggrastr::rasterize(pl, layers='Point', dpi=.raster.dpi)
  } else {
    pl
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA Volcano
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dgea_volcano = function(
    dgea.res = NULL,
    dgea.res.sign = NULL,
    nbr.tops = 7,
    cl.label = NULL, # e.g.  cl.label = setNames(c("CD4", "CD8"), c("CD4 T-Cell", "CD8 T-Cell"))
    sort.by.p = T,
    pval.column = "p_val",
    logFC.column = "avg_log2FC",
    facet.scales = "free_y",
    facet.rows = 1,
    nudge_x = 2,
    x.axis.sym = F,
    x.axis.ext = 0,
    geom.hline = log2(1.5),
    box.padding = 0.3,
    label.padding = .12,
    label.size = 2,
    leg.title = "FDR <0.05",
    cut.y.thresh = NULL
){

  dgea.res$ID = paste0(dgea.res$cluster, "_", dgea.res$feature)
  dgea.res.sign$ID = paste0(dgea.res.sign$cluster, "_", dgea.res.sign$feature)
  dgea.res = dgea.res[dgea.res$cluster %in% unique(dgea.res.sign$cluster), ]

  if(sort.by.p == T) {
    tops_up = dgea.res.sign %>%
      dplyr::filter(.data[[logFC.column]] > 0) %>%
      dplyr::group_by(cluster) %>%
      dplyr::arrange(.data[[pval.column]], -abs(.data[[logFC.column]])) %>%
      dplyr::slice_head(n=nbr.tops)
    tops_up$GENE_CL = paste0(tops_up$feature, "_", tops_up$cluster)

    tops_down = dgea.res.sign %>%
      dplyr::filter(.data[[logFC.column]] < 0) %>%
      dplyr::group_by(cluster) %>%
      dplyr::arrange(.data[[pval.column]], -abs(.data[[logFC.column]])) %>%
      dplyr::slice_head(n=nbr.tops)
    tops_down$GENE_CL = paste0(tops_down$feature, "_", tops_down$cluster)
  } else {
    tops_up = dgea.res.sign %>%
      dplyr::filter(.data[[logFC.column]] > 0) %>%
      dplyr::group_by(cluster) %>%
      dplyr::slice_max(.data[[logFC.column]], n = nbr.tops)
    tops_up$GENE_CL = paste0(tops_up$feature, "_", tops_up$cluster)

    tops_down = dgea.res.sign %>%
      dplyr::filter(.data[[logFC.column]] < 0) %>%
      dplyr::group_by(cluster) %>%
      dplyr::top_n(n = nbr.tops, wt = -.data[[pval.column]])
    tops_down$GENE_CL = paste0(tops_down$feature, "_", tops_down$cluster)
  }

  dgea.res$GENE_CL = paste0(dgea.res$feature, "_", dgea.res$cluster)
  dgea.res = dgea.res %>% mutate(label_up = ifelse(GENE_CL %in% tops_up$GENE_CL, feature, ""))
  dgea.res = dgea.res %>% mutate(label_down = ifelse(GENE_CL %in% tops_down$GENE_CL, feature, ""))
  dgea.res$SIGNIFICANT = dgea.res$ID %in% dgea.res.sign$ID
  dgea.res$SIGNIFICANT = factor(dgea.res$SIGNIFICANT, levels = c(TRUE, FALSE))

  axis.max = max(abs(dgea.res[[logFC.column]]))
  dgea.res = dgea.res %>%
    mutate(p_val = ifelse(p_val == 0, 1e-300, dgea.res[[pval.column]]))
  if(!is.null(cut.y.thresh)){
    dgea.res$p_val[dgea.res$p_val < cut.y.thresh] = cut.y.thresh
  }

  if(is.null(cl.label)) {
    cl.label = setNames(unique(dgea.res$cluster), unique(dgea.res$cluster))
  }

  g <- make_gradient(
    deg = 180, n = 500,
    cols = scico::scico(
      9, palette = 'vik', begin = .3, end = .7, direction = -1,
      )
  )
  set.seed(42)
  pl =
    ggplot(dgea.res, aes(x = .data[[logFC.column]], y = -log10(p_val))) +
    annotation_custom(
      grob = g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    ) +
    geom_point(data = subset(dgea.res, SIGNIFICANT == F), aes(color = SIGNIFICANT), size = .5) +
    geom_point(data = subset(dgea.res, SIGNIFICANT == T), aes(color = SIGNIFICANT), size = .5) +
    # scattermore::geom_scattermore(data = subset(dgea.res, SIGNIFICANT == F), aes(color = SIGNIFICANT), pointsize = 3) +
    # scattermore::geom_scattermore(data = subset(dgea.res, SIGNIFICANT == T), aes(color = SIGNIFICANT), pointsize = 3) +
    theme(
      panel.spacing = unit(1.5, "lines"),
      panel.grid.minor = element_blank(),
      legend.position="bottom",
      strip.text = element_text(size = rel(1), face = "plain")
    ) +
    facet_wrap(
      ~ cluster, scales = facet.scales, labeller = labeller(cluster = cl.label),
      nrow = facet.rows
    ) +
      geom_label_repel(
      data = subset(dgea.res, SIGNIFICANT == T & label_up != ""),
      label = subset(dgea.res, SIGNIFICANT == T & label_up != "")$label_up,
      segment.colour = "black",
      size = label.size,
      direction = "y",
      hjust = .5,
      nudge_x = nudge_x,
      nudge_y = -2,
      segment.size = .2,
      box.padding = box.padding,
      label.padding = label.padding,
      min.segment.length = 0,
      max.overlaps = 50
    ) +
      geom_label_repel(
      data = subset(dgea.res, SIGNIFICANT == T & label_down != ""),
      label = subset(dgea.res, SIGNIFICANT == T & label_down != "")$label_down,
      segment.colour = "black",
      size = label.size,
      direction = "y",
      hjust = .5,
      # xlim = c(NA, 0),
      nudge_x = -nudge_x,
      segment.size = .2,
      box.padding = box.padding,
      label.padding = label.padding,
      min.segment.length = 0,
      max.overlaps = 50
    ) +
    geom_vline(xintercept = -geom.hline, linetype = "dashed", linewidth = .2) +
    geom_vline(xintercept = geom.hline, linetype = "dashed", linewidth = .2) +
    scale_color_manual(values = c("TRUE" = "#555555", "FALSE" = "#BBBBBB")) +
    labs(y = "-Log10(p-value)", x = "Log2 fold change", colour = leg.title)  +
    guides(colour = guide_legend(override.aes = list(size=3)))
  if(x.axis.sym == T) {
    pl = pl + xlim(-axis.max - x.axis.ext, axis.max + x.axis.ext)
  }
  pl
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Bubble plot for cell identity markers
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
mrkr_bubble = function(
    dgea.res = NULL,
    se.obj = NULL,
    ref.group = "celltype",
    effect.size = "avg_log2FC",
    positive.effect.size = T,
    order.by.effect.size = T,
    nbr.tops = 5,
    quantile.cut = .99,
    features = NULL,
    features.rm = NULL,
    features.grep.rm = NULL,
    font.size = 10,
    aspectRatio = 2,
    pt.size = 1,
    rm.rp = T,
    range = c(1, 4),
    export.table = F,
    rename.ftrs = F,
    rename.ftrs.col = "feature_name",
    barwidth = unit(.5, 'lines'),
    barheight = unit(4, 'lines')
) {

  library(scico)
  library(circlize)
  library(ComplexHeatmap)

  DefaultAssay(se.obj) = "RNA"

  dgea.res.all = dgea.res
  dgea.res = dgea.res[dgea.res$significant == T, ]

  if (!is.null(features.rm)) {
    dgea.res = dgea.res[!dgea.res$gene %in% features.rm, ]
  }

  if (!is.null(features.grep.rm)) {
    dgea.res = dgea.res[!grepl(features.grep.rm, dgea.res$gene), ]
  }

  if (rm.rp == T) {
    dgea.res = dgea.res[!grepl("^RPL|^RPS", dgea.res$gene), ]
  }

  if(is.null(features)) {

    dge.res = dgea.res[!is.na(dgea.res[[effect.size]]), ]
    if (positive.effect.size == T) {
      dge.res = dge.res[dge.res[[effect.size]] > 0, ]
    }

    if (order.by.effect.size == T) {
      dge.res = dge.res %>% dplyr::arrange(cluster, desc(abs(!!as.name(effect.size))))
    } else {
      dge.res = dge.res %>% dplyr::arrange(cluster, p_val_adj)
    }

    ct.keep = unique(as.character(dge.res[["cluster"]]))
    se.obj = se.obj[, se.obj@meta.data[[ref.group]] %in% ct.keep]
    se.obj@meta.data = droplevels(se.obj@meta.data)

    tops = dge.res %>%
      dplyr::group_by(cluster) %>%
      dplyr::slice_head(n = nbr.tops) %>% data.frame()

    lvls = levels(se.obj@meta.data[[ref.group]])
    lvls.order = setNames(seq(1, length(lvls)), lvls)
    tops$order = lvls.order[match(tops$cluster, names(lvls.order))]
    tops = tops[naturalsort::naturalorder(tops$order), ]
    export.tops = tops
    tops = tops$gene

  } else {
    tops = features
  }

  se.obj = se.obj[unique(tops), ]
  mat = AverageExpression(se.obj, group.by = ref.group, assays = "RNA", slot = "data")[[1]]
  mat = t(scale(t(mat)))
  df = reshape2::melt(mat)
  colnames(df) = c("GENE", "CT", "AVE")
  df$CT = gsub("^g", "", df$CT)
  # df$CT = gsub("-", "_", df$CT)

  v.max = quantile(df$AVE, quantile.cut)
  df$AVE[abs(df$AVE) > v.max] = v.max

  df$PERC = dgea.res.all$pct.1[match(paste0(df$GENE, "_", df$CT), paste0(dgea.res.all$gene, "_", dgea.res.all$cluster))]
  df$PERC = df$PERC * 100
  df$PERC[df$PERC < 0.01] = NA
  if(!is.null(levels(se.obj@meta.data[[ref.group]]))){
    df$CT = factor(df$CT, levels = levels(se.obj@meta.data[[ref.group]]))
  }

  if(rename.ftrs == T) {
    ftrs = se.obj@assays$RNA@meta.features
    lvls = as.character(ftrs[[rename.ftrs.col]])[match(levels(df$GENE), rownames(ftrs))]
    df$GENE = ftrs[[rename.ftrs.col]][match(df$GENE, rownames(ftrs))]
    df$GENE = factor(df$GENE, levels = lvls)
  }

  pl =
    ggplot(df, aes(x = CT,y = GENE, size = PERC)) +
    geom_point(aes(fill = AVE), shape = 21, stroke = .3) +
    scale_size(range = range) +
    scale_fill_scico(palette = "vik", midpoint = 0, limits = c(-v.max,v.max)) +
    mytheme(base_size = font.size) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(.5, "lines"),
      axis.text.x = element_text(angle=45, hjust=1, vjust = 1),
      axis.text.y = element_text(size = rel(1)),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      # axis.ticks.x = element_blank(),
      aspect.ratio = aspectRatio
    ) +
    guides(
      size = guide_legend(title = "Percent\nExpressed"),
      fill = guide_colorbar(
        title = "Scaled\nAverage\nExpression", order = 1,
        title.hjust = 0, barwidth = barwidth, barheight = barheight
      )
    )

  if(export.table == T) {
    export.tops
  } else {
    pl
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Gene set enrichment for most abundant clonotypes
# Heatmap: Result from UCell; ucell_enrich()
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
enrich_heatmap = function(
    se.w = NULL,
    meta.vars = c("SAMPLE", "TopClones_2"),
    .slot = "custom_UCell_score",
    pl.title = NULL
){

  var.1 = meta.vars[[1]]
  var.2 = meta.vars[[2]]
  se.w@meta.data$VAR_1 = se.w@meta.data[[var.1]]
  se.w@meta.data$VAR_2 = se.w@meta.data[[var.2]]

  DefaultAssay(se.w) = .slot
  df = FetchData(se.w, vars = c(c("VAR_1", "VAR_2"), rownames(se.w)), layer = "data")

  df.ftrs =  df[ , -which(names(df) %in% c("VAR_1", "VAR_2"))]
  df.ftrs[df.ftrs < .1] = 0 # Methods !!!
  df = cbind(df[ , which(names(df) %in% c("VAR_1", "VAR_2"))], df.ftrs)

  df.fltr = reshape2::melt(df, ids = c("VAR_1", "VAR_2"))
  df.fltr = df.fltr %>%
    dplyr::group_by(VAR_1, VAR_2, variable) %>%
    dplyr::mutate(total = n()) %>%
    dplyr::mutate(bool = value > 0) %>%
    dplyr::mutate(bool_sum = sum(bool)) %>%
    dplyr::mutate(frac = bool_sum/total) %>%
    data.frame()
  df.fltr = df.fltr[!duplicated(paste0(df.fltr$VAR_1, df.fltr$VAR_2, df.fltr$variable)), ]
  df.fltr = df.fltr %>%
    dplyr::filter(frac < .25) %>%
    dplyr::group_by(VAR_2, variable) %>%
    dplyr::mutate(total = n())  %>%
    dplyr::filter(total == 2)
  df.fltr$variable = gsub("-", "\\.", df.fltr$variable)

  # l = split(df, df[, "VAR_2"])
  # l = lapply(l, function(x){
  #   mat.scaled = base::scale(x[ , -which(names(x) %in% c("VAR_1", "VAR_2"))])
  #   mat.scaled = data.frame(mat.scaled)
  #   x = cbind(
  #     x[ , which(names(x) %in% c("VAR_1", "VAR_2"))],
  #     mat.scaled
  #   )
  #   x
  # })
  # df = do.call("rbind", l)

  mat = as.data.frame(base::scale(df[ , -which(names(df) %in% c("VAR_1", "VAR_2"))]))
  df = cbind(
    df[ , which(names(df) %in% c("VAR_1", "VAR_2"))],
    mat
  )

  df.m = reshape2::melt(df, ids = c("VAR_1", "VAR_2"))

  wlx.res = df.m %>%
    dplyr::group_by(VAR_2, variable) %>%
    rstatix::wilcox_test(
      data = .,
      as.formula(paste0("value ~ ", "VAR_1"))
    ) %>%
    rstatix::adjust_pvalue(method = "bonferroni") %>%
    rstatix::add_significance("p.adj")
  # wlx.res %>% data.frame()

  df.m = df.m %>%
    dplyr::group_by(VAR_1, VAR_2, variable) %>%
    dplyr::summarise(ave = mean(value))

  # Remove features with negative values accros all samples
  df.m$ave[is.na(df.m$ave)] = 0
  rm.ftrs = df.m %>%
    dplyr::group_by(VAR_2, variable) %>%
    dplyr::filter(ave <= 0) %>%
    dplyr::group_by(variable) %>%
    dplyr::count()
  l = length(unique(df.m[["VAR_1"]])) * length(unique(df.m[["VAR_2"]]))
  df.m = df.m[!df.m$variable %in% rm.ftrs[rm.ftrs$n == l, ]$variable, ]

  # Cluster rows
  df.dc = df.m %>% reshape2::dcast(VAR_1 + VAR_2 ~ variable)
  colnames(df.dc)[1:2] = c("VAR_1", "VAR_2")
  full_dend <- as.dendrogram(
    hclust(dist( t(df.dc[ , -which(names(df.dc) %in% c("VAR_1", "VAR_2"))])  ))
  )
  df.m$variable = factor(df.m$variable, levels = labels(full_dend))

  wlx.res$ID = paste0(wlx.res[["VAR_2"]], "_", wlx.res$variable)
  df.m$ID = paste0(df.m[["VAR_2"]], "_", df.m$variable)
  df.m$p.adj = wlx.res$p.adj[match(df.m$ID, wlx.res$ID)]
  df.m$p.adj.signif = wlx.res$p.adj.signif[match(df.m$ID, wlx.res$ID)]
  # df.m$ave[df.m$p.adj >= 0.05] = NA
  # df.m$ave[paste0(df.m$variable, df.m[["VAR_2"]]) %in% paste0(df.fltr$variable, df.fltr[["VAR_2"]])] = NA
  df.m$p.adj.signif[paste0(df.m$variable, df.m[["VAR_2"]]) %in% paste0(df.fltr$variable, df.fltr[["VAR_2"]])] = "ns"
  r = table(df.fltr$variable)[table(df.fltr$variable) == 4]
  df.m = df.m[!df.m$variable %in% names(r), ]
  r = table(df.m$variable, df.m$p.adj.signif == "ns")
  if(any(colnames(r) == "TRUE")) {
    r = as.data.frame.matrix(r)
    df.m = df.m[!df.m$variable %in% rownames(r[r[, 2] == 4,]), ]
  }
  lbls = setNames(unique(gsub("\\.", " ", df.m$variable)), unique(df.m$variable))

  df.m$p.adj.signif[is.na(df.m$ave)] = "ns"
  df.m$p.adj.signif[df.m$VAR_1 == unique(df.m$VAR_1)[1]] = NA
  df.m$ave[df.m$ave > quantile(df.m$ave, .99)] = quantile(df.m$ave, .99)

  ggplot(df.m) +
    geom_tile(aes(x = VAR_1, y =variable, fill=ave)) +
    geom_text(aes(x = VAR_1, y =variable, label = p.adj.signif), nudge_x = -.5, size = 2) +
    geom_tile(data=df.m[is.na(df.m$ave), ], aes(x=VAR_1, y=variable, col="#BBBBBB"), fill = "#BBBBBB") +
    scale_color_manual(name="FDR<0.05", labels=NULL, values="#BBBBBB") +
    scale_fill_scico(palette = "vik", midpoint = 0, na.value = "#BBBBBB", begin = .1, end = .9) +
    facet_grid(~ VAR_2) +
    scale_y_discrete(labels= lbls) +
    guides(
      fill = guide_colorbar(
        title = "Scaled\nAverage\nEnrichment", order = 1,
        title.hjust = 0, barwidth = unit(.4, 'lines'), barheight = unit(4, 'lines'),
        ticks.linewidth = 1.5/.pt
      )
    ) + xlab(NULL) + ylab(NULL) +
    theme(
      axis.text.x = element_text(angle=45, hjust=1, vjust = 1.05),
      plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1)),
      axis.ticks.y = element_blank(),
      legend.key.size = unit(9,"pt")
    ) +
    ggtitle(pl.title)

}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Make grandient for ggplot (background)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
make_gradient <- function(deg = 45, n = 100, cols = blues9) {

  .cran_packages = c("grid", "ggplot2","RColorBrewer")
  .inst = .cran_packages %in% installed.packages()
  if (any(!.inst)) {
    install.packages(.cran_packages[!.inst], repos = "http://cran.rstudio.com/")
  }
  for (pack in .cran_packages) {
    suppressMessages(library(
      pack,
      quietly = TRUE,
      verbose = FALSE,
      character.only = TRUE
    ))
  }

  cols <- colorRampPalette(cols)(n + 1)
  rad <- deg / (180 / pi)
  mat <- matrix(
    data = rep(seq(0, 1, length.out = n) * cos(rad), n),
    byrow = TRUE,
    ncol = n
  ) +
    matrix(
      data = rep(seq(0, 1, length.out = n) * sin(rad), n),
      byrow = FALSE,
      ncol = n
    )
  mat <- mat - min(mat)
  mat <- mat / max(mat)
  mat <- 1 + mat * n
  mat <- matrix(data = cols[round(mat)], ncol = n)
  grid::rasterGrob(
    image = mat,
    width = unit(1, "npc"),
    height = unit(1, "npc"),
    interpolate = TRUE
  )
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Bubble plot: TOP DE genes
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
de_tops_bubble = function(
    df,
    se.obj,
    pl.title = NULL,
    nbr.tops = 15,
    group = "CT_L1_COARSE",
    ctrst = "SAMPLE",
    ctrst.lvl = NULL
){
  tops.up =
    df %>%
    dplyr::filter(avg_log2FC > 0) %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>%
    slice_head(n = nbr.tops)
  tops.dwn = df %>%
    dplyr::filter(avg_log2FC < 0) %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(p_val_adj, avg_log2FC) %>%
    slice_head(n = nbr.tops)
  tops = rbind(
    tops.dwn[order(tops.dwn$avg_log2FC), ],
    tops.up[order(tops.up$avg_log2FC, decreasing = F), ]
  )

  se.obj$group = se.obj@meta.data[[group]]
  se.obj$ctrst = se.obj@meta.data[[ctrst]]

  mat = GetMatrixFromSeuratByGroupMulti(
    obj= se.obj, features = tops$feature,
    CT_L1_COARSE, ctrst
  )
  exp_mat = reshape2::melt(mat$exp_mat)
  colnames(exp_mat) = c("FTR", "GROUP", "AVE")
  exp_mat$GROUP = gsub(".+\\|", "", exp_mat$GROUP)
  percent_mat = reshape2::melt(mat$percent_mat)
  df = cbind(exp_mat, PERC = percent_mat$value)
  df$FTR = factor(df$FTR, levels = tops$feature)
  g <- make_gradient(
    deg = 180, n = 500,
    cols = scico::scico(
      9, palette = 'vik', begin = .2, end = .8, direction = -1,
    )
  )

  if(!is.null(ctrst.lvl)) {
    df$GROUP = factor(df$GROUP, levels = ctrst.lvl)
  }

  ggplot(df, aes(FTR, GROUP, fill = AVE, size = PERC * 100)) +
    annotation_custom(
      grob = g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    ) +
    geom_point(colour="black", pch=21, stroke = .2) +
    geom_vline(xintercept = nbr.tops + .5, linewidth = .1, linetype = "dashed") +
    scale_size(range = c(1, 3.5), breaks = c(10, 25, 50, floor(max(df$PERC) * 10) * 10)) +
    mytheme() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle=45, hjust=1, vjust = 1, size = rel(.9)),
      legend.text = element_text(margin = margin(l = -5, unit = "pt"), size = rel(1)),
      legend.title  = element_text(margin = margin(r = -2, l = 10, unit = "pt")),
      legend.margin = margin(t = -3),
      plot.title = element_text(hjust = 0.5, face = "plain")
    ) +
    scale_fill_scico(palette = "bilbao", begin =  0, end = 1, direction = -1) +
    guides(
      fill = guide_colorbar(
        title = "Average\nExpression", order = 1,
        title.hjust = 0, title.vjust = .75, barwidth = unit(4, 'lines'),
        barheight = unit(.4, 'lines'), ticks.linewidth = 1.5/.pt
      ),
      size = guide_legend(title = "Percent\nExpressed")
    ) + xlab(NULL) + ylab(NULL) +
    ggtitle(pl.title)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Violin plots for adt genes.
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
plot_marker = function(
    obj = NULL,
    obj.assay = "ADT",
    ftrs = NULL,
    col = NULL,
    col.lvl = NULL,
    label = NULL,
    max.y.val = NULL
){

  obj@meta.data = droplevels(obj@meta.data)

  DefaultAssay(obj) = obj.assay
  expr.data = FetchData(obj, ftrs, slot = "data")
  stopifnot(identical(rownames(expr.data), rownames(obj@meta.data)))
  expr.data = expr.data %>%
    bind_cols(obj@meta.data) %>%
    pivot_longer(cols =!c(colnames(obj@meta.data)) , names_to='FEATURE', values_to='EXPRS') %>%
    data.frame()

  expr.data = expr.data %>%
    filter(.data[[col]] %in% !!col.lvl)

  keep = names(table(expr.data$orig.ident)[table(expr.data$orig.ident) > 1])
  expr.data = expr.data[expr.data$orig.ident %in% keep, ]
  expr.data = droplevels(expr.data)

  newlevels = expr.data %>%
    group_by(orig.ident, RESPONSE_CONSENSUS) %>%
    summarise(AVE_EXPRS = median(EXPRS)) %>%
    arrange(-AVE_EXPRS) %>%
    ungroup %>% select(orig.ident) %>% unlist
  expr.data = expr.data %>% mutate(orig.ident = factor(orig.ident, levels = newlevels))

  s = levels(expr.data$orig.ident)
  p = as.character(expr.data$PATIENT_ID[match(s, expr.data$orig.ident)])
  p = paste0(gsub("Patient 0", "P", p))
  x.labels = setNames(p, s)

  print(max(expr.data$EXPRS))
  if(is.null(max.y.val)) {
    max.y.val = max(expr.data$EXPRS)
  }

  pl =
    ggplot(expr.data, aes(x=orig.ident, y=EXPRS,fill=RESPONSE_CONSENSUS))+
    # geom_violin(linewidth = 0.2) +
    geom_violin(linewidth = 0.2, scale = "width", width = .5) +
    geom_jitter(shape = ".", alpha = 0.5, height=0, width = .1, show.legend = F, color='#555555')+
    scale_x_discrete(labels = x.labels) +
    stat_summary(fun= "median",geom = "crossbar", width = 0.4, show.legend = F, lwd = .3)+
    scale_fill_manual(values = c("CR" = "#6699CC", "nonCR" = "#997700")) +
    mytheme(base_size = 5) +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
      axis.title.x = element_blank()
    ) +
    annotate("text",  x=Inf, y = Inf, label = label, vjust=1.6, hjust=1.4, size = 2) +
    ylab("Expression") + labs(fill = NULL) +
    ylim(0, max.y.val)

  pl
}


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Plot for complete copy number profile
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
plotProfile_mark<-function (pos,chr,gains,losses, chrominfo = human.chrom.info.May04,
                            ylim = c(-1,1), ylb = "Fraction gained or lost", centromere = FALSE, chromosomes = NULL, markchr,markpos,main="",...)
{
  chromosomes <- unique(chr)
  chrom <- chr
  mb <- pos/1e+06
  chrominfo <- chrominfo[chromosomes, ]
  chrom.start <- c(0, cumsum(chrominfo$length)[1:(length(chromosomes) -
                                                    1)])/1000
  chrom.centr <- chrom.start + chrominfo$centr/1000
  chrom.mid <- chrom.start + chrominfo$length/2000
  op <- par(cex = 0.6, pch = 18, lab = c(1, 6, 7), cex.axis = 2,
            xaxs = "i", xaxt = "n", mai = c(0.4, 0.6, 0.8, 0.3),mar=c(2,3,3.5,1),
            cex.lab = 2, cex.main = 2)
  genomepos <- rep(0, length(mb))
  for (i in 1:length(chromosomes)) genomepos[chrom == chromosomes[i]] <- mb[chrom ==
                                                                              chromosomes[i]] + chrom.start[i]
  data <- gains
  plot(genomepos, data, ylim = ylim, xlab = "", ylab = ylb,
       xlim = c(0, sum(chrominfo$length)/1000), col = "darkred", type = "h",lwd=3,
       ...)
  title(main=main,adj=0,line=2)
  data <- (-losses)
  lines(genomepos, data, type = "h", col = "darkblue",lwd=3)
  for (i in 1:ceiling(length(chromosomes)/2)) mtext(paste("",
                                                          chromosomes[(i - 1) * 2 + 1]), side = 1, at = (chrom.mid[(i -
                                                                                                                      1) * 2 + 1]), line = 0.8, col = "black", cex = 1.0)
  if (length(chromosomes) > 1)
    for (i in 1:floor(length(chromosomes)/2)) mtext(paste("",
                                                          chromosomes[i * 2]), side = 3, at = chrom.mid[i *
                                                                                                          2], line = 0.4, col = "black", cex = 1.0)
  abline(v = c(chrom.start, (chrom.start[nrow(chrominfo)] +
                               chrominfo$length[nrow(chrominfo)])), lty = 1)
  abline(h = seq(-1, 1, b = 0.2), lty = 2, col = "grey")
  abline(h = 0, col = 2)
  if (centromere == TRUE) {
    abline(v = (chrom.centr), lty = 3, col = 4)
  }
  par(op)
}
