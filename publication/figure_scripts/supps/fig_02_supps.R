.cran_packages = c(
  "yaml", "ggplot2","reshape2", "dplyr", "naturalsort", "devtools", "scales",
  "stringr", "Seurat", "tibble", "tidyr", "HGNChelper", "forcats", "cowplot",
  "rlang", "remotes", "scGate", "patchwork", "openxlsx", "scCustomize", "ggpubr",
  "tidyverse", "scCustomize", "ggh4x", "ggrepel", "anndata", "scico", "ggVennDiagram",
  "semla", "ggfittext"
)
.bioc_packages = c("dittoSeq", "scRepertoire", "UCell","aCGH","infercnv")

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
  remotes::install_github('satijalab/seurat-wrappers')
}
library(SeuratWrappers)

if (any(!"SeuratData" %in% installed.packages())) {
  devtools::install_github('satijalab/seurat-data')
}
library(SeuratData)

if (any(!"ProjecTILs" %in% installed.packages())) {
  remotes::install_github("carmonalab/ProjecTILs")
}
library(ProjecTILs)

if (any(!"SeuratDisk" %in% installed.packages())) {
  remotes::install_github("mojaveazure/seurat-disk")
}
library(SeuratDisk)

source("code/helper/styles.R")
source("code/helper/functions_plots.R")
source("code/helper/functions.R")
source("code/helper/ora.R")

theme_set(mytheme(base_size = 8))

theme.custom = mytheme(base_size = 8) &
  theme(
    legend.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.title = element_blank(),
    panel.border = element_blank(),
    plot.title = element_text(hjust = .5, size = rel(1))
  )

base.size = 8

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA | single-cell
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")

se.meta.t = readRDS(paste0(manifest$multi_omics$work, "integration/seurat_harmony_t.Rds"))
se.meta.cd4 = readRDS(paste0(manifest$multi_omics$work, "integration/seurat_harmony_cd4.Rds"))
se.meta.cd8 = readRDS(paste0(manifest$multi_omics$work, "integration/seurat_harmony_cd8.Rds"))

se.meta.cd4$celltype = as.factor(gsub("_EOMES", "", se.meta.cd4$celltype))

sample.lbls = setNames(c("AphNB", "P2248", "P2249", "P2264"), c("Aph", "BM", "PB", "PB+Dexa"))
se.meta.t$SAMPLE = names(sample.lbls)[match(se.meta.t$orig.ident, sample.lbls)]

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

se.meta.t@meta.data = se.meta.t@meta.data %>% mutate(
  TopClones_2 = case_when(
    TopClones == "Cl_1" | TopClones == "Cl_3" ~ "Cl_1_3",
    TRUE ~ TopClones
  )
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA | spatial
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sp.work = paste0(manifest$spatial$cellranger, "C24_13076_one_probe/outs/") # Haut Köln_normaler Sondenpool
sp.k = get_special_se(sp.work)
sp.work = paste0(manifest$spatial$cellranger, "H21088-21_Sonden_SpikeIn_one_probe/outs/") # 21088/21_spike-in Carvyctysonden
sp.m = get_special_se(sp.work)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# T-cells: DimReduc by samples
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pd = get_metadata(se.meta.t)
t.ct.pl =
  ggplot(data = pd, aes(x = umap_1, y = umap_2, col = celltype)) +
  scattermore::geom_scattermore(pointsize = 5, color="black")+
  scattermore::geom_scattermore(pointsize = 4, color="white") +
  geom_point(size = .01) +
  # geom_point(shape = ".") +
  theme(
    legend.position = "bottom",
    legend.margin = margin(t=-5),
    legend.title = element_text(margin = margin(r = 3)),
    legend.spacing.y = unit(2, 'mm'),
    legend.key.size = unit(4, "mm"),
    legend.text = element_text(margin = margin(l = -2.5, unit = "pt")),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    plot.title = element_blank()
  ) +
  guides(colour = guide_legend(
    title = NULL, ncol = 5, override.aes = list(shape = 16, size = 3)
  )) +
  scale_color_manual(values = til.col, na.value = "red")


ggsave2(
  filename="publication/figures_tables/supps/fig_02_t_reduc_by_sample.png",
  t.ct.pl +
    facet_wrap(~ orig.ident, nrow = 1) +
    theme(
      panel.spacing = unit(1, "lines"),
      legend.position = "bottom"
    ),
  width = 160,
  height = 50,
  dpi = 200,
  bg = "white",
  units = "mm",
  scale = 1.6
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Marker for CD4 and CD8 cells
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.work.cd8 = se.meta.cd8
table(se.work.cd8$celltype)
res.wlx.cd8 = wlx_test_dor(
  obj = se.work.cd8,
  se.grp = "celltype",
  downsample = F,
  # downsample.nbr = 200,
  only.pos = F
)

res.wlx.cd8$significant = res.wlx.cd8$p_val_adj < 0.05 & res.wlx.cd8$avg_log2FC > log2(1.5) &
  res.wlx.cd8$logDOR > 0  & res.wlx.cd8$pct.1 > 0.25
table(res.wlx.cd8[res.wlx.cd8$significant == T, ]$cluster)

mrkr.cd8 =
  mrkr_bubble(
    dgea.res = res.wlx.cd8,
    se.obj = se.work.cd8,
    ref.group = "celltype",
    nbr.tops = 8,
    order.by.effect.size = F,
    positive.effect.size = T,
    range = c(1, 5),
    font.size = base.size,
    features.grep.rm = "^TRA|^TRB|IGL|MT-",
    aspectRatio = NULL,
    barwidth = unit(.5, 'lines'),
    barheight = unit(3, 'lines')
  ) + coord_flip()

###

se.work.cd4 = se.meta.cd4
table(se.work.cd4$celltype)
res.wlx.cd4 = wlx_test_dor(
  obj = se.work.cd4,
  se.grp = "celltype",
  downsample = F,
  # downsample.nbr = 100,
  only.pos = F
)

res.wlx.cd4$significant = res.wlx.cd4$p_val_adj < 0.05 & res.wlx.cd4$avg_log2FC > log2(1.5) &
  res.wlx.cd4$logDOR > 0 & res.wlx.cd4$pct.1 > 0.25
table(res.wlx.cd4[res.wlx.cd4$significant == T, ]$cluster)

mrkr.cd4 =
  mrkr_bubble(
    dgea.res = res.wlx.cd4,
    # dgea.res = res.wlx.cd4,
    se.obj = se.work.cd4,
    ref.group = "celltype",
    nbr.tops = 10,
    order.by.effect.size = F,
    positive.effect.size = T,
    range = c(1, 5),
    font.size = base.size,
    features.grep.rm = "^TRA|^TRB|IGL|MT-",
    aspectRatio = NULL,
    barwidth = unit(.5, 'lines'),
    barheight = unit(3, 'lines')
  ) + coord_flip()

ggsave2(
  filename="publication/figures_tables/supps/fig_02_t_marker.png",
  plot_grid(
    NULL,
    mrkr.cd4,
    NULL,
    mrkr.cd8, align = "vh",
    nrow = 4, rel_heights = c(.25, 1, .35, 1),
    labels = c("A)", "", "B)"),
    label_fontface = "bold", label_size = 10
  ),
  width = 160, height = 80, dpi = 200, bg = "white", units = "mm", scale = 1.6
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# CD8 DimReduc + CD45RO/RA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
DefaultAssay(se.meta.cd8) = "ADT"
a = FeaturePlot_scCustom(
  se.meta.cd8, #split.by = "SAMPLE",
  reduction = "umap",
  features = c("CD45RO", "CD45RA"),
  pt.size = .5, na_cutoff = 1,  # max.cutoff = "q99",
  colors_use = scico::scico(30, palette = "navia", direction = -1),
  raster = F
) & mytheme() & theme(aspect.ratio = 1) &
  guides(
    color = guide_colorbar(
      title = NULL, order = 1, barwidth = unit(.4, 'lines'),
      barheight = unit(4, 'lines'), ticks.linewidth = 1.5/.pt
    ),
    size = guide_legend(title = "Percent\nExpressed")
  )
DefaultAssay(se.meta.cd8) = "RNA"
b = DimPlot_scCustom(se.meta.cd8, group.by = "celltype", pt.size = .5, colors_use = til.col) +
  mytheme() + theme(aspect.ratio = 1) + ggtitle("CD8 T")

ggsave2(
  filename="publication/figures_tables/supps/fig_02_CD45RO_RA_marker.png",
  plot_grid(
    b, a, rel_widths = c(1, 1.6), align = "vh",
    labels = c("A)", "B)"),  label_fontface = "bold", label_size = 10
  ),
  width = 160, height = 50, dpi = 200, bg = "white", units = "mm", scale = 1.6
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA: CAR+ vs. CAR-
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.work = se.meta.t
se.work = subset(se.work, orig.ident != "AphNB")

fc.threshold = 1.5
res.wlx = run_wilx(
  obj = se.work, target = "SAMPLE",
  min.pct.thres =  0.25, fc.thresh = fc.threshold,
  contrast.group = "CAR_BY_EXPRS", contrast = c("TRUE", "FALSE")
)
res.wlx.sign = subset(res.wlx, significant == T)

vol.dgea.car =
  dgea_volcano(
    dgea.res = res.wlx,
    dgea.res.sign = res.wlx.sign,
    sort.by.p =  T,
    nbr.tops = 5,
    x.axis.sym = T,
    x.axis.ext = 4,
    nudge_x = 4,
    geom.hline = log2(fc.threshold),
    cl.label = setNames(
      paste0(names(table(res.wlx.sign$cluster)), ": CAR+ vs. CAR-"),
      names(table(res.wlx.sign$cluster))
    ),
    leg.title = paste0("FDR <0.05 & abs(Fold-change) >", fc.threshold),
    cut.y.thresh = 1e-50
  ) & theme(aspect.ratio = 1)

ggsave2(
  filename="publication/figures_tables/supps/fig_02_volcano_dgea_car.png",
  vol.dgea.car,
  width = 160, height = 70, dpi = 200, bg = "white", units = "mm", scale = 1
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA: Clone 1 vs 3
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.work = se.meta.t
se.work = subset(se.work, orig.ident != "AphNB")

fc.threshold = 1.5
res.wlx = run_wilx(
  obj = se.work, target = "SAMPLE",
  min.pct.thres =  0.25, fc.thresh = fc.threshold,
  contrast.group = "TopClones", contrast = c("Cl_3", "Cl_1")
)
res.wlx.sign = subset(res.wlx, significant == T)
table(res.wlx.sign$cluster)

vol.dgea.car =
  dgea_volcano(
    dgea.res = res.wlx,
    dgea.res.sign = res.wlx.sign,
    sort.by.p =  T,
    nbr.tops = 10,
    x.axis.sym = T,
    x.axis.ext = 4,
    nudge_x = 4,
    geom.hline = log2(fc.threshold),
    cl.label = setNames(
      paste0(names(table(res.wlx.sign$cluster)), ": Clone 3 vs. Clone 1"),
      names(table(res.wlx.sign$cluster))
    ),
    leg.title = paste0("FDR <0.05 & abs(Fold-change) >", fc.threshold),
    cut.y.thresh = 1e-50
  ) & theme(aspect.ratio = 1)

ggsave2(
  filename="publication/figures_tables/supps/fig_02_volcano_dgea_clone_3_1.png",
  vol.dgea.car,
  width = 160, height = 70, dpi = 200, bg = "white", units = "mm", scale = 1
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA: BM vs. PB
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.work = se.meta.t
fc.threshold = 1.5
res.wlx = run_wilx(
  obj = se.work, target = "CT_L1_COARSE",
  min.pct.thres =  0.25, fc.thresh = fc.threshold,
  contrast.group = "SAMPLE", contrast = c("BM", "PB")
)
res.wlx.sign = subset(res.wlx, significant == T)
table(res.wlx.sign$cluster)

tbl = as.data.frame(table(res.wlx.sign$avg_log2FC > 0))
if(nrow(tbl) == 1) {
  if(tbl$Var1 == "TRUE"){
    tbl = rbind(tbl, data.frame(Var1 = "FALSE", Freq = 0))
  } else {
    tbl = rbind(data.frame(Var1 = "TRUE", Freq = 0), tbl)
  }
}

lbls = paste0(
  "BM vs. PB (up: ",
  tbl[1, ][[2]],
  ", down: ",
  tbl[2, ][[2]],
  ")"
)

vol.dgea.bmpb =
  dgea_volcano(
  dgea.res = res.wlx,
  dgea.res.sign = res.wlx.sign,
  sort.by.p =  T,
  nbr.tops = 10,
  x.axis.sym = T,
  x.axis.ext = 2,
  geom.hline = log2(fc.threshold),
  cl.label = setNames(lbls, names(table(res.wlx.sign$cluster))),
  leg.title = paste0("FDR <0.05 & abs(Fold-change) >", fc.threshold),
) & theme(aspect.ratio = 1)

ggsave2(
  filename="publication/figures_tables/supps/fig_02_volcano_dgea_bm_pb.png",
  vol.dgea.bmpb,
  width = 180, height = 70, dpi = 200, bg = "white", units = "mm", scale = 1
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Gene set enrichment for most abundant clonotypes
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
mir = readRDS("data/metadata/signatures/mir_gene_set_collection.Rds")
# unique(mir$TermID)

cust.gene.sets = mir[grepl("CD8", mir$TermID), ]
cust.gene.sets = cust.gene.sets[!grepl("Cytopus", cust.gene.sets$TermID), ]
cust.gene.sets$TermID = gsub("CD8-", "", cust.gene.sets$TermID)
cust.gene.sets$TermID = gsub(" .+", "", cust.gene.sets$TermID)
cust.gene.sets = split(cust.gene.sets, cust.gene.sets$TermID)
cust.gene.sets = lapply(cust.gene.sets, function(x){x$GeneID})

se.cl = subset(se.meta.t, TopClones_2 != "Cl_Other")
se.cl = subset(se.cl, orig.ident != "P2248" & orig.ident != "AphNB")
se.cl = ucell_enrich(se.w = se.cl, cust.gene.sets)

# DefaultAssay(se.cl) = "custom_UCell_score"
# g.s.redu = dimreduc_features(
#   se.cl, rownames(se.cl), .reduc = "umap", .quantile.fltr = T, .title.size = 1, pl.points = T, pt.size = .2,
#   .assay = "custom_UCell_score", .na_cutoff = 0.1, plot.grid = F, #legend.size = 6,
#   .x.title = NULL, .y.title = NULL, base.size = 8, legend.wh = c(.6, 5)
# )
# DefaultAssay(se.cl) = "RNA"
# plot_grid(plotlist = g.s.redu, ncol = 3, scale = .98)
# DimPlot_scCustom(se.cl, group.by = "TopClones_2")
# DimPlot_scCustom(se.cl, group.by = "orig.ident")

# Alle Scores <.1 werden auf 0 gesetzt -> Methods !!!
enrich.cl.hm = enrich_heatmap(se.w = se.cl, meta.vars = c("TopClones_2", "SAMPLE"))
ggsave2(
  filename="publication/figures_tables/supps/fig_02_ucell.png",
  plot_grid(NULL, enrich.cl.hm, NULL, nrow = 1, rel_widths = c(.4, 1, .4)),
  width = 180, height = 70, dpi = 200, bg = "white", units = "mm", scale = 1
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA: PB+Dexa vs. PB
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.work = se.meta.t
se.work = subset(se.work, orig.ident != "P2248" & orig.ident != "AphNB")
se.work = subset(se.work, TopClones_2 != "Cl_Other")
se.work$TopClones_2 = gsub("Cl_", "Clone ", se.work$TopClones_2)
fc.threshold = 1.5
res.wlx = run_wilx(
  obj = se.work, target = "TopClones_2",
  min.pct.thres =  0.25, min.cells = 20, fc.thresh = fc.threshold,
  contrast.group = "SAMPLE", contrast = c("PB+Dexa", "PB")
)
res.wlx.sign = subset(res.wlx, significant == T)
table(res.wlx.sign$cluster, res.wlx.sign$avg_log2FC > 0)

# de_tops_bubble(
#   df = res.wlx.sign[res.wlx.sign$cluster == "Cl_1_3", ],
#   se.obj = subset(se.work, TopClones_2 == "Cl_1_3"),
#   pl.title = "Clone 1_3: PB+Dexa vs. PB"
# )
#
# de_tops_bubble(
#   df = res.wlx.sign[res.wlx.sign$cluster == "Cl_2", ],
#   se.obj = subset(se.work, TopClones_2 == "Cl_2"),
#   pl.title = "Clone 2: PB+Dexa vs. PB"
# )
vol.pl =
  dgea_volcano(
  dgea.res = res.wlx,
  dgea.res.sign = res.wlx.sign,
  sort.by.p =  T,
  nbr.tops = 10,
  x.axis.sym = T,
  x.axis.ext = 2,
  label.size = 1.75,
  label.padding = .1,
  cl.label = setNames(
    c("Clone 1_3: PB+Dexa vs. PB\nup: 68, down: 554", "Clone 2: PB+Dexa vs. PB\nup: 47, down: 694"),
    c("Clone 1_3", "Clone 2")
  ),
  leg.title = paste0("FDR <0.05 & abs(Fold-change) >", fc.threshold),
) &
  theme(
    legend.title = element_text(margin = margin(r = 15))
  )

dendro.pl = ggVennDiagram(
  split(f = res.wlx.sign$cluster, x = res.wlx.sign$feature),
  label_size = 3, set_size = 3, edge_size = .5
) +
  theme(legend.position = "none") +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + coord_flip()

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ORA: PB+Dexa vs. PB
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ftrs.l = res.wlx.sign
ftrs.l = split(ftrs.l, ftrs.l$cluster)
ftrs.all.l = lapply(ftrs.l, function(x){
  ftrs = x$feature
  names(ftrs) = x$avg_log2FC
  ftrs
})
lengths(ftrs.all.l)

ora.go.all = parallel::mclapply(ftrs.all.l, function(x){
  run_nmf_ora(
    genes = x,
    universe = rownames(se.work),
    # category = "C2", db.sub = "REACTOME"
    category = "C5", subcategory = "BP"
  )
}, mc.cores = length(ftrs.all.l))

ora.clono.all =
  ora_bubble(
    gsea.res = ora.go.all,
    ftrs.list = ftrs.all.l,
    nbr.tops = 15,
    font.size = 8,
    min.genes = 5,
    max.value = 1.5,
    term.length = 70,
    sort.by.padj = F,
    barwidth = unit(.4, 'lines'),
    barheight = unit(4, 'lines'),
    y.text.size = rel(.8)
  ) +
  ggtitle("PB+Dexa vs. PB")


ggsave2(
  filename="publication/figures_tables/supps/fig_02_dgea_pbdexa_pb.png",
  plot_grid(
    plot_grid(
      vol.pl, NULL,
      plot_grid(NULL, dendro.pl, NULL, nrow = 3, rel_heights = c(.2, 1, .4)),
      ncol = 3, rel_widths = c(1, .1, .5),
      labels = c("A)", "", "B)"),
      label_fontface = "bold", label_size = 10
    ),
    NULL,
    plot_grid(
      NULL,
      ora.clono.all,
      NULL,
      labels = c("C)", "", ""), ncol = 3, rel_widths = c(.1, 1, .15),
      label_fontface = "bold", label_size = 10
    ),
    nrow = 3, rel_heights = c(1.3, .1, 2)
  ),
  width = 180, height = 185, dpi = 200, bg = "white", units = "mm", scale = 1
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Spatial
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
theme.custom = mytheme(base_size = 8) &
  theme(
    axis.text = element_blank(),
    legend.position = "none",
    axis.ticks = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.title = element_blank(),
    panel.border = element_blank(),
    plot.subtitle = element_blank(),
    # plot.subtitle = element_text(hjust = 0.5, size = rel(1), face = "plain"),
    plot.title = element_blank(),
    legend.title = element_blank(),
  )

ftrs =  c("CD3E", "CD4", "CD8B", "ciltacel")

global.max = max(log10(unlist(FetchData(sp.k, ftrs, layer = "counts")) + 1))
ftr.l = spatial_features(
  .obj = sp.k, features = ftrs, pt.size = .1, min.max = c(0, global.max),
  legend.wh = c(4, .4), base.size = 8, slot = "counts", log.counts = T,
  .colors = scico::scico(30, palette = "glasgow", direction = -1), plot.grid = F,
  aspect.ratio = 1.1
)
leg.k = get_legend(
  ftr.l[[1]] + theme(legend.position = "bottom") + guides(
    color = guide_colorbar(
      title = "log10(counts + 1)", title.vjust = 1.05, barwidth = 4,
      barheight = .4, ticks.linewidth = 1.5/.pt
    )
  ) + theme(
    legend.title = element_text(margin = margin(r = 10, unit = "pt")),
  )
)
ftr.l = lapply(ftr.l, function(x){
  x = x + theme(legend.position='none')
})
ftr.l.k = ftr.l

global.max = max(log10(unlist(FetchData(sp.m, ftrs, layer = "counts")) + 1))
ftr.l = spatial_features(
  .obj = sp.m, features = ftrs, pt.size = .1, min.max = c(0, global.max),
  legend.wh = c(4, .4), base.size = 8, slot = "counts", log.counts = T,
  .colors = scico::scico(30, palette = "glasgow", direction = -1), plot.grid = F,
  aspect.ratio = .8
)
leg.m = get_legend(
  ftr.l[[1]] + theme(legend.position = "bottom") + guides(
    color = guide_colorbar(
      title = "log10(counts + 1)", title.vjust = 1.05, barwidth = 4,
      barheight = .4, ticks.linewidth = 1.5/.pt
    )
  ) + theme(
    legend.title = element_text(margin = margin(r = 10, unit = "pt")),
  )
)
ftr.l = lapply(ftr.l, function(x){
  x = x + theme(legend.position='none')
})
ftr.l.m = ftr.l

ggsave2(
  filename="publication/figures_tables/supps/fig_02_spatial.png",
  plot_grid(
    plot_grid(
      MapFeatures(
        sp.k, features = "nFeature_Spatial",pt_alpha = 0,
        image_use = "raw", override_plot_dims = TRUE
      ) & theme.custom & theme(plot.margin = margin(b = 15, t = 10, unit = "pt")),
      plot_grid(
        plot_grid(plotlist = ftr.l.k, nrow = 1),
        leg.k, nrow = 2, rel_heights = c(1, .2)
      ),
      nrow = 1, rel_widths = c(1, 4.5)
    ),
    NULL,
    plot_grid(
      MapFeatures(
        sp.m, features = "nFeature_Spatial",pt_alpha = 0,
        image_use = "raw", override_plot_dims = TRUE
      ) & theme.custom & theme(plot.margin = margin(b = 15, t = 10, r = 1, unit = "pt")),
      plot_grid(
        plot_grid(plotlist = ftr.l.m, nrow = 1),
        leg.m, nrow = 2, rel_heights = c(1, .2)
      ),
      nrow = 1, rel_widths = c(1, 4.3)
    ),
    nrow = 3, rel_heights = c(1, .0, 1),
    labels = c("A)", "", "B)"), label_fontface = "bold", label_size = 10, label_x = -.01
  ),
  width = 180, height = 105, dpi = 200, bg = "white", units = "mm", scale = 1
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# InferCNV copy number profile plots
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
thresh = 0.05

# Clone 1_3
infercnv_c1_3 <- readRDS("./analysis/multi_omics/inferCNV/infercnv_obj_Cl_1_3.rds")

exprs.complete = infercnv_c1_3@expr.data
exprs.complete[exprs.complete>(1+thresh)] <- 2
exprs.complete[exprs.complete<(1-thresh)] <- 0
exprs.complete[exprs.complete>=(1-thresh) & exprs.complete<=(1+thresh)] <- 1
exprs.complete = exprs.complete - 1
gain_perc_1_3 = apply(exprs.complete,1,function(x){sum(x==1)/length(x)})
loss_perc_1_3 = apply(exprs.complete,1,function(x){sum(x==-1)/length(x)})

# Clone 2
infercnv_c2 <- readRDS("./analysis/multi_omics/inferCNV/infercnv_obj_Cl_2.rds")
exprs.complete = infercnv_c2@expr.data
exprs.complete[exprs.complete>(1+thresh)] <- 2
exprs.complete[exprs.complete<(1-thresh)] <- 0
exprs.complete[exprs.complete>=(1-thresh) & exprs.complete<=(1+thresh)] <- 1
exprs.complete = exprs.complete - 1
gain_perc_2 = apply(exprs.complete,1,function(x){sum(x==1)/length(x)})
loss_perc_2 = apply(exprs.complete,1,function(x){sum(x==-1)/length(x)})

# Other Clones
infercnv_Cl_Other <- readRDS("./analysis/multi_omics/inferCNV/infercnv_obj_Cl_Other.rds")
exprs.complete = infercnv_Cl_Other@expr.data
exprs.complete[exprs.complete>(1+thresh)] <- 2
exprs.complete[exprs.complete<(1-thresh)] <- 0
exprs.complete[exprs.complete>=(1-thresh) & exprs.complete<=(1+thresh)] <- 1
exprs.complete = exprs.complete - 1
gain_perc_o = apply(exprs.complete,1,function(x){sum(x==1)/length(x)})
loss_perc_o = apply(exprs.complete,1,function(x){sum(x==-1)/length(x)})

png("publication/figures_tables/supps/fig_02_copy_number_profile.png", width = 280, height = 180, res = 400, bg = "white", units = "mm", type = c("cairo-png"))
  par(mfrow=c(3,1))
  plotProfile_mark(infercnv_c1_3@gene_order$start,as.numeric(gsub(pattern = "chr","",infercnv_c1_3@gene_order$chr)),gain_perc_1_3,loss_perc_1_3,main="Clone 1_3",ylb="")
  plotProfile_mark(infercnv_c2@gene_order$start,as.numeric(gsub(pattern = "chr","",infercnv_c2@gene_order$chr)),gain_perc_2,loss_perc_2,main="Clone 2",ylb="")
  plotProfile_mark(infercnv_Cl_Other@gene_order$start,as.numeric(gsub(pattern = "chr","",infercnv_Cl_Other@gene_order$chr)),gain_perc_o,loss_perc_o,main="Other clones",ylb="")
dev.off()

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# InferCNV copy number heatmaps
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
plot_cnv(infercnv_c1_3,title = "Clone 1_3",png_res=400,cluster_references=F,output_filename="publication/figures_tables/supps/fig_02_copy_number_heatmap_A",output_format="png",obs_title = "Clone 1_3 (cells)")
plot_cnv(infercnv_c2,title = "Clone 2",png_res=400,cluster_references=F,output_filename="publication/figures_tables/supps/fig_02_copy_number_heatmap_B",output_format="png",obs_title = "Clone 2 (cells)")
plot_cnv(infercnv_Cl_Other,title = "Other clones",png_res=400,cluster_references=F,output_filename="publication/figures_tables/supps/fig_02_copy_number_heatmap_C",output_format="png",obs_title = "Other clones (cells)")




