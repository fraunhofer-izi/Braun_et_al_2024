.cran_packages = c(
  "yaml", "ggplot2","reshape2", "dplyr", "naturalsort", "devtools", "scales",
  "stringr", "Seurat", "tibble", "tidyr", "HGNChelper", "forcats", "cowplot",
  "rlang", "remotes", "scGate", "patchwork", "openxlsx", "scCustomize", "ggpubr",
  "tidyverse", "scCustomize", "ggh4x", "ggrepel", "anndata", "scico", "semla",
  "ggalluvial", "ggfittext"
)
.bioc_packages = c("dittoSeq", "scRepertoire", "UCell","GenVisR")

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

if (any(!"semla" %in% installed.packages())) {
  Sys.unsetenv("GITHUB_PAT")
  remotes::install_github("ludvigla/semla")
}
library(semla)

source("code/helper/styles.R")
source("code/helper/functions_plots.R")
source("code/helper/functions.R")
source("code/helper/ora.R")

theme_set(mytheme(base_size = 8))
base.size = 8

leg.text.l = -2.5

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")

se.meta = readRDS(paste0(manifest$multi_omics$work, "integration/seurat_harmony.Rds"))
se.meta.t = readRDS(paste0(manifest$multi_omics$work, "integration/seurat_harmony_t.Rds"))

sample.lbls = setNames(c("AphNB", "P2248", "P2249", "P2264"), c("Aph", "BM-preDexa", "PB-preDexa", "PB-postDexa"))
se.meta$SAMPLE = names(sample.lbls)[match(se.meta$orig.ident, sample.lbls)]
se.meta$SAMPLE = factor(se.meta$SAMPLE, levels = names(sample.lbls))
se.meta.t$SAMPLE = names(sample.lbls)[match(se.meta.t$orig.ident, sample.lbls)]
se.meta.t$SAMPLE = factor(se.meta.t$SAMPLE, levels = names(sample.lbls))

keep.ct = table(se.meta$celltype_short_3) > 50
se.meta = se.meta[, se.meta$celltype_short_3 %in% names(keep.ct[keep.ct])]
se.meta@meta.data = droplevels(se.meta@meta.data)
se.meta$celltype_short_3 = as.character(se.meta$celltype_short_3)
se.meta$celltype_short_3[se.meta$celltype_short_3 == "Progenitor"] = "HSPC"
se.meta$celltype_short_3 = as.factor(se.meta$celltype_short_3)

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
clone.col = setNames(c("#555555", "#004488", "#DDAA33", "#5AB2FF"), clone.col)

se.meta.t@meta.data = se.meta.t@meta.data %>% mutate(
  TopClones_2 = case_when(
    TopClones == "Cl_1" | TopClones == "Cl_3" ~ "Cl_1_3",
    TRUE ~ TopClones
  )
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DimReduc and composition: all celltypes
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pd = get_metadata(se.meta)

all.ct.pl  =
  ggplot(data = pd, aes(x = umap_1, y = umap_2, col = celltype_short_3)) +
  scattermore::geom_scattermore(pointsize = 7, color="black")+
  scattermore::geom_scattermore(pointsize = 6, color="white") +
  scattermore::geom_scattermore(pointsize = 2) +
  theme(
    legend.position = "none",
    legend.margin = margin(t=-5),
    legend.title = element_text(margin = margin(r = 3)),
    legend.spacing.y = unit(2, 'mm'),
    legend.key.size = unit(4, "mm"),
    legend.text = element_text(margin = margin(l = leg.text.l, unit = "pt")),
    axis.title = element_text(size = 6),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    plot.title = element_blank()
  ) +
  guides(colour = guide_legend(
    title = NULL, ncol = 3, override.aes = list(shape = 16, size = 2.5)
  )) +
  scale_color_manual(values = ct.col, na.value = "red") +
  xlab("UMAP 1") + ylab("UMAP 2")

df = dittoBarPlot(
  se.meta, "celltype_short_3", group.by = "SAMPLE",
  data.out = T
)[[2]]

df$grouping = factor(df$grouping, levels = names(sample.lbls))
comp.pl =
  ggplot(df,  aes(grouping, count, fill = label)) +
  geom_col(aes(fill = label),  width = .9) +
  # geom_col(aes(fill = label), position = "fill", width = .9) +
  scale_fill_manual(values = ct.col) +
  theme(
    legend.position = "right",
    axis.title.x = element_blank(),
    legend.key.size = unit(9,"pt"),
    legend.margin = margin(l = 1),
    strip.text = element_text(size = rel(1), colour = "black"),
    axis.text.x = element_text(angle=45, hjust=1, vjust = 1.05)
  ) +
  # scale_y_continuous(breaks = c(0,0.5,1)) +
  ylab("Number of cells") + labs(fill = NULL)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DimReduc and composition: T-Cells
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pd = get_metadata(se.meta.t)
lbl = as.character(unique(pd$celltype))
lbl = setNames(gsub("\\.", " ", lbl), lbl)
t.reduc..pl =
  ggplot(data = pd, aes(x = umap_1, y = umap_2, col = celltype)) +
  scattermore::geom_scattermore(pointsize = 7, color="black")+
  scattermore::geom_scattermore(pointsize = 6, color="white") +
  scattermore::geom_scattermore(pointsize = 3) +
  theme(
    legend.position = "none",
    legend.margin = margin(l = 1),
    legend.title = element_text(margin = margin(b = 3), face = "plain"),
    legend.spacing.y = unit(2, 'mm'),
    legend.key.size = unit(4, "mm"),
    legend.text = element_text(margin = margin(l = leg.text.l, unit = "pt")),
    axis.title = element_text(size = 6),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    plot.title = element_blank()
  ) +
  guides(colour = guide_legend(
    title = NULL, ncol = 1, override.aes = list(shape = 16, size = 2.5)
  )) +
  scale_color_manual(values = til.col, na.value = "red", labels = lbl) +
  xlab("UMAP 1") + ylab("UMAP 2")

df = dittoBarPlot(
  se.meta.t, "celltype", group.by = "SAMPLE",
  split.by = "celltype_short_3", data.out = T
)[[2]]
df$grouping = factor(df$grouping, levels = names(sample.lbls))
df$celltype_short_3 = gsub("-Cell", "", df$celltype_short_3)
lvls = c(
  (levels(df$label)[grepl("CD4", levels(df$label))]),
  (levels(df$label)[grepl("CD8", levels(df$label))])
)
df$label = factor(df$label, levels = lvls)
lbls = as.character(unique(se.meta.t$celltype))
lbls = setNames(
  gsub("\\.", " ", gsub("_EOMES", "", lbls)),
  lbls
)
t.comp.pl =
  ggplot(df,  aes(grouping, count, fill = label)) +
  geom_col(aes(fill = label), position = "fill", width = .9) +
  facet_grid(~ celltype_short_3, scales="free", space = "free") +
  scale_fill_manual(values = til.col, labels = lbls) +
  theme(
    legend.position = "right",
    axis.title.x = element_blank(),
    legend.key.size = unit(9,"pt"),
    legend.margin = margin(l = 1),
    strip.text = element_text(size = rel(1), colour = "black"),
    axis.text.x = element_text(angle=45, hjust=1, vjust = 1.05)
  ) +
  scale_y_continuous(breaks = c(0,0.5,1)) +
  ylab("Percent of cells") + labs(fill = NULL)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DimReduc and composition: CAR+
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pd = get_metadata(se.meta.t)
pd = pd[order(pd$CAR_BY_EXPRS, decreasing = F), ]
car.reduc..pl =
  ggplot(data = pd, aes(x = umap_1, y = umap_2, col = CAR_BY_EXPRS)) +
  scattermore::geom_scattermore(pointsize = 7, color="black")+
  scattermore::geom_scattermore(pointsize = 6, color="white") +
  scattermore::geom_scattermore(pointsize = 3) +
  #  geom_point()
  theme(
    legend.position = "none",
    legend.margin = margin(l = 1),
    legend.title = element_text(margin = margin(b = 3), face = "plain"),
    legend.spacing.y = unit(2, 'mm'),
    legend.key.size = unit(4, "mm"),
    legend.text = element_text(margin = margin(l = leg.text.l, unit = "pt")),
    axis.title = element_text(size = 6),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    plot.title = element_blank()
  ) +
  guides(colour = guide_legend(
    title = "CAR+", rnow = 1, override.aes = list(shape = 16, size = 3)
  )) +
  scale_color_manual(
    values = c("#4477AA", "#DDDDDD"), breaks = c("TRUE", "FALSE"),
    na.value = "#FFFFFF"
  ) +
  xlab("UMAP 1") + ylab("UMAP 2")

df = dittoBarPlot(
  se.meta.t, "CAR_BY_EXPRS", group.by = "SAMPLE",
  split.by = "celltype_short_3", data.out = T
)[[2]]

df$grouping = factor(df$grouping, levels = names(sample.lbls))
df$celltype_short_3 = gsub("-Cell", "", df$celltype_short_3)
car.comp.pl =
  ggplot(df,  aes(grouping, count, fill = label)) +
  geom_col(aes(fill = label), position = "fill", width = .9) +
  facet_grid(~ celltype_short_3, scales="free", space = "free") +
  scale_fill_manual(
    values = c("#4477AA", "#DDDDDD"), breaks = c("TRUE", "FALSE"),
    na.value = "#FFFFFF"
  ) +
  theme(
    legend.position = "right",
    axis.title.x = element_blank(),
    legend.key.size = unit(9,"pt"),
    legend.margin = margin(l = 1),
    strip.text = element_text(size = rel(1), colour = "black"),
    axis.text.x = element_text(angle=45, hjust=1, vjust = 1.05)
  ) +
  scale_y_continuous(breaks = c(0,0.5,1)) +
  ylab("Percent of cells") + labs(fill = "CAR+")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DimReduc and composition: Clonotype groups
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pd = get_metadata(se.meta.t)
clono.col <- length(levels(pd$cloneSize))
clono.col = setNames(
  c("#003285", "#5AB2FF", "#DCDFA2", "#F08A5D", "#B83B5E", "#BBBBBB"),
  c(rev(levels(pd$cloneSize)), NA)
)
pd = pd[order(pd$clonalFrequency, decreasing = F, na.last = F), ]

clono.reduc.pl =
  ggplot(data = pd, aes(x = umap_1, y = umap_2, col = cloneSize)) +
  scattermore::geom_scattermore(pointsize = 7, color="black")+
  scattermore::geom_scattermore(pointsize = 6, color="white") +
  scattermore::geom_scattermore(pointsize = 3) +
  theme(
    legend.position = "none",
    legend.margin = margin(t=-5),
    legend.title = element_text(margin = margin(r = 3)),
    legend.spacing.y = unit(2, 'mm'),
    legend.key.size = unit(4, "mm"),
    legend.text = element_text(margin = margin(l = leg.text.l, unit = "pt")),
    axis.title = element_text(size = 6),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    plot.title = element_blank()
  ) +
  guides(colour = guide_legend(
    title = "Clonotype\ngroups", ncol = 2, override.aes = list(shape = 16, size = 3)
  )) +
  scale_color_manual(values = clono.col) +
  xlab("UMAP 1") + ylab("UMAP 2")


df = dittoBarPlot(
  se.meta.t, "cloneSize", group.by = "SAMPLE",
  split.by = "celltype_short_3", data.out = T
)[[2]]

df$grouping = factor(df$grouping, levels = names(sample.lbls))
df$celltype_short_3 = gsub("-Cell", "", df$celltype_short_3)
df$label = factor(df$label, levels = levels(se.meta.t$cloneSize))
lbls = as.character(unique(df$label))
lbls = setNames(gsub(" \\(", "\n\\(", lbls), lbls)
cl.comp.pl =
  ggplot(df,  aes(grouping, count, fill = label)) +
  geom_col(aes(fill = label), position = "fill", width = .9) +
  facet_grid(~ celltype_short_3, scales="free", space = "free") +
  scale_fill_manual(values = clono.col, labels = lbls) +
  theme(
    legend.position = "right",
    axis.title.x = element_blank(),
    legend.key.size = unit(8,"pt"),
    legend.margin = margin(l = 1),
    legend.text = element_text(margin = margin(t = 1.5, b = 1.5, unit = "pt"), size = rel(1)),
    strip.text = element_text(size = rel(1), colour = "black"),
    axis.text.x = element_text(angle=45, hjust=1, vjust = 1.05)
  ) +
  scale_y_continuous(breaks = c(0,0.5,1)) +
  ylab("Percent of cells") + labs(fill = "Clonotype groups")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Top Clonotypes
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cl.c = c(
  "#004488", "#DDAA33", "#5AB2FF",
  # "#004488", "#DDAA33", "#BB5566",
  "#117733", "#999933", "#CC6677", "#882255", "#AA4499", "#DDDDDD"
)

con.df = clonalCompare(
  se.meta.t,
  top.clones = 4, group.by = "SAMPLE",
  relabel.clones = T, cloneCall="strict",
  graph = "alluvial", exportTable = T,
)
con.df$clones = gsub("one:", "", con.df$clones)
con.df$clones = gsub("Clonotype\\:", "Cl", con.df$clones)

t = naturalsort(unique(con.df$clones))
con.df = con.df[!con.df$clones %in% t[c(length(t)-1, length(t))],  ]
con.df$Sample = factor(con.df$Sample, levels = names(sample.lbls))

top.clones.pl =
  ggplot(con.df, aes(
    x = Sample, fill = clones, group = clones,  stratum = clones,
    alluvium = clones, y = Proportion,  label = clones)
  ) +
  scale_fill_manual(values = cl.c) +
  theme(
    axis.title.x = element_blank(),
    panel.border = element_rect(colour = NA),
    axis.line = element_line(colour="black"),
    axis.text.x = element_text(angle=45, hjust=1, vjust = 1.05, size = rel(.9)),
    # legend.text = element_text(size = rel(1)),
    legend.key.size = unit(9,"pt"),
    legend.margin = margin(l = -1),
  ) +
  geom_stratum(linewidth = .1) +
  geom_flow(stat = "alluvium") +
  labs(fill = "Top clones") +
  ylab("Percent of cells")


pd = get_metadata(se.meta.t)
top.clones.reduc.pl  =
  ggplot(data = pd, aes(x = umap_1, y = umap_2, col = TopClones)) +
  scattermore::geom_scattermore(pointsize = 7, color="black")+
  scattermore::geom_scattermore(pointsize = 6, color="white") +
  geom_point(shape = ".") +
  # scattermore::geom_scattermore(pointsize = 2.8) +
  theme(
    legend.position = "none",
    legend.margin = margin(t=-5),
    legend.title = element_text(margin = margin(r = 3)),
    legend.spacing.y = unit(2, 'mm'),
    legend.key.size = unit(4, "mm"),
    legend.text = element_text(margin = margin(l = leg.text.l, unit = "pt")),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    plot.title = element_blank()
  ) +
  guides(colour = guide_legend(
    title = NULL, ncol = 3, override.aes = list(shape = 16, size = 2.5)
  )) +
  scale_color_manual(values = clone.col, na.value = "red")

ftrs = c("CD8A", "CD8B", "CD5", "CD28", "CD7", "TRAV29DV5", "TRBV5-1", "TRAV21", "TRBV27", "TRAV39")
mat = GetMatrixFromSeuratByGroupMulti(
  obj= se.meta.t, features = ftrs,
  CT_L1_COARSE, TopClones
)
exp_mat = reshape2::melt(mat$exp_mat)
colnames(exp_mat) = c("FTR", "GROUP", "AVE")
exp_mat$GROUP = gsub(".+\\|", "", exp_mat$GROUP)
percent_mat = reshape2::melt(mat$percent_mat)
colnames(percent_mat) = c("FTR", "GROUP", "PERC")
percent_mat$GROUP = gsub(".+\\|", "", percent_mat$GROUP)
df = cbind(exp_mat, PERC = percent_mat$PERC)
df$FTR = factor(df$FTR, levels = ftrs)
df$GROUP = factor(df$GROUP, levels = (naturalsort(unique(df$GROUP))))

p = scico(palette = "hawaii", n = 10, direction = -1)
top.clones.cd4cd8.pl =
  ggplot(df, aes(FTR, GROUP, fill = AVE, size = PERC * 100)) +
  geom_point(colour="black", pch=21, stroke = .2) +
  scale_size(range = c(.1, 3.5), breaks = c(20, 50, 100)) +
  mytheme() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle=45, hjust=1, vjust = 1.05, size = rel(.8)),
    axis.text.y = element_text(size = rel(.8)),
    # legend.text = element_text(margin = margin(l = 0, unit = "pt"), size = rel(1)),
    legend.title  = element_text(margin = margin(r = 5, l = 5, unit = "pt")),
    legend.margin = margin(t = -3, r = 40, l = -40),
    plot.title = element_text(hjust = 0.5, face = "plain"),
    legend.ticks.length = unit(0.05, 'cm')

  ) +
  scale_fill_gradientn(colors = c(p, rep(p[length(p)], 5))) +
  # scale_fill_scico(palette = "hawaii", direction = -1) +
  guides(
    fill = guide_colorbar(
      title = "Average\nExpression", order = 1,
      title.hjust = 0, title.vjust = .56, barwidth = unit(4, 'lines'),
      barheight = unit(.4, 'lines'), ticks.linewidth = 1.5/.pt
    ),
    size = guide_legend(title = "Percent\nExpressed")
  ) +
  xlab(NULL) + ylab(NULL) +
  coord_flip()

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA | Clone_1_3 vs Clone 2
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.work = se.meta.t
se.work = subset(se.work, orig.ident != "P2248" & orig.ident != "AphNB")
se.work = subset(se.work, TopClones_2 != "Cl_Other")
# Idents(se.work) = "TopClones_2"
# se.work = subset(se.work, downsample = 1000)
# Idents(se.work) = "orig.ident"
res.pb = run_wilx(
  obj = se.work, target = "celltype_short_3",
  min.pct.thres =  0.25, min.cells = 20, fc.thresh = 1.5,
  contrast.group = "TopClones_2", contrast = c("Cl_1_3", "Cl_2")
)
res.pb.sign = subset(res.pb, significant == T)
table(res.pb.sign$cluster)

dot.clono.de =
  de_tops_bubble(
    df = res.pb.sign[!grepl("^TRA|^TRB", res.pb.sign$feature), ],
    se.obj = se.work,
    ctrst.lvl = c("Cl_2", "Cl_1_3"),
    pl.title = paste0(
      "Clone 1_3 vs. Clone 2 (up: ",
      table(res.pb.sign$avg_log2FC > 0)[2],
      ", down: ",
      table(res.pb.sign$avg_log2FC > 0)[1],
      ")"
    ),
    ctrst = "TopClones_2"
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ORA | Clone_1_3 vs Clone 2
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ftrs.l = res.pb.sign
# ftrs.l$cluster = ftrs.l$avg_log2FC > 0
# ftrs.l$cluster = ifelse(ftrs.l$cluster == T, "LFC>0", "LFC<0")
ftrs.l = split(ftrs.l, ftrs.l$cluster)
ftrs.l = lapply(ftrs.l, function(x){
  ftrs = x$feature
  names(ftrs) = x$avg_log2FC
  ftrs
})

ora.go = parallel::mclapply(ftrs.l, function(x){
  run_nmf_ora(
    genes = x,
    universe = rownames(se.meta.t),
    # category = "C2", db.sub = "REACTOME"
    category = "C5", subcategory = "BP"
  )
}, mc.cores = length(ftrs.l))

# Redundant term (most gene overlap)
ora.go[[1]] = ora.go[[1]][!ora.go[[1]]$pathway %in%
                            c(
                              "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II" ,
                              "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_OR_POLYSACCHARIDE_ANTIGEN_VIA_MHC_CLASS_II"
                            ), ]

ora.clono =
  ora_barpl(
    gsea.res = ora.go,
    ftrs.list = ftrs.l,
    nbr.tops = 10,
    font.size = 8,
    min.genes = 5,
    max.value = 1.5,
    term.length = 70,
    sort.by.padj = T,
    barwidth = unit(.4, 'lines'),
    barheight = unit(4, 'lines'),
    bar.width = .8
  ) +
  ggtitle("Cl 1_3 vs. Cl 2") +
  theme(axis.text.x = element_text(size = rel(1)))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Spatial
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
theme.custom = mytheme(base_size = 8) &
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.title = element_blank(),
    panel.border = element_blank(),
    plot.subtitle = element_text(hjust = 0.5, size = rel(1), face = "plain"),
    plot.title = element_blank(),
    legend.title = element_blank(),
  )

work.base = manifest$spatial$cellranger

sp.work = paste0(work.base, "C24_13076_Sonden_SpikeIn_one_probe/outs/") # Köln_spike-in Carvyctysonden
samples =  file.path(sp.work, "raw_feature_bc_matrix.h5")
imgs = file.path(sp.work,  "spatial", "tissue_lowres_image.png")
spotfiles = file.path(sp.work,"spatial", "tissue_positions.csv")
json = file.path(sp.work, "spatial", "scalefactors_json.json")
infoTable <- tibble(samples, imgs, spotfiles, json, sample_id = c("Köln_spike-in"))
se.1 <- ReadVisiumData(infoTable)
se.1 <- LoadImages(se.1) %>% NormalizeData()


ftrs =  c("CD3E", "CD4", "CD8B", "ciltacel")
global.max = max(log10(unlist(FetchData(se.1, ftrs, layer = "counts")) + 1))
ftr.l = spatial_features(
  .obj = se.1, features = ftrs, pt.size = .01, min.max = c(0, global.max),
  legend.wh = c(4, .4), base.size = 8, slot = "counts", log.counts = T,
  .colors = scico::scico(30, palette = "glasgow", direction = -1), plot.grid = F
)
leg.1 = get_legend(
  ftr.l[[4]] + theme(legend.position = "bottom") + guides(
    color = guide_colorbar(
      title = "log10(counts + 1)", title.vjust = 1.05, barwidth = 4,
      barheight = .4, ticks.linewidth = 1.5/.pt
    )
  ) + theme(
    legend.title = element_text(margin = margin(r = 3, unit = "pt")),
    legend.ticks.length = unit(0.05, 'cm')
  )
)
ftr.l = lapply(ftr.l, function(x){
  x = x + theme(legend.position='none')
})
ftr.l.1 = ftr.l

se_masked <- se.1 |>
  MaskImages(minPixels = 0, method = "blurSeg")
se_masked = ImagePlot(se_masked, return_as_gg = T) & theme.custom
# se_masked = se_masked + plot_annotation(subtitle = 'Skin')
se_masked.1 = se_masked

ftrs =  c("TRAV29DV5", "TRBV5-1", "TRAV21", "TRBV27", "TRAV39")
global.max = max(log10(unlist(FetchData(se.1, ftrs, layer = "counts")) + 1))
ftr.l = spatial_features(
  .obj = se.1, features = ftrs, pt.size = .01, min.max = c(0, global.max),
  legend.wh = c(4, .4), base.size = 8, slot = "counts", log.counts = T,
  .colors = scico::scico(30, palette = "glasgow", direction = -1), plot.grid = F
)
leg.2 = get_legend(
  ftr.l[[1]] + theme(legend.position = "bottom") + guides(
    color = guide_colorbar(
      title = "log10(counts + 1)", title.vjust = 1.05, barwidth = 4,
      barheight = .4, ticks.linewidth = 1.5/.pt
    )
  ) + theme(
    legend.title = element_text(margin = margin(r = 3, unit = "pt")),
    legend.ticks.length = unit(0.05, 'cm')
  )
)
ftr.l = lapply(ftr.l, function(x){
  x = x + theme(legend.position='none')
})
ftr.l.2 = ftr.l

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Final Plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ggsave2(
  filename="publication/figures_tables/main/fig_02.png",
  plot_grid(
    # NULL,
    plot_grid(
      plot_grid(
        NULL,
        plot_grid(NULL, all.ct.pl, NULL, ncol = 3, rel_widths = c(.25, 1, .25)),
        NULL,
        comp.pl, nrow = 4, rel_heights = c(1.1, .8, .0, 1)
      ),
      NULL,
      plot_grid(
        plot_grid(NULL, t.reduc..pl, NULL, ncol = 3, rel_widths = c(.1, 1, .1)),
        NULL,
        t.comp.pl, nrow = 3, rel_heights = c(1, 0, 1.3)
      ),
      NULL,
      plot_grid(
        plot_grid(NULL, car.reduc..pl, NULL, ncol = 3, rel_widths = c(.07, 1, .07)),
        NULL,
        car.comp.pl, nrow = 3, rel_heights = c(1, 0, 1.3)
      ),
      NULL,
      plot_grid(
        plot_grid(NULL, clono.reduc.pl, NULL, ncol = 3, rel_widths = c(.12, 1, .12)),
        NULL,
        cl.comp.pl, nrow = 3, rel_heights = c(1, 0, 1.3)
      ),
      NULL,
      ncol = 8, rel_widths = c(.75, .15, 1, .1, .85, .1, 1.05, .03),
      labels = c("A)", "", "B)", "", "C)", "", "D)"),
      label_fontface = "bold", label_size = 10
    ),
    NULL,
    # NULL,
    plot_grid(
      plot_grid(
        top.clones.pl + theme(plot.margin = margin(15, -5, 4, 4, unit = "pt")),
        NULL,
        plot_grid(
          plot_grid(
            NULL,
            plot_grid(top.clones.reduc.pl, NULL, nrow = 2, rel_heights = c(1, .6)),
            top.clones.cd4cd8.pl,
            ncol = 3, rel_widths = c(0, 1, 1.1),
            labels = c("F)", "", "G)"), label_y = 1.07,
            label_fontface = "bold", label_size = 10
          ),
          NULL, nrow = 2, rel_heights = c(1, 0)
        ),
        nrow = 3, rel_heights = c(1, .11, 1.1)
      ),
      NULL,
      plot_grid(
        dot.clono.de,
        NULL,
        ora.clono,
        nrow = 3, rel_heights = c(1, .1 , 1.25),
        labels = c("H)", "", "I)"),
        label_fontface = "bold", label_size = 10
      ),
      NULL,
      ncol = 4, rel_widths = c(.382, .075, .618, .0),
      labels = c("E)", "", "", ""),
      label_fontface = "bold", label_size = 10
    ),
    NULL,
    plot_grid(
      plot_grid(
        se_masked.1 & theme(plot.margin = margin(2, -10, 2, 0, unit = "pt")),
        plot_grid(
          plot_grid(plotlist = ftr.l.1, nrow = 1),
          leg.1, nrow = 2, rel_heights = c(1, .125)
        ),
        nrow = 1, rel_widths = c(1, 3.4)
      ),
      NULL,
      plot_grid(
        plot_grid(plotlist = ftr.l.2, nrow = 1),
        leg.2, nrow = 2, rel_heights = c(1, .125)
      ),
      ncol = 3, rel_widths = c(1.05, 0.1, 1),
      labels = c("J)","","K)"),label_fontface = "bold", label_size = 10
    ),
    nrow = 5, rel_heights = c(.95, .05, .95, .05, .425)
  ),
  width = 180,
  height = 165,
  dpi = 400,
  bg = "white",
  units = "mm",
  scale = 1.6,
  device = png, type = "cairo"
)

