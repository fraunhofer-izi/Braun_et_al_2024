.cran_packages = c(
  "yaml", "ggplot2", "reshape2", "dplyr", "ggpubr", "tidyverse",
   "cowplot", "patchwork", "ggh4x", "ggbeeswarm"
)
.bioc_packages = c()

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

Read_Metrics_10X_Multi <- function(
  base_path,
  secondary_path = "outs/per_sample_outs/",
  lib_list = NULL,
  col_values = "Metric.Value"
) {

  library(cli)
  library(dplyr)

  # Confirm directory exists
  if (dir.exists(paths = base_path) == FALSE) {
    cli_abort(message = "Directory: {.val {base_path}} specified by {.code base_path} does not exist.")
  }
  # Detect libraries if lib_list is NULL
  if (is.null(x = lib_list)) {
    lib_list <- list.dirs(path = base_path, full.names = F, recursive = F)
  }

  # Check if full directory path exists
  for (i in 1:length(x = lib_list)) {
    full_directory_path <- file.path(base_path, lib_list[i], secondary_path)
    if (dir.exists(paths = full_directory_path) == FALSE) {
      cli_abort(message = "Full Directory does not exist {.val {full_directory_path}} was not found.")
    }
  }

  raw_data_list <- lapply(1:length(x = lib_list), function(x) {
    file_path <- file.path(base_path, lib_list[x], paste0(secondary_path, "/", lib_list[x], "/"))
    raw_data <- read.csv(file = paste0(file_path, "metrics_summary.csv"), stringsAsFactors = F)
    return(raw_data)
  })
  names(raw_data_list) <- lib_list
  full_data <- dplyr::bind_rows(raw_data_list, .id = "sample_id")

  # make value numeric
  full_data$quantity = ifelse(grepl("%", full_data$Metric.Value), "PERC", "COUNT" )
  row_numbers <- grep(pattern = ",", x = full_data[, col_values])
  full_data[row_numbers, ][, col_values] = as.numeric(gsub(",", "", full_data[row_numbers, ][, col_values]))
  full_data[, col_values] = gsub("%", "", full_data[, col_values])
  full_data[, col_values] = as.numeric(full_data[, col_values])
  full_data$Metric.Name = paste0(full_data$Metric.Name, "_", full_data$quantity)

  raw.gex = full_data[full_data$Library.Type == "Gene Expression", ]
  raw.gex.cells = raw.gex[raw.gex$Category == "Cells", ]
  raw.gex.cells = raw.gex.cells[!grepl("Confidently", raw.gex.cells$Metric.Name), ]
  raw.gex.seq = raw.gex[raw.gex$Grouped.By == "Fastq ID", ]
  raw.gex.mapping = raw.gex[grepl("Confidently", raw.gex$Metric.Name), ]
  raw.gex.mapping = subset(raw.gex.mapping, Category != "Cells")
  raw.gex.physical = raw.gex[!rownames(raw.gex) %in% rownames(rbind(raw.gex.cells, raw.gex.seq, raw.gex.mapping)), ]
  raw.gex.physical = subset(raw.gex.physical, Metric.Name != "Estimated number of cells_COUNT") # Redundant

  raw.vdj.t = full_data[full_data$Library.Type == "VDJ T", ]
  raw.vdj.t.cells = raw.vdj.t[raw.vdj.t$Category == "Cells", ]
  raw.vdj.t.cells.vdj.anno =  raw.vdj.t.cells[grepl("^Cells|Paired clonotype", raw.vdj.t.cells$Metric.Name), ]
  raw.vdj.t.cells.texprs = raw.vdj.t.cells[!rownames(raw.vdj.t.cells) %in% rownames(raw.vdj.t.cells.vdj.anno), ]
  raw.vdj.t.seq = raw.vdj.t[raw.vdj.t$Grouped.By == "Fastq ID", ]
  raw.vdj.t.ph = raw.vdj.t[!rownames(raw.vdj.t) %in% rownames(rbind(raw.vdj.t.cells, raw.vdj.t.seq)), ]
  raw.vdj.t.ph = subset(raw.vdj.t.ph, Metric.Name != "Estimated number of cells_COUNT") # Redundant
  raw.vdj.t.mapping = raw.vdj.t.ph[grepl("Reads mapped", raw.vdj.t.ph$Metric.Name), ]
  raw.vdj.t.physical = raw.vdj.t.ph[!grepl("Reads mapped", raw.vdj.t.ph$Metric.Name), ]

  res = list(
    "GEX" = list(
      "Cell Metrics" = raw.gex.cells[, c("sample_id", "Metric.Name", "Metric.Value")] %>% tidyr::spread(Metric.Name,  Metric.Value),
      "Sequencing Metrics" = raw.gex.seq[, c("sample_id", "Metric.Name", "Metric.Value")] %>% tidyr::spread(Metric.Name,  Metric.Value),
      "Mapping Metrics" = raw.gex.mapping[, c("sample_id", "Metric.Name", "Metric.Value")] %>% tidyr::spread(Metric.Name,  Metric.Value),
      "Metrics Per Physical Library" = raw.gex.physical[, c("sample_id", "Metric.Name", "Metric.Value")] %>% tidyr::spread(Metric.Name,  Metric.Value)
    ),
    "VDJ_T" = list(
      "T Cell Expression" = raw.vdj.t.cells.texprs[, c("sample_id", "Metric.Name", "Metric.Value")] %>% tidyr::spread(Metric.Name,  Metric.Value),
      "V(D)J Annotation" = raw.vdj.t.cells.vdj.anno[, c("sample_id", "Metric.Name", "Metric.Value")] %>% tidyr::spread(Metric.Name,  Metric.Value),
      "Sequencing Metrics" = raw.vdj.t.seq[, c("sample_id", "Metric.Name", "Metric.Value")] %>% tidyr::spread(Metric.Name,  Metric.Value),
      "Mapping Metrics" = raw.vdj.t.mapping[, c("sample_id", "Metric.Name", "Metric.Value")] %>% tidyr::spread(Metric.Name,  Metric.Value),
      "Metrics Per Physical Library" = raw.vdj.t.physical[, c("sample_id", "Metric.Name", "Metric.Value")] %>% tidyr::spread(Metric.Name,  Metric.Value)
    )
  )

  return(res)
}

plot_func = function(
  df,
  pl.title = "Title",
  dot.col = NULL,
  group.by = NULL,
  shape.by = NULL
) {

  tmp = unlist(lapply(colnames(df), function(x){paste(strwrap(x, width = 20), collapse = "\n")}))
  colnames(df) = tmp

  df.pl = reshape2::melt(df)

  if(!is.null(group.by)){

    if(!is.null(dot.col)){
      d.col = dot.col
    } else {
      d.col = group.by
    }

    if(!is.null(shape.by)){
      d.shape = shape.by
    } else {
      df.pl$TMP = "1"
      d.shape = "TMP"
    }

    pl =
      ggplot(df.pl, aes(x = .data[[group.by]], y = value, col = .data[[d.col]], shape = .data[[d.shape]])) +
      geom_jitter(position=position_jitter(0.1, seed = 1234), size = 2) +
      stat_summary(
        aes(group = .data[[group.by]] ),
        fun.min = function(z) { quantile(z,0.25) },
        fun.max = function(z) { quantile(z,0.75) },
        fun = median, colour = "black", size = .2, linewidth = .2, show.legend = F
      ) +
      # geom_beeswarm(cex = 5.5, method = "center") +
      scale_color_manual(values = ggthemes::tableau_color_pal("Tableau 10")(10)) +
      facet_wrap2(~ variable, scales = "free", nrow = 1)
  } else {

    if(is.null(dot.col)) {
      pl = ggplot(df.pl, aes(x = sample_id, y = value)) +
        geom_point()
    } else {
      pl = ggplot(df.pl, aes(x = sample_id, y = value, col = .data[[dot.col]])) +
        geom_point()
    }

  }


  pl = pl +
    theme(
      aspect.ratio = 1,
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.position = "bottom",
      # legend.title = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = rel(1)),
      axis.ticks = element_blank(),
      panel.spacing = unit(1.5, "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold", colour = "black", size = rel(1.5)),
      strip.text = element_text(size = rel(1), face = "plain"),
      strip.background = element_blank()
    ) +
    facet_wrap(~ variable, scales = "free", nrow = 1) +
    labs(title = pl.title) + ylab(NULL)

  pl
}

source("code/helper/styles.R")
theme_set(mytheme_grid(base_size = 8))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")
metrics = Read_Metrics_10X_Multi(base_path = paste0(manifest$multi_omics$cellranger))

metrics.gex = lapply(metrics$GEX, function(x){
  x$sample_id = gsub("multi_", "", x$sample_id)
  x$sample_id = gsub("AphNB", "Aph", x$sample_id)
  x$sample_id = gsub("P2248", "BM", x$sample_id)
  x$sample_id = gsub("P2249", "PB", x$sample_id)
  x$sample_id = gsub("P2264", "PB+Dexa", x$sample_id)
  colnames(x) = gsub("_PERC", " (%)", gsub("_COUNT", "", colnames(x)))
  colnames(x) = gsub("Confidently", "Conf.", colnames(x))
  x
})

metrics.gex[[3]]$`Conf. mapped antisense (%)` = NULL
metrics.gex[[2]]$`Number of short reads skipped` = metrics.gex[[4]]$`Sequencing saturation (%)`
colnames(metrics.gex[[2]])[3] = "Sequencing saturation (%)"

a = plot_func(df = metrics.gex[[1]], pl.title = names(metrics.gex)[1])
b = plot_func(metrics.gex[[2]], names(metrics.gex)[2])
c = plot_func(metrics.gex[[3]], names(metrics.gex)[3])
# d = plot_func(metrics.gex[[4]], names(metrics.gex)[4])

# ggsave2(
#   plot_grid(a, b, c, ncol = 1, align = "vh", scale = .95),
#   filename = "analysis/multi_omics_mir/figures/secondary_analysis/gex.png",
#   units = "mm", width = 160, height = 110, scale = 1.7
# )

ggsave2(
  plot_grid(a, b, c, ncol = 1, align = "vh", scale = .95),
  filename = "analysis/multi_omics_mir/figures/secondary_analysis/gex.pdf",
  units = "mm", width = 160, height = 140, scale = 1.6
)

metrics.vdj = lapply(metrics$VDJ_T, function(x){
  x$sample_id = gsub("multi_", "", x$sample_id)
  colnames(x) = gsub("_PERC", " (%)", gsub("_COUNT", "", colnames(x)))
  x
})

a = plot_func(metrics.vdj[[1]], names(metrics.vdj)[1])
b = plot_func(metrics.vdj[[2]], names(metrics.vdj)[2])
c = plot_func(metrics.vdj[[3]], names(metrics.vdj)[3])
d = plot_func(metrics.vdj[[4]], names(metrics.vdj)[4])
e = plot_func(metrics.vdj[[5]], names(metrics.vdj)[5])

ggsave2(
  plot_grid(a, b, ncol = 1, align = "vh", scale = .95),
  filename = "analysis/multi_omics_mir/figures/secondary_analysis/vdj_tcell_1.pdf",
  units = "mm", width = 230, height = 180, scale = 1.1
)

ggsave2(
  plot_grid(c, d, e, ncol = 1, align = "vh", scale = .95),
  filename = "analysis/multi_omics_mir/figures/secondary_analysis/vdj_tcell_2.pdf",
  units = "mm", width = 220, height = 230, scale = 1.1
)


