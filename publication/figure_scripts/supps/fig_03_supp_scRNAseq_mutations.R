# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Read packages
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c(
  "tidyverse", "stringr","stringi", "readr"
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

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Styles
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
source("./code/helper/styles.R")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Read and transform data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
df=readRDS("./data/sc_mutation_data.RDS")
df=df[!(df$sample=="AphNB" & df$cells=="Cl_Other"),]

df = aggregate(cbind(n_total, n_mut) ~ cells, data = df, FUN = sum)
df <- df %>%
  mutate(Sample = recode(cells, "Cl_1_3" = "Clone 1_3", "Cl_2" = "Clone 2", "Cl_Other" = "Other clones", "complete" = "Apheresis"))
plot.df=df
plot.df$Sample=ordered(plot.df$Sample,levels=c("Apheresis","Other clones","Clone 2","Clone 1_3"))
plot.df=plot.df[order(plot.df$Sample),]
plot.df$VAF=plot.df$n_mut/plot.df$n_total
plot.df$VAF_txt = paste0(plot.df$n_mut,"/",plot.df$n_total)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Plot Fig S23
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
scrna_mut = ggplot(plot.df, aes(y="",x = VAF, fill = Sample, group = Sample)) +
  geom_vline(xintercept = 0.5) +
  geom_vline(xintercept = 1.) +
  geom_col(position = position_dodge(), color = "black", linewidth=0.2) +
  ylab(NULL) +
  xlab("Variant Allele Fraction") +
  scale_x_continuous(limits = c(0., 1.), expand = c(0., 0.)) +
  scale_fill_manual(values = c("#DDDDDD", "#AACC88", "#DDAA33", "#6699CC"),
                    guide = guide_legend(nrow = 1, reverse = T,
                                         title = "")) +
  mytheme_grid(base_size = 10) +
  theme(legend.position = "bottom",
        legend.margin = margin(0,0,0,-3, "cm"),
        plot.margin = margin(0.3,0.6,0.3,0.3, "cm"),
        plot.background = element_rect(linewidth = 0),
        panel.grid.major.y = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.text.y = element_text(face = "italic"))+
  ggtitle(label="TET2 p.R544*") +
  theme(plot.title = element_text(size=14)) +
  geom_text(aes(label=VAF_txt,x=0.08,y=c(seq(0.66,1.35,by=0.221))))
scrna_mut

ggsave2(
  filename="publication/figures_tables/supps/fig_S23.png",
  scrna_mut,
  width = 100,
  height = 40,
  dpi = 400,
  bg = "white",
  units = "mm",
  scale = 1.6
)
