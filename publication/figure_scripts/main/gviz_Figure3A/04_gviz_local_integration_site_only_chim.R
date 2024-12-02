# install necessary packages
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

BiocManager::install("Gviz")
BiocManager::install("GenomicRanges")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

# load necessary packages
library(Gviz)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(stringr)


# setwd('segemehl_realigned/') #this can be changed
# setwd('IntegrationSite_coverage/')

# set global environment objects that can be used in the functions

gen = 'hg38'
chr01 = 'chr3'
chr02 = 'chr19'
chr03 = 'chr20'
color1 = '#6699CC'
color2 = '#DDAA33'
color3 = '#AACC88'
color4 = '#DDDDDD'
pop1 <- 'CARpos_CD5dim'
pop2 <- 'CARpos_CD5pos'
pop3 <- 'CARneg_CD5pos'
pop4 <- 'Apherese'
label.1 <- 'CAR+CD5-'
label.2 <- 'CAR+CD5+'
label.3 <- 'CAR-CD5+'
label.4 <- 'Apheresis'
chim <- 'unique' 
all <- '5k'
chim_coverage = 27
chim_tick = 5
coverage = 90
all_tick = 25
folder <- 'Braun_et_al_2024/data/vis_reads/combine_runs/'

file_names <- read.table(paste0(c(folder,
                                  'bam_list'),
                                  collapse = "")
                         ) %>%
  unlist
names(file_names) <- NULL
print(file_names)

## load bam files

myBam <- function(bam, sample, chr, gr_color){
  DataTrack(range = bam, 
            genome = gen, 
            type = "l", 
            name = sample, 
            window = -1,
            # windowSize = 4,
            chromosome = chr,
            rotation.title = 90,
            col.axis = "black",
            cex.axis = 1,
            showTitle = TRUE,
            col.title = 'black',
            cex.title = 0.9,
            frame = FALSE,
            background.title = gr_color,
            )
}

# example

# population <- 'CARpos_CD5dim'
# chr <- chr01
# chim <- all # or chim
# label <- 'CD5-CAR+' # y-axis
# bgcol <- color1

myChoice <- function(input_list, folder, population, chr, chim, label, bgcol){
  for (i in 1:length(input_list)){
    if (grepl(population, input_list[i]) && grepl(chr, input_list[i]) && grepl(chim, input_list[i])){
      dtrack <- myBam(paste0(c(folder, input_list[i]),
                              collapse = ""),
                       label,
                       chr,
                       bgcol)
      if (grepl('unique', input_list[i])){
        displayPars(dtrack) <- list(ylim = c(0,chim_coverage),
                                     yTicksAt = c(seq(0, chim_coverage, chim_tick))
        )
      }else{
        displayPars(dtrack) <- list(ylim = c(0,coverage), yTicksAt = c(seq(0, coverage, all_tick)))
      }
    }
  }
  return(dtrack)
}

############

S2_chr3 <- myChoice(file_names, folder, pop1, chr01, all, label.1, color1)
S2_chr19 <- myChoice(file_names, folder, pop1, chr02, all, label.1, color1)
S2_chr20 <- myChoice(file_names, folder, pop1, chr03, all, label.1, color1)

S2_chr3_chim <- myChoice(file_names, folder, pop1, chr01, chim, label.1, color1)
S2_chr19_chim <- myChoice(file_names, folder, pop1, chr02, chim, label.1, color1)
S2_chr20_chim <- myChoice(file_names, folder, pop1, chr03, chim, label.1, color1)

S3_chr3 <- myChoice(file_names, folder, pop2, chr01, all, label.2, color2)
S3_chr19 <- myChoice(file_names, folder, pop2, chr02, all, label.2, color2)
S3_chr20 <- myChoice(file_names, folder, pop2, chr03, all, label.2, color2)

S3_chr3_chim <- myChoice(file_names, folder, pop2, chr01, chim, label.2, color2)
S3_chr19_chim <- myChoice(file_names, folder, pop2, chr02, chim, label.2, color2)
S3_chr20_chim <- myChoice(file_names, folder, pop2, chr03, chim, label.2, color2)

S4_chr3 <- myChoice(file_names, folder, pop3, chr01, all, label.3, color3)
S4_chr19 <- myChoice(file_names, folder, pop3, chr02, all, label.3, color3)
S4_chr20 <- myChoice(file_names, folder, pop3, chr03, all, label.3, color3)

S4_chr3_chim <- myChoice(file_names, folder, pop3, chr01, chim, label.3, color3)
S4_chr19_chim <- myChoice(file_names, folder, pop3, chr02, chim, label.3, color3)
S4_chr20_chim <- myChoice(file_names, folder, pop3, chr03, chim, label.3, color3)
## Aphrese not included yet
# S1_chr3 <- myChoice(file_names, folder, pop4, chr01, all, label.4, color4)
# S1_chr19 <- myChoice(file_names, folder, pop4, chr02, all, label.4, color4)
# S1_chr20 <- myChoice(file_names, folder, pop4, chr03, all, label.4, color4)
# 
# S1_chr3_chim <- myChoice(file_names, folder, pop4, chr01, chim, label.4, color4)
# S1_chr19_chim <- myChoice(file_names, folder, pop4, chr02, chim, label.4, color4)
# S1_chr20_chim <- myChoice(file_names, folder, pop4, chr03, chim, label.4, color4)

## create object: genomic coordinates as reference, custom build
### highlight regions
Integration_chr3 <- GenomeAxisTrack(name = c("Integration Site"),
                                    showTitle = TRUE,
                                    range=IRanges(start = c(160533000,
                                                            160536003,
                                                            160541696),
                                                  end = c(160536002,
                                                          160541695,
                                                          160545000),
                                                  names = c("KPNA4:intron 2",
                                                            "cilta-cel vector",
                                                            "KPNA4:intron 2")))

Integration_chr19 <- GenomeAxisTrack(name = c("Integration Site"),
                                     showTitle = TRUE,
                                     range=IRanges(start = c(1243700,
                                                             1246861,
                                                             1252554),
                                                   end = c(1246860,
                                                           1252553,
                                                           1255700),
                                                   names = c("lncRNA:intron 3",
                                                             "cilta-cel vector",
                                                             "lncRNA:intron 3")))

Integration_chr20 <- GenomeAxisTrack(name = c("Integration Site"),
                                     showTitle = TRUE,
                                     range=IRanges(start = c(63716000,
                                                             63719666,
                                                             63725359),
                                                  end = c(63719665,
                                                          63725358,
                                                          63728000),
                                                  names = c("ZGPAT:intron 2",
                                                            "cilta-cel vector",
                                                            "ZGPAT:intron 2")))

# highlight boxes: need new chr-Carvykti annotation
ht1 <- HighlightTrack(trackList=list(Integration_chr3,
                                     # S2_chr3,
                                     S2_chr3_chim,
                                     # S3_chr3,
                                     S3_chr3_chim,
                                     # S4_chr3,
                                     S4_chr3_chim
                                     ),
                      start = c(160535753, 160541445),
                      width = 500,
                      chromosome = chr01)

ht2 <- HighlightTrack(trackList=list(Integration_chr19,
                                     # S2_chr19,
                                     S2_chr19_chim,
                                     # S3_chr19,
                                     S3_chr19_chim,
                                     # S4_chr19,
                                     S4_chr19_chim
                                     ),
                      start = c(1246611, 1252303),
                      width = 500,
                      chromosome = chr02)

ht3 <- HighlightTrack(trackList=list(Integration_chr20,
                                     # S2_chr20,
                                     S2_chr20_chim,
                                     # S3_chr20,
                                     S3_chr20_chim,
                                     # S4_chr20,
                                     S4_chr20_chim
                                     ), 
                      start = c(63719466, 63725158),
                      width = 500,
                      chromosome = chr03)

##############
myPNG <- function(chr, section, fr, t){
  svg(paste0(chr, "_Coverage.svg", collapse = " "), width = 10, height = 10)
  # Code of the plot
  plotTracks(section,
             from = fr,
             to = t,
             add53 = TRUE,
             add35 = TRUE, 
             littleTicks = FALSE, 
             showId = TRUE,
             showTitle = TRUE,
             cex = 0.9, 
             cex.axis = 1.3,
             cex.id = 1.1, 
             col.id = c('blue', 'red', 'blue'),
             fill.range = c('grey', 'white', 'grey'), 
             # col.range = 'black', 
             exponent = 6, 
             fontcolor = 'black', 
             fontsize = 20,
             rotation.title = 90,
             cex.title = 1.2,
             lwd = 2,
             labelPos = 'alternating',
             baseline = 0,
             col.baseline = 'black'
             )
# Close the graphics device
dev.off()
}

myPNG("chr3_chim", ht1, 160533000, 160545000)
myPNG("chr19_chim", ht2, 1243500, 1255500)
myPNG("chr20_chim", ht3, 63717000, 63729000)
