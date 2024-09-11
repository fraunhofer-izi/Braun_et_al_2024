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
library(dplyr)
library(stringr)
library(grid)

## known genome region of KPNA4 = chr3:160,495,007-160,565,571
# find table information on https://genome.ucsc.edu/cgi-bin/hgTables
# obtain output: All fields from selected table

# load table

## testing
table <- '~/uranos_home/workspace/2024-BCMA-CAR-Koeln/gviz/KPNA4'
transcrip <- 2
chr <- "chr3"
gene <- "KPNA4_test"
IS <- 63719663
         

          
myGene <- function(table, transcrip, chr, flank1, flank2, gene, IS){
  knowngene <- read.delim(as.character(table), header = FALSE)
  Exons <- knowngene[transcrip,8] %>%
    as.numeric()
  exonStarts <- knowngene[transcrip,9] %>%
    str_split_1(pattern = ',') %>%
    as.numeric() %>% 
    na.omit() 
  exonEnds <- knowngene[transcrip,10] %>%
    str_split_1(pattern = ',') %>%
    as.numeric() %>% 
    na.omit()
  exonWidths <- exonEnds - exonStarts
  strand <- knowngene[transcrip,3]
  # Ex <- c(rep("Exon", Exons))
  number <- 1:Exons 
  if (strand == "-"){
    number <- rev(number)
    }
  atrack <- AnnotationTrack(start = exonStarts,
                            width = exonWidths,
                            chromosome = chr,
                            strand = rep(strand,Exons), 
                            genome = "hg38",
                            # name = paste0(chr, ": ", gene, collapse = ""),
                            group = rep("Exon",Exons),
                            fontsize = 25
                            )
  feature(atrack) <- as.character(number)
  if (strand == "-"){
    f1 <- Exons - flank2 +1
    f2 <- Exons - flank1 +1
  }else{
    f1 <- flank1
    f2 <- flank2
  }
  ht <- HighlightTrack(trackList=atrack, 
                       start = exonStarts[f1], 
                       width = exonEnds[f2]- exonStarts[f1], 
                       chromosome = chr
                       )
  svg(paste0(gene, "_gene_v3.svg", collapse = " "),width = 12, height = 3)
  # Code of the plot
  plotTracks(ht, 
             # featureAnnotation = "feature",
             groupAnnotation  = "group",
             # fontcolor.feature = 1, 
             # cex.feature = 0.7,
             lwd = 1.5,
             fontsize = 20,
             cex.group = 3,
             fontcolor.group = 'black',
             lineheight = 1,
             # col.title = 'lightgrey',
             # rotation.title = 0,
             showTitle = FALSE,
             collapse = FALSE,
             # main = paste0(chr, ": ", gene, collapse = " "),
             cex.main = 2.5
             )
  # Close the graphics device
  dev.off()
}

myGene('~/uranos_home/workspace/2024-BCMA-CAR-Koeln/gviz/KPNA4', 2, "chr3", 2, 3, "KPNA4", 160535998)
myGene('~/uranos_home/workspace/2024-BCMA-CAR-Koeln/gviz/ZGPAT', 10, "chr20", 2, 3, "ZGPAT", 63719663) 
myGene('~/uranos_home/workspace/2024-BCMA-CAR-Koeln/gviz/Polycomb-Associated_Non-Coding_RNAs', 4, "chr19", 3, 4, "MIDN", 1246860)
myGene('~/uranos_home/workspace/2024-BCMA-CAR-Koeln/gviz/Polycomb-Associated_Non-Coding_RNAs', 2, "chr19", 3, 4, "lncRNA", 1246860)
myGene('~/uranos_home/workspace/2024-BCMA-CAR-Koeln/gviz/Polycomb-Associated_Non-Coding_RNAs', 2, "chr19", 3, 4, "Polycomb-Associated Non-Coding RNAs", 1246860)
