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

## known genome region of KPNA4 = chr3:160,495,007-160,565,571
# find table information on https://genome.ucsc.edu/cgi-bin/hgTables
# obtain output: All fields from selected table

# load table

## testing
# table <- 'Polycomb-Associated_Non-Coding_RNAs'
# transcrip <- 3
# chr <- "chr20"
# gene <- "ZGPAT_test"
# IS <- 63719663
# flank1 = 2
# flank2 = 3

          
myZoom <- function(table, transcrip, chr, flank1, flank2, label, gene, IS){
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
  if (strand == "-"){
    f1 <- Exons - flank1 + 1
    f2 <- Exons - flank2 + 1
    exonStarts <- exonStarts[f2:f1]
    # exonEnds <- exonEnds[f2:f1]
    exonWidths <- exonWidths[f2:f1]
  } else{
    f1 <- flank1
    f2 <- flank2
    exonStarts <- exonStarts[f1:f2]
    # exonEnds <- exonEnds[f1:f2]
    exonWidths <- exonWidths[f1:f2]
  }
  Exons <- length(exonStarts)
  atrack <- AnnotationTrack(start = exonStarts,
                            width = exonWidths,
                            chromosome = chr,
                            strand = rep(strand,Exons), 
                            genome = "hg38",
                            # name = paste0(chr, ": ", gene, collapse = ""),
                            group = rep(label,Exons),
                            fontsize = 20
                            )
  if (strand == "-"){ 
    feature(atrack) <- as.character(c(flank2, flank1))
  }else{
    feature(atrack) <- as.character(c(flank1, flank2))
  }
  ht <- HighlightTrack(trackList=atrack, 
                       start = IS, 
                       width = 4, 
                       chromosome = chr
                       )
  svg(paste0(gene, "_gene_zoom.svg", collapse = " "),width = 12, height = 3)
  # Code of the plot
  plotTracks(ht, 
             featureAnnotation = "feature",
             groupAnnotation  = "group",
             fontcolor.feature = 1, 
             cex.feature = 1.8,
             lwd = 1.5,
             fontsize = 20,
             cex.group = 3,
             fontcolor.group = 'black',
             lineheight = 1,
             # col.title = 'black',
             # rotation.title = 0,
             showTitle = FALSE,
             collapse = FALSE
             )
  # Close the graphics device
  dev.off()
}

myZoom('KPNA4', 2, "chr3", 2, 3, "Exon", "KPNA4", 160535998)
myZoom('ZGPAT', 10, "chr20", 2, 3, "Exon", "ZGPAT", 63719663) # why is there suddenly an error?
myZoom('Polycomb-Associated_Non-Coding_RNAs', 4, "chr19", 3, 4, "Exon", "MIDN", 1246860)
myZoom('Polycomb-Associated_Non-Coding_RNAs', 2, "chr19", 3, 4, "Exon", "lncRNA", 1246860)
myZoom('Polycomb-Associated_Non-Coding_RNAs', 2, "chr19", 3, 4, "Exon", "Polycomb-Associated_Non-Coding_RNAs", 1246860)
