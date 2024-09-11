Sys.setenv( TZ="Etc/GMT+1" )

setwd('~/gviz/IntegrationSite_coverage/') #this can be changed

library(Gviz)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(stringr)

## create object: chromosome ideogram; can a gene be displayed instead of a chromosome? does not seem like it
### current genomic location will be indicated on the chromosome by a red box/line

gen = 'hg38'
chr01 = 'chr3'
chr02 = 'chr19'
chr03 = 'chr20'

IS01 <- 160536002
IS02 <- 1246860 
IS03 <- 63719665 

myIdeo <- function(chr, gene, IS){
  itrack <- IdeogramTrack(genome = gen,
                          chromosome = chr,
                          showBandId = TRUE,
                          cex.bands = 1.5,
                          fontcolor = "black",
                          cex = 1
                          )
  svg(paste0(c(chr, "_", gene, "_Ideogram.svg"),
             collapse = ""),
      width = 10,
      height = 2)
  plotTracks(itrack,
             chromosome = chr,
             from = IS,
             to = IS+4,
             fontsize = 40,
             main = gene,
             cex.main = 4
             )
  dev.off()
}

myIdeo(chr01, "KPNA4", IS01)
myIdeo(chr02, "Polycomb-Associated non-coding RNAs", IS02)
myIdeo(chr03, "ZGPAT", IS03)
