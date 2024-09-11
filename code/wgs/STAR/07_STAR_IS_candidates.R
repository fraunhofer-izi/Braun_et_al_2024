library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)

# for unique mappers MAPQ = 255... with filtering for unique mappers, they are all 255 anyway
# instead of MAPQ show Alignment score "is of max possible"
# fil is the AS-NH-HI filtered chimeric.out file e.g. '/CARposCD5pos/CARposCD5pos.Chimeric.Carvykti_vector.NH1HI1AS100_filtered.txt'

## load AS-NH-HI filtered chimeric.out file
myDat <- function(fil){
  j <- read.table(fil, header = FALSE)
  j <- select(j, 1:7, 10, 12, 14, 16, 18)
  col_header <- c("donorA", 
                  "positionA", 
                  "strandA", 
                  "acceptorB", 
                  "positionB", 
                  "strandB", 
                  "junction_type", 
                  "read_name", 
                  "cigarA", 
                  "cigarB", 
                  "max_possible_AS", 
                  "AS")
  colnames(j) <- col_header
  return(j)
}

## swap column for "donor" and "acceptor" so that "donor" is always the host chromosome and remove rows with Carvykti-Carvykti chimera

myOrder <- function(dat){
  j <- dat
  for (i in 1:nrow(dat)){
    if (dat[i,1] == "Carvykti_vector"){
      j[i,1] <- dat[i,4]
      j[i,4] <- dat[i,1]
      j[i,2] <- dat[i,5]
      j[i,5] <- dat[i,2]
      j[i,3] <- dat[i,6]
      j[i,6] <- dat[i,3]
      j[i,9] <- dat[i,10]
      j[i,10] <- dat[i,9]
    }
  }
  j <- filter(j, j[,1] != "Carvykti_vector")
  return(j)
}

## filter host-Carvykti reads only

myChim <- function(dat){
  j <- filter(dat, dat$junction_type != -1) 
  return(j)
}

## filter integration site flanking reads only

myFlank <- function(dat){
  j <- filter(dat, dat$junction_type == -1) 
  return(j)
}

## restriction: insertion site candidates must have a host-Carvykti chimeric read to be nominated as such
## removes those flanking reads indicating suggesting an insertion site that is not supported by a host-Carvykti chimeric read

myChimsupported <- function(chim, flank){
  x <- chim[, 1] %>% unique()
  y <- flank[, 1] %>% unique()
  z <- intersect(x,y)
  j <- filter(flank, flank[,1] %in% z)
  return(j)
}

## merge back the filtered flanking reads with the chimeric reads

myMerge <- function(chim,fflank){
  j <- rbind(chim,fflank)
  return(j)
}

## compare suggested insertion positions within a chromosome. distance between them
## set a range of 200 bases around insertion site position that is accepted as one group

myRange <- function(chim){
  x <- chim[, 1] %>% unique()
  candidates <- data.frame(matrix(, nrow = 0, ncol = 5))
  for (i in 1:length(x)){
    positions <- vector()
    for (j in 1:nrow(chim)){
      if (chim[j,1] == x[i]){
        positions <- c(positions, chim[j,2]) 
      }
    }
    minPos <- min(positions)
    maxPos <- max(positions)
    diffPos <- maxPos - minPos
    minRange <- minPos - 200
    maxRange <- minPos +200
    app <- c(x[i], minPos, diffPos, minRange, maxRange)
    candidates <- rbind(candidates, app)
    colnames(candidates) <- c("chr", "min_pos", "diff", "min_range", "max_range")
    candidates$min_pos <- as.numeric(candidates$min_pos)
    candidates$diff <- as.numeric(candidates$diff)
    candidates$min_range <- as.numeric(candidates$min_range)
    candidates$max_range <- as.numeric(candidates$max_range)
  }
  return(candidates)
}

## gives a warning if differences between suggested insertion sites is bigger than 10. It could suggest that there are multiple insertion points within the same chromosome.

myWarning <- function(candidates){
  message <- vector()
  for (i in 1:nrow(candidates)){
    if (candidates[i,3] > 10){
      message <- c(message, candidates[i,1])
    }
  }
  print(paste0(c("Possibly more than one Insertion Site in chromosome:", message), collapse = " "))
}


## assign each read to a insertion site candidate and add it as an additional first column

myTags <- function(filt, candidates){
  tag <- data.frame(matrix(NA, nrow = nrow(filt), ncol = 1))
  for (i in 1:nrow(filt)){
    for (j in 1:nrow(candidates)){
      if ((filt[i,1] == candidates[j,1]) && (candidates[j,4] <= filt[i, 2]) && (filt[i, 2]<= candidates[j,5])){
        tag[i,1] <- paste0(c(candidates[j,1], candidates[j,2]), collapse = " ")
      }
    }
  }
  colnames(tag) <- "Integration Site candidate"
  tagged <- cbind(tag, filt)
  return(tagged)
}

# sort table and omit rows without candidate assignment (NA, flank reads unsupported by chimeric read)
mySort <- function(tagged){
  sorted <- tagged[order(tagged$`Integration Site candidate`),]
  sorted <- na.omit(sorted)
  return(sorted)
}

# count occurrences per candidate
myCount <- function(sorted){
  sorted %>%
    group_by(sorted$`Integration Site candidate`) %>%
    mutate(count = n()) %>%
    ungroup() %>%
    select(!`sorted$\`Integration Site candidate\``)
}

# filter occurrences of candidates >=2

myCandidateSupports <- function(counted){
  counted %>% filter(counted$'count' >=2)
}

# as bar-plot (add auto-saving as png later)

# myPlot <- function(finals){
#   uniq <- finals %>%
#     count(`Integration Site candidate`, sort = FALSE, name = "count")
#   ggplot(uniq, aes(x=`Integration Site candidate`, y=count)) +
#     geom_bar(stat = "identity", fill = "blue", color = "black") + 
#     ylim(0,25) +
#     labs(title = "Support per candidate", x = "Candidate", y = "Reads") +
#     theme_minimal()
# }

## as part of myPlot. alternative coding instead of count()
# uniq <- finals %>%
#   distinct(`Integration Site candidate`, .keep_all = TRUE) %>%
#   select(`Integration Site candidate`,'count')



# ########## test each function
# CARposCD5pos <- myDat('CARposCD5pos/CARposCD5pos.Chimeric.Carvykti_vector.NH1HI1AS100_filtered.txt')
# CARposCD5pos_swapped <- myOrder(CARposCD5pos)
# CARposCD5pos_chim <- myChim(CARposCD5pos_swapped)
# CARposCD5pos_flank <- myFlank(CARposCD5pos_swapped)
# CARposCD5pos_flank_filtered <- myChimsupported(CARposCD5pos_chim, CARposCD5pos_flank)
# CARposCD5pos_filtered <- myMerge(CARposCD5pos_chim, CARposCD5pos_flank_filtered)
# CARposCD5pos_candidates <- myRange(CARposCD5pos_chim)
# myWarning(CARposCD5pos_candidates)
# tagged_CARposCD5pos <- myTags(CARposCD5pos_filtered,CARposCD5pos_candidates)
# sorted_tag_CARposCD5pos <- mySort(tagged_CARposCD5pos)
# counted_CARposCD5pos <- myCount(sorted_tag_CARposCD5pos)
# above2_CARposCD5pos <- myCandidateSupports(counted_CARposCD5pos)
# myPlot(above2_CARposCD5pos)


##########  COMBINE THE FUNCTIONS INTO A FLOW ##############

myCombination <- function(STARoutput){
  dat <- myDat(STARoutput) %>%
    myOrder()
  chim <- myChim(dat)
  flank <- myFlank(dat)
  candidates <- myRange(chim)
  filter_flank <-  myChimsupported(chim, flank)
  filt <- myMerge(chim, filter_flank)
  final <- myTags(filt, candidates) %>%
    mySort() %>%
    myCount() %>%
    myCandidateSupports()
  print(myWarning(candidates))
  return(final)
}

############# DATA OUTPUT: final tables with all information and counted supporting reads

`CD5-CAR+` <- myCombination('CARposCD5dim/CARposCD5dim.Chimeric.Carvykti_vector.NH1HI1AS100_filtered.txt')
`CD5+CAR+` <- myCombination('CARposCD5pos/CARposCD5pos.Chimeric.Carvykti_vector.NH1HI1AS100_filtered.txt')
`CD5+CAR-` <- myCombination('CARnegCD5pos/CARnegCD5pos.Chimeric.Carvykti_vector.NH1HI1AS100_filtered.txt')
`Apherese` <- myCombination('only_arribaSetting_multimapper/S1.Chimeric.Carvykti_vector.NH1HI1AS100_filtered.txt')

##### output as csv to import into word document as supplemental tables

write.csv2(file = "CD5negCARpos_candidates_support.csv", x = `CD5-CAR+`)
write.csv2(file = "CD5posCARpos_candidates_support.csv", x = `CD5+CAR+`)
write.csv2(file = "CD5posCARneg_candidates_support.csv", x = `CD5+CAR-`)
write.csv2(file = "Apherese_candidates_support.csv", x = `Apherese`)

#### PLOT data in barplot with all samples
## optimize function with for loop, putting arguments in a list?

myGroupPlot <- function(sample1, sample2, sample3, sample4){
  uniq1 <- sample1 %>%
    count(`Integration Site candidate`, sort = FALSE, name = "count")
  uniq1$Source <- deparse(substitute(sample1))
  uniq2 <- sample2 %>%
    count(`Integration Site candidate`, sort = FALSE, name = "count")
  uniq2$Source <- deparse(substitute(sample2))
  uniq3 <- sample3 %>%
    count(`Integration Site candidate`, sort = FALSE, name = "count")
  uniq3$Source <- deparse(substitute(sample3))
  uniq4 <- sample4 %>%
    count(`Integration Site candidate`, sort = FALSE, name = "count")
  uniq4$Source <- deparse(substitute(sample4))
  graph <- rbind(uniq1, uniq2, uniq3, uniq4)
  graph <- cbind((str_remove(graph$`Integration Site candidate`, " .*")),select(graph, 2:3))
  colnames(graph) <- c("Integration Site candidate", "count", "source")
  graph <- graph %>%
    complete(`Integration Site candidate`, source, fill = list(count= 0))
  ggplot(graph, aes(x=`Integration Site candidate`, y=count, fill = source)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("CD5-CAR+" = "#6699CC", "CD5+CAR+" = "#DDAA33")) +
    labs(fill = "Sample") +
    ylim(0,25) +
    labs(x = "Candidate", y = "Reads") +
    theme_minimal()
}

myGroupPlot(`CD5-CAR+`, `CD5+CAR+`, `CD5+CAR-`, `Apherese`)
