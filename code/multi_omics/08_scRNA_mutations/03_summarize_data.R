# Select Files
files =list.files(path="./data/scRNA_BAMs_sorted/split_bams/",pattern="*.txt",recursive = T)

# Create data.frame
sc_mut_data=data.frame(matrix(ncol=6,nrow=0))

# Read pileup output copy data of TET2 mutation into data.frame
for (i in 1:length(files)){
  data = read.table(paste0("/mnt/ribolution/user_worktmp/markus.kreuz/Projects/2024_07_16_Fusions_Köln/scRNA_BAMs_sorted/split_bams/",files[[i]]),sep="\t")
  names(data)= c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SUMMARY")
  ref=data$REF[1]
  alt=strsplit(data$ALT[1],",")[[1]]
  s = strsplit(data$SUMMARY[[1]],":")[[1]][2]
  s = as.numeric(strsplit(s,",")[[1]])
  n_total=sum(s)
  n_alt=ifelse("T" %in% alt,s[which(alt=="T")+1],0)
  sample = strsplit(files[[i]],"/")[[1]][1]
  cells = strsplit(files[[i]],"_")[[1]]
  cells = cells[-c(1,length(cells))]
  cells = paste(cells,collapse = "_")

  sc_mut_data=rbind(sc_mut_data,
  c(sample,cells,ref,paste(alt,collapse=","),n_total,n_alt)
  )
}

# Set colnames and change allele counts to numeric
names=c("sample","cells","ref","mut","n_total","n_mut")
colnames(sc_mut_data)=names
sc_mut_data$n_total=as.numeric(sc_mut_data$n_total)
sc_mut_data$n_mut=as.numeric(sc_mut_data$n_mut)

# Save results
saveRDS(sc_mut_data,file = "./data/sc_mutation_data.RDS")

