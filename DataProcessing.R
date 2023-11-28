##.....Data Processing.....##

#load the library
library(data.table ,lib.loc="/mnt/ceph/oluw5072/Rpackages/MRGN_R")
library(splitstackshape, lib.loc="/mnt/ceph/oluw5072/Rpackages/MRGN_R")

#load methylation data of LUAD
load("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/probes_LUAD.RData")
dim(probes_LUAD)

#split ID name 
for(i in 5:ncol(probes_LUAD)){
  
  print(i)
  
  #split the name in the TCGA data
  split.ID <- unlist(strsplit(as.character(colnames(probes_LUAD)[i]), '-'))
  split.ID
  
  #now merge the first 4 parts
  new.rowname <- paste(split.ID[1],"-",split.ID[2],"-",split.ID[3], "-", substr (split.ID[4], 1,2), sep="")
  new.rowname
  
  #replace the old column name to the new column name
  colnames(probes_LUAD)[i] <- new.rowname
}

#save it to a new text file
write.table(probes_LUAD, file = "/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/split.names.probes_LUAD.txt", sep = "\t", row.names = FALSE,
            col.names = TRUE,quote=FALSE)


##logit transformation of methylation data##

LUAD.meth <- as.data.frame(fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/split.names.probes_LUAD.txt"))
dim(LUAD.meth)

data <- LUAD.meth[rowSums(is.na(LUAD.meth[,-(1:4)])) != ncol(LUAD.meth[,-(1:4)]), ]
dim(data)

# split data into numbers and probe info
data.info <- data[,1:4]

data.nona <- log(data[,-(1:4)]/(1-data[,-(1:4)]))

final <- cbind(data.info, data.nona)

#save it to a new text file
write.table(final, file = "/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/split.names.LUAD.meth.logit.txt", sep = "\t", row.names = FALSE,
            col.names = TRUE,quote=FALSE)
