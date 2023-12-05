##.....Formatting Patient ID.....##
library(MRGN,lib.loc="/mnt/ceph/oluw5072/Rpackages/MRGN_R/")
library(MRTrios,lib.loc="/mnt/ceph/oluw5072/Rpackages/MRGN_R/")
#library(MRPC,lib.loc="/mnt/ceph/oluw5072/Rpackages/MRGN_R/")
library(data.table,lib.loc="/mnt/ceph/oluw5072/Rpackages/MRGN_R/")
library(na.tools,lib.loc="/mnt/ceph/oluw5072/Rpackages/MRGN_R/")

#Load the inital Methylation dataset
LUAD.meth<- as.data.frame(fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/split.names.LUAD.meth.logit.txt"))
dim(LUAD.meth)

#Gene Expression dataset
LUAD.gene<- fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/LUAD_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt")
dim(LUAD.gene)

#CNA dataset
LUAD.cna<- fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/LUAD_tcga_pan_can_atlas_2018/data_cna.txt")
dim(LUAD.cna)

#split ID name for METH
for(i in 5:ncol(LUAD.meth)){
  
  print(i)
  
  #split the name in the TCGA data
  split.ID <- unlist(strsplit(as.character(colnames(LUAD.meth)[i]), '-'))
  split.ID
  
  
  #now merge the first 4 parts
  new.rowname <- paste(split.ID[1],split.ID[2],split.ID[3],  sep="-")
  new.rowname
  colnames(LUAD.meth)[i] <- new.rowname
} 
#save it to a new text file
write.table(LUAD.meth, file = "/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/split.names.LUAD_meth.new.txt", sep = "\t", row.names = FALSE,
            col.names = TRUE,quote=FALSE)

#split ID name for cna data set
for(i in 3:ncol(LUAD.cna)){
  
  print(i)
  
  #split the name in the TCGA data
  split.ID <- unlist(strsplit(as.character(colnames(LUAD.cna)[i]), '-'))
  split.ID
  
  
  #now merge the first 4 parts
  new.rowname <- paste(split.ID[1],split.ID[2],split.ID[3],  sep="-")
  new.rowname
  colnames(LUAD.cna)[i] <- new.rowname
}
#save it to a new text file
write.table(LUAD.cna, file = "/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/split.names.LUAD_cna.new.txt", sep = "\t", row.names = FALSE,
            col.names = TRUE,quote=FALSE)

#split ID name for gene data set
for(i in 3:ncol(LUAD.gene)){
  
  print(i)
  
  #split the name in the TCGA data
  split.ID <- unlist(strsplit(as.character(colnames(LUAD.cna)[i]), '-'))
  split.ID
  
  
  #now merge the first 4 parts
  new.rowname <- paste(split.ID[1],split.ID[2],split.ID[3],  sep="-")
  new.rowname
  colnames(LUAD.gene)[i] <- new.rowname
}
#save it to a new text file
write.table(LUAD.gene, file = "/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/split.names.LUAD_gene.new.txt", sep = "\t", row.names = FALSE,
            col.names = TRUE,quote=FALSE)

#Clinical dataset

LUAD.clinical<- fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/LUAD_tcga_pan_can_atlas_2018/data_clinical_patient.txt")

LUAD.cdata<-LUAD.clinical[-(1:4),]



#save it to a new text file
write.table(LUAD.cdata, file = "/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/split.names.LUAD_clinical.new.txt", sep = "\t", row.names = FALSE,
            col.names = TRUE, quote=FALSE)

#Load the new Methylation dataset
LUAD.meth<- as.data.frame(fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/split.names.LUAD_meth.new.txt"))
dim(LUAD.meth)

#new Gene Expression dataset
LUAD.gene<- fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/split.names.LUAD_gene.new.txt")
dim(LUAD.gene)

#new CNA dataset
LUAD.cna<- fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/split.names.LUAD_cna.new.txt")
dim(LUAD.cna)

#new clinical dataset
clinical.LUAD<-fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/split.names.LUAD_clinical.new.txt")
dim(clinical.LUAD)
