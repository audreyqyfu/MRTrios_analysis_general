##.....Main analyzeTrios.....##
#To install MRTrios package from github
library(usethis,lib.loc="/mnt/ceph/oluw5072/Rpackages/MRGN_R/")
library(devtools,lib.loc="/mnt/ceph/oluw5072/Rpackages/MRGN_R/")
install_github("audreyqyfu/MRTrios",lib="/mnt/ceph/oluw5072/Rpackages/MRGN_R/")

#Load library's
library(MRGN,lib.loc="/mnt/ceph/oluw5072/Rpackages/MRGN_R/")
library(MRTrios,lib.loc="/mnt/ceph/oluw5072/Rpackages/MRGN_R/")
#library(MRPC,lib.loc="/mnt/ceph/oluw5072/Rpackages/MRGN_R/")
library(data.table,lib.loc="/mnt/ceph/oluw5072/Rpackages/MRGN_R/")
library(na.tools,lib.loc="/mnt/ceph/oluw5072/Rpackages/MRGN_R/")


#Load the new Methylation dataset
LUAD.meth<- as.data.frame(fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/split.names.LUAD_meth.new.txt"))
dim(LUAD.meth)

#new Gene Expression dataset
LUAD.gene<- fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/split.names.LUAD_gene.new.txt")
dim(LUAD.gene)

#new CNA dataset
LUAD.cna<- fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/split.names.LUAD_cna.new.txt")
dim(LUAD.cna)

#Clinical dataset
LUAD.clinical<- fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/LUAD_tcga_pan_can_atlas_2018/data_clinical_patient.txt")
#Remove the first 4 rows that are not needed for the analysis
LUAD.cdata<-LUAD.clinical[-(1:4),]

#save it to a new text file
write.table(LUAD.cdata, file = "/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/split.names.LUAD_clinical.new.txt", sep = "\t", row.names = FALSE,
            col.names = TRUE, quote=FALSE)
#load the new clinical dataset
clinical.LUAD<-fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/split.names.LUAD_clinical.new.txt")
dim(clinical.LUAD)



#Read in the PC score matrix
pc.meth<- read.table("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/PCA.meth.txt", header = TRUE)

pc.gene<- read.table("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/PCA.gene.exp.txt", header=TRUE)


#reading in the Trios data
trios <- data.frame(fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/trio.final.protein.coding.txt"))

#read in the indices table
meth.table<- fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/meth.table.txt", drop = 1)
gene.table<- fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/gene.exp.table.txt", drop = 1)


#read in the sig pcs data
meth.sig.asso.pcs<- readRDS("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/meth.sig.asso.pcs.RData")
gene.sig.asso.pcs<- readRDS("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/gene.exp.sig.asso.pcs.RData")
 


final.result = analyzeTrios(LUAD.meth, LUAD.gene, LUAD.cna, trios[1:100000,], pc.meth, pc.gene, meth.sig.asso.pcs[[1]], gene.sig.asso.pcs[[1]],clinical.LUAD, meth.table, gene.table,5,26, writeToFile =TRUE, file= "/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/analyze.trios.LUAD_1.txt")

##Write to file
write.table(final.result, file = paste("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/analyze.trios.LUAD_1.txt", sep = ""), sep = "\t", row.names = TRUE,
            col.names = TRUE,  quote=FALSE)

##read in the analyzeTrios 
#LUAD.analyzeTrios <- read.table("mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/analyze.trios.LUAD.txt", header=TRUE)
