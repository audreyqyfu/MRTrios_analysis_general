##.....Main find PCs.....##
library(MRGN,lib.loc="/mnt/ceph/fern5249/Rpackages/")
library(MRTrios,lib.loc="/mnt/ceph/fern5249/Rpackages/")
#library(MRPC,lib.loc="/mnt/ceph/oluw5072/Rpackages/MRGN_R/")
library(data.table,lib.loc="/mnt/ceph/oluw5072/Rpackages/MRGN_R/")
library(na.tools,lib.loc="/mnt/ceph/oluw5072/Rpackages/MRGN_R/")


#Load the Methylation dataset
LUAD.meth<- as.data.frame(fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/split.names.LUAD.meth.logit.txt"))
dim(LUAD.meth)

#Gene Expression dataset
LUAD.gene<- fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/luad_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt")
dim(LUAD.gene)

#CNA dataset
LUAD.cna<- fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/luad_tcga_pan_can_atlas_2018/data_cna.txt")
dim(LUAD.cna)

#clinical dataset

LUAD.clinical<- fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/luad_tcga_pan_can_atlas_2018/data_clinical_patient.txt")
dim(LUAD.clinical)

LUAD.cdata<-LUAD.clinical[-(1:4),]

for (i in 1:nrow(LUAD.cdata)) {
  LUAD.cdata[i, 1] <- paste(LUAD.cdata[i, 1], "-01", sep = "")
}

#save it to a new text file
write.table(LUAD.cdata, file = "mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/new_data_clinical_patient.txt", sep = "\t", row.names = FALSE,
            col.names = TRUE,quote=FALSE)

clinical.LUAD<-fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/new_data_clinical_patient.txt")
dim(clinical.LUAD)


#finding common individuals among the 3 datasets
com.ind = intersect(colnames(LUAD.gene)[3:ncol(LUAD.gene)], colnames(LUAD.meth)[5:ncol(LUAD.meth)])

#Use the function to get PC score matrix and significantly associated PCs.
pc.gene = findPCsGeneral(as.data.frame(LUAD.gene), 3, 1, com.ind, 100)
write.table(pc.gene[[1]], file = paste("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/PCA.gene.exp.txt", sep = ""), sep = "\t", row.names = TRUE,
            col.names = NA,  quote=FALSE)
saveRDS(pc.gene[[2]], file="/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/gene.exp.sig.asso.pcs.RData")
#load(file="/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/gene.exp.sig.asso.pcs.RData")

pc.meth = findPCsGeneral(as.data.frame(LUAD.meth), 5, 2, com.ind, 100)
write.table(pc.meth[[1]], file = paste("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/PCA.meth.txt", sep = ""), sep = "\t", row.names = TRUE,
            col.names = NA, quote=FALSE)
saveRDS(pc.meth[[2]], file="/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/meth.sig.asso.pcs.RData")
    

#use the function to get the indices matrix 
gene.table = findIndex(LUAD.gene, 3, 1, com.ind, com.ind)
write.table(gene.table, file = paste("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/gene.exp.table.txt", sep = ""), sep = "\t", row.names = TRUE,
            col.names = NA,  quote=FALSE)

meth.table = findIndex(LUAD.meth, 5, 2, com.ind, com.ind)
write.table(meth.table, file = paste("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/meth.table.txt", sep = ""), sep = "\t", row.names = TRUE,
            col.names = NA,  quote=FALSE)




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
 


final.result = analyzeTrios(LUAD.meth, LUAD.gene, LUAD.cna, trios[1:5,], pc.meth, pc.gene, meth.sig.asso.pcs[[1]], gene.sig.asso.pcs[[1]],clinical.LUAD, meth.table, gene.table,5,26,5,3)
write.table(final.result, file = paste("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/analyze.trios.LUAD.txt", sep = ""), sep = "\t", row.names = TRUE,
            col.names = TRUE,  quote=FALSE)
LUAD.analyzeTrios <- read.table("mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/analyze.trios.LUAD.txt", header=TRUE)


#The PC matrix
dim (pc.gene[[1]])
pc.gene[[1]]

dim (pc.meth[[1]])
pc.meth[[1]]

#List of significantly associated PCs
pc.gene[[2]]$sig.asso.covs
pc.meth[[2]]$sig.asso.covs

write.table(pc.meth[[1]], file = paste("/mnt/ceph/oluw5072/TGCA-LUAD/Analysis/PCA.meth.txt, sep = ""), sep = "\t", row.names = TRUE,
            col.names = NA, append = TRUE, quote=FALSE)

saveRDS(pc.meth[[2]], file="/mnt/ceph/oluw5072/TCGA-LUAD/Analysis/meth.sig.asso.pcs.RData")

write.table(meth.table, file = paste("/mnt/ceph/oluw5072/TGCA-LUAD/Analysis/meth.table.txt, sep = ""), sep = "\t", row.names = TRUE,
            col.names = NA, append = TRUE, quote=FALSE)



write.table(pc.gene[[1]], file = paste("/mnt/ceph/oluw5072/TGCA-LUAD/Analysis/PCA.gene.exp.txt, sep = ""), sep = "\t", row.names = TRUE,
            col.names = NA, append = TRUE, quote=FALSE)

saveRDS(pc.gene[[2]], file="/mnt/ceph/oluw5072/TCGA-LUAD/Analysis/gene.exp.sig.asso.pcs.RData")

write.table(gene.table, file = paste("/mnt/ceph/oluw5072/TGCA-LUAD/Analysis/gene.exp.table.txt, sep = ""), sep = "\t", row.names = TRUE,
            col.names = NA, append = TRUE, quote=FALSE)
