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




