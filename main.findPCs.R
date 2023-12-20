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



#splitting ID function
IDsplit <- function(data, delimiter = '-', start_columns = NULL) {
  
  for (col_index in start_columns) {
    
    split_names <- unlist(strsplit(colnames(data)[col_index], delimiter))
    
    # Ensure at least four parts exist before merging
    if (length(split_names) >= 4) {
      new_rowname <- paste(split_names[1:3], collapse = delimiter)
      colnames(data)[col_index] <- new_rowname
    }
  }
  return(data)
}


#split ID name for METH, start splitting from column 5
LUAD.meth_1 <- IDsplit(LUAD.meth, delimiter = '-', start_columns = c(5:ncol(LUAD.meth)))
#save it to a new text file
write.table(LUAD.meth_1, file = "/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/split.names.LUAD_meth.new.txt", sep = "\t", row.names = FALSE,
            col.names = TRUE,quote=FALSE)

#split ID name for cna data seT, start splitting from column 3
LUAD.cna_1 <- IDsplit(LUAD.cna, delimiter = '-', start_columns = c(3:ncol(LUAD.cna)))

#save it to a new text file
write.table(LUAD.cna_1, file = "/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/split.names.LUAD_cna.new.txt", sep = "\t", row.names = FALSE,
            col.names = TRUE,quote=FALSE)

#split ID name for gene data set, start splitting from column 3
LUAD.gene_1 <- IDsplit(LUAD.gene, delimiter = '-', start_columns = c(3:ncol(LUAD.gene)))

#save it to a new text file
write.table(LUAD.gene_1, file = "/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/split.names.LUAD_gene.new.txt", sep = "\t", row.names = FALSE,
            col.names = TRUE,quote=FALSE)



#Load the Methylation dataset
LUAD.meth<- as.data.frame(fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/split.names.LUAD_meth.new.txt"))
dim(LUAD.meth)

#Gene Expression dataset
LUAD.gene<- fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/split.names.LUAD_gene.new.txt")
dim(LUAD.gene)

#CNA dataset
LUAD.cna<- fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/split.names.LUAD_cna.new.txt")
dim(LUAD.cna)

#clinical dataset
clinical.LUAD<-fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/split.names.LUAD_clinical.new.txt")
dim(clinical.LUAD)

#clinical dataset
#read in the clinical dataset
#formatting clinical to match meth and gene ID's



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




