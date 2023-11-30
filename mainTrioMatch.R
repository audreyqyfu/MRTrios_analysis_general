##.....MainTrioMatch.....##

#load the library
library(MRGN,lib.loc="/mnt/ceph/oluw5072/Rpackages/MRGN_R/")
library(MRTrios,lib.loc="/mnt/ceph/oluw5072/Rpackages/MRGN_R/")
library(data.table,lib.loc="/mnt/ceph/oluw5072/Rpackages/MRGN_R/")
library(splitstackshape,lib.loc="/mnt/ceph/oluw5072/Rpackages/MRGN_R/")
library(org.Hs.eg.db,lib.loc="/mnt/ceph/oluw5072/Rpackages/MRGN_R/")


#read in the datasets

#Load the Methylation dataset
LUAD.meth<- as.data.frame(fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/split.names.LUAD.meth.logit.txt"))
#Gene Expression dataset
LUAD.gene<- fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/luad_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt")
dim(LUAD.gene)
#CNA dataset
LUAD.cna<- fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/luad_tcga_pan_can_atlas_2018/data_cna.txt")
dim(LUAD.cna)



################################### first round of matching trios ##############

#go to the findTrioAll function and run the code and get the trios dataset
final<- findTrioAll(LUAD.meth,LUAD.cna, LUAD.gene, nStartMeth=5, nStartGene=3, nStartCNA=3)
print(final[1:5,])

colnames(final) <- c("Gene name", "meth.row", "cna.row", "gene.row")
final.colnames <- colnames(final)


#save the trios dataset to the path in txt format
#some methylation probes in this dataset are not matched to any gene
write.table(t(final.colnames), file = "/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/pre.final.txt", sep = "\t", row.names = FALSE,
            col.names = TRUE,append = FALSE,quote=FALSE)

write.table(final, file = "/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/pre.final.txt", sep = "\t", row.names = FALSE,
            col.names = TRUE,quote=FALSE)
#........................................................................................................................................#



################################### second round of matching trios ##############

#find the duplicates for genes with multiple entrez ids in CNA data
final.tmp <- addDupsCNA(final, cna)

if (!is.null (final.tmp[[2]])) {
  final <- final.tmp[[1]]
  print(dim(final))
  
  dup.final <- final.tmp[[2]]
  print(dim(dup.final))
  
  #find the duplicates for genes with multiple entrez ids in Gene Exp data
  final.res <- addDupsGENE(final, dup.final, gene)
  print(final.res[1:10,])
  print(dim(final.res))
} else {
  final.res <- final
}

#save the final trio data matrix to text file
write.table(t(colnames(final.res)), file = "/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/Trios.final.txt", sep = "\t", row.names = FALSE,
            col.names = TRUE,quote=FALSE)

write.table(final.res, file = "/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/Trios.final.txt", sep = "\t", row.names = FALSE,
            col.names = TRUE,quote=FALSE)
