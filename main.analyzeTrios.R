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


#Load the Methylation dataset
LUAD.meth<- as.data.frame(fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/split.names.LUAD.meth.logit.txt"))
dim(LUAD.meth)

# Apply partial matching to column names #I am using this
colnames(LUAD.meth)[5:ncol(LUAD.meth)] <- sapply(strsplit(colnames(LUAD.meth)[5:ncol(LUAD.meth)], "-"), function(parts) paste(parts[1:3], collapse="-"))
LUAD.meth[1, 1:10]


#Gene Expression dataset
LUAD.gene<- fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/luad_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt")
dim(LUAD.gene)

# Apply partial matching to column names
colnames(LUAD.gene)[3:ncol(LUAD.gene)] <- sapply(strsplit(colnames(LUAD.gene)[3:ncol(LUAD.gene)], "-"), function(parts) paste(parts[1:3], collapse="-"))
LUAD.gene[1, 1:10]

#CNA dataset
LUAD.cna<- fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/luad_tcga_pan_can_atlas_2018/data_cna.txt")
dim(LUAD.cna)

# Apply partial matching to column names
colnames(LUAD.cna)[3:ncol(LUAD.cna)] <- sapply(strsplit(colnames(LUAD.cna)[3:ncol(LUAD.cna)], "-"), function(parts) paste(parts[1:3], collapse="-"))
LUAD.cna[1, 1:10]

 #Clinical Dataset                                              
LUAD.clinical<- fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/luad_tcga_pan_can_atlas_2018/data_clinical_patient.txt")
#Remove the first 4 rows that are not needed for the analysis
LUAD.cdata<-LUAD.clinical[-(1:4),]
#save it to a new text file
write.table(LUAD.cdata, file = "/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/split.names.LUAD_clinical.new.txt", sep = "\t", row.names = FALSE,
            col.names = TRUE, quote=FALSE)
#load the new clinical dataset
clinical.LUAD<-fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/split.names.LUAD_clinical.new.txt")
dim(clinical.LUAD)

# Define removeDupsPCs function
removeDupsPCs <- function(pc, duplicated_rows) {
  # Extract the first three elements from the row names 
  new_row_names <- sapply(strsplit(rownames(pc), "-"), function(parts) paste(parts[1:3], collapse = "-"))
  # Identify duplicated row names
  duplicated_rows <- duplicated(new_row_names) | duplicated(new_row_names, fromLast = TRUE)
  # Get the row names of pc that are duplicated
  duplicated_row_names <- new_row_names[duplicated_rows]
  # Remove the duplicated rows
  pc <- pc[!duplicated_rows,]
  # Assign the modified row names to the unique data frame
  rownames(pc) <- new_row_names[!duplicated_rows]
  return(list(pc = pc, duplicated_row_names = duplicated_row_names))
}

# Read in PC score matrix
pc.meth <- read.table("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/PCA.meth.txt", header = TRUE)
pc.gene <- read.table("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/PCA.gene.exp.txt", header = TRUE)

# Remove duplicated rows from pc.meth
result_meth <- removeDupsPCs(pc.meth, NULL)
pc.meth <- result_meth$pc
duplicated_row_names_meth <- result_meth$duplicated_row_names

# Remove duplicated rows from pc.gene
result_gene <- removeDupsPCs(pc.gene, NULL)
pc.gene <- result_gene$pc
duplicated_row_names_gene <- result_gene$duplicated_row_names

# Remove corresponding columns from LUAD.cna
LUAD.cna <- LUAD.cna[, !colnames(LUAD.cna) %in% duplicated_row_names_meth, with = FALSE]

# Remove corresponding columns from LUAD.gene
LUAD.gene <- LUAD.gene[, !colnames(LUAD.gene) %in% duplicated_row_names_gene, with = FALSE]



                             
#reading in the Trios data
trios <- data.frame(fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/trio.final.protein.coding.txt"))

#read in the indices table
meth.table<- fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/meth.table.txt", drop = 1)
gene.table<- fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/gene.exp.table.txt", drop = 1)


#read in the sig pcs data
meth.sig.asso.pcs<- readRDS("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/meth.sig.asso.pcs.RData")
gene.sig.asso.pcs<- readRDS("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/gene.exp.sig.asso.pcs.RData")

final.result = analyzeTrios(LUAD.meth, LUAD.gene, LUAD.cna, trios[1:100000,], pc.meth, pc.gene, meth.sig.asso.pcs[[1]], gene.sig.asso.pcs[[1]],clinical.LUAD, meth.table, gene.table,age.col=5, race.col=26, sex.col=6, writeToFile =TRUE, file= "/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/analyze.trios.LUAD_4.txt")


