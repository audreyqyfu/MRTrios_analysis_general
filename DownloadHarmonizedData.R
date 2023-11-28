

####installed packages####
BiocManager::install("sesameData", lib="/mnt/ceph/oluw5072/Rpackages/Work/")
BiocManager::install("sesame", lib="/mnt/ceph/oluw5072/Rpackages/Work/")
BiocManager::install("SummarizedExperiment", lib="/mnt/ceph/oluw5072/Rpackages/Work/")
BiocManager::install("ExperimentHub", lib="/mnt/ceph/oluw5072/Rpackages/Work/")
BiocManager::install("BiocGenerics", lib="/mnt/ceph/oluw5072/Rpackages/Work/")
BiocManager::install("AnnotationHub", lib="/mnt/ceph/oluw5072/Rpackages/Work/")
BiocManager::install("BiocFileCache", lib="/mnt/ceph/oluw5072/Rpackages/Work/")
install.packages("dbplyr", lib="/mnt/ceph/oluw5072/Rpackages/Work/")
BiocManager::install("MatrixGenerics", lib="/mnt/ceph/oluw5072/Rpackages/Work/")
BiocManager::install("GenomicRanges", lib="/mnt/ceph/oluw5072/Rpackages/Work/")
BiocManager::install("S4Vectors", lib="/mnt/ceph/oluw5072/Rpackages/Work/")
BiocManager::install("IRanges", lib="/mnt/ceph/oluw5072/Rpackages/Work/")
BiocManager::install("GenomeInfoDb", lib="/mnt/ceph/oluw5072/Rpackages/Work/")
BiocManager::install("Biobase", lib="/mnt/ceph/oluw5072/Rpackages/Work/")
BiocManager::install("splitstackshape", lib="/mnt/ceph/oluw5072/Rpackages/Work/")

####load library####
library (TCGAbiolinks, lib.loc="/mnt/ceph/oluw5072/Rpackages/Work/")
library (BiocGenerics, lib.loc="/mnt/ceph/oluw5072/Rpackages/Work/")
library (dbplyr, lib.loc="/mnt/ceph/oluw5072/Rpackages/Work/")
library (BiocFileCache, lib.loc="/mnt/ceph/oluw5072/Rpackages/Work/")
library (AnnotationHub, lib.loc="/mnt/ceph/oluw5072/Rpackages/Work/")
library (ExperimentHub, lib.loc="/mnt/ceph/oluw5072/Rpackages/Work/")
library (sesameData, lib.loc="/mnt/ceph/oluw5072/Rpackages/Work/")
library (BiocGenerics, lib.loc="/mnt/ceph/oluw5072/Rpackages/Work/")
library (MatrixGenerics, lib.loc="/mnt/ceph/oluw5072/Rpackages/Work/")
library (sesame, lib.loc="/mnt/ceph/oluw5072/Rpackages/Work/")
library (S4Vectors, lib.loc="/mnt/ceph/oluw5072/Rpackages/Work/")
library (IRanges, lib.loc="/mnt/ceph/oluw5072/Rpackages/Work/")
library (GenomeInfoDb, lib.loc="/mnt/ceph/oluw5072/Rpackages/Work/")
library (GenomicRanges, lib.loc="/mnt/ceph/oluw5072/Rpackages/Work/")
library (Biobase, lib.loc="/mnt/ceph/oluw5072/Rpackages/Work/")
library (SummarizedExperiment, lib.loc="/mnt/ceph/oluw5072/Rpackages/Work/")
library(splitstackshape, lib.loc="/mnt/ceph/oluw5072/Rpackages/Work/")

# download harmonized data
query.met <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "DNA Methylation",
  data.format = "TXT",
  platform = "Illumina Human Methylation 450"
)
GDCdownload(query.met)

# check for the dimension
dim(query.met$results[[1]])

methyl <- GDCprepare(query.met, directory = "GDCdata")

#extract the genomic data matrix
methyl.data <- assay(methyl)
methyl %>% assay %>% head %>% as.data.frame

## extract methylation probes
probes <- methyl %>% rowRanges %>% as.data.frame
probes[1:10, c(1,2,21)]

## extract only chr, start, gene, gene_HGNC
probes.small <- probes[, c(1, 2, 20, 21)]

## Check if it is data frame
is.data.frame(probes.small)

## sorting data
# sort data using the row.names function of the dataframe 
methyl.data.sorted<-methyl.data[order(row.names(methyl.data)),]

## TCGA-BRCA
# loading from fu_lab path
load("/mnt/ceph/fu_lab/TCGA/GDCdata/TCGA-BRCA/TCGA.meth.RData")
#TCGA.meth[,1:4]


#combine TCGA.meth(4 columns) with all methyl.data
probes_LUAD <- cbind(TCGA.meth[, 1:4],methyl.data.sorted) 

##Next step : save 
save(probes_LUAD,file = "/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/probes_LUAD.RData")

##load data
load("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/probes_LUAD.RData")


