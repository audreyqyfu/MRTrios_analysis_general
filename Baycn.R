library(data.table,lib.loc="/mnt/ceph/fern5249/Rpackages")
library(na.tools,lib.loc="/mnt/ceph/fern5249/Rpackages")
library(MRGN,lib.loc="/mnt/ceph/fern5249/Rpackages")
library(MRTrios,lib.loc="/mnt/ceph/fern5249/Rpackages")
library(baycn,lib.loc="/mnt/ceph/fern5249/Rpackages")

#read in the original datasets
#Gene Expression dataset
BLCA.gene<- fread("/mnt/ceph/fern5249/GDCdata/TCGA-BLCA/Analysis/split.names.BLCA.gene.new.test.txt")

#CNA dataset
BLCA.cna<- fread("/mnt/ceph/fern5249/GDCdata/TCGA-BLCA/Analysis/split.names.BLCA.cna.new.test.txt")


#Methylation data set
BLCA.meth<- as.data.frame(fread("/mnt/ceph/fern5249/GDCdata/TCGA-BLCA/Analysis/split.names.BLCA.meth.new.test.txt"))


#clinical data
clinical.BLCA<-fread("/mnt/ceph/fern5249/GDCdata/TCGA-BLCA/Analysis/split_new_data_clinical_patient.txt")


#Read in the PC score matrix
pc.meth<- read.table("/mnt/ceph/fern5249/GDCdata/TCGA-BLCA/Analysis/PCA.meth.new.txt", header = TRUE)

pc.gene<- read.table("/mnt/ceph/fern5249/GDCdata/TCGA-BLCA/Analysis/PCA.gene.exp.new.txt",header=TRUE)


#reading in the Trios data
trios <- data.frame(fread("/mnt/ceph/fern5249/GDCdata/TCGA-BLCA/Analysis/trio.final.protein.coding.txt"))

#read in the indices table
meth.table<- fread("/mnt/ceph/fern5249/GDCdata/TCGA-BLCA/Analysis/meth.table.new.txt", drop = 1)
gene.table<- fread("/mnt/ceph/fern5249/GDCdata/TCGA-BLCA/Analysis/gene.exp.table.new.txt", drop = 1)


#read in the sig pcs data
meth.sig.asso.pcs<- readRDS("/mnt/ceph/fern5249/GDCdata/TCGA-BLCA/Analysis/meth.sig.asso.pcs.new.RData")
gene.sig.asso.pcs<- readRDS("/mnt/ceph/fern5249/GDCdata/TCGA-BLCA/Analysis/gene.exp.sig.asso.pcs.new.RData")


datamatrix<- function(TCGA.meth, gene.exp, cna, trios, pc.meth, pc.gene, meth.sig.asso.pcs, gene.sig.asso.pcs, clinical, meth.table, gene.table, age.col=5, race.col=26, sex.col=6, nObs = 30, nPCs = 50, writeToFile = FALSE, file){
  # find the common individuals between the 3 datasets
  # pc matrix has common individuals from meth and gene exp
  com.ind <- intersect(rownames(pc.meth),colnames(cna))
  
  #find the rows in clinical data for the common individuals
  rows.clinical <- match(com.ind, unlist(clinical[,1]))
  
  #extract the age, race and sex for those individuals
  age <- as.data.frame(clinical)[rows.clinical,age.col]
  race <- as.data.frame(clinical)[rows.clinical,race.col]
  sex <- as.data.frame(clinical)[rows.clinical,sex.col]
  
  #find the rows for the common individuals in the resp datasets
  ind.col.cna <- match(com.ind, colnames(cna))
  ind.col.gene <- match(com.ind, colnames(gene.exp))
  ind.col.meth <- match(com.ind, colnames(TCGA.meth))
  
  # extract data for the common individuals
  cna.common <- t(as.data.frame(cna)[,ind.col.cna])
  gene.exp.common <- t(as.data.frame(gene.exp)[,ind.col.gene])
  meth.common <- t(as.data.frame(TCGA.meth)[,ind.col.meth])
  
  #find the row numbers for the common individuals in the pc score matrix
  com.ind.pc <- match(com.ind, rownames(pc.meth))
  
  #initialize result
  result <- NULL
  
  # write to file if writeToFile is TRUE
  #if (writeToFile) {
  #result.colnames <- c("V1","V2","V3","PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "Age", "Race", "Sex")
  #write.table(t(result.colnames), file = file, sep = "\t", row.names = FALSE, col.names = FALSE, append = FALSE, quote = FALSE)
  #}
  
  #begin the loop for rows in trios
  for(i in 1:nrow(trios)){
    
    #check if the values for cna and gene exp are NA or not
    if(!(is.na(trios[i,3]))  & !(is.na(trios[i,4]))) {
      
      #print the row number
      print(i)
      
      #extract data for each dataset and create the trio
      trio.cna = cna.common[,as.numeric(trios[i,3])]
      trio.gene = gene.exp.common[,as.numeric(trios[i,4])]
      trio.meth = meth.common[,as.numeric(trios[i,2])]
      
      #combine the matrix
      trio.mat = as.data.frame(cbind(as.numeric(trio.cna), as.numeric(trio.gene), as.numeric(trio.meth)))
      
      # check whether there are enough observations 
      na.counts <- rowSums(is.na (trio.mat))
      
      # since a row with NAs is ignored in regression,
      # need to make sure that there are enough observations left
      if(sum (na.counts == 0) > nObs){
        
        #extract the updated row number from the indices table
        row.sig.pcs.meth <- unlist(meth.table[as.numeric(trios[i,2]),2])
        row.sig.pcs.gene <- unlist(gene.table[as.numeric(trios[i,4]),2])
        
        #extract the column numbers for sig pcs
        sig.pcs.cols.meth <- as.integer(unlist(meth.sig.asso.pcs[row.sig.pcs.meth]))
        sig.pcs.cols.gene <- as.integer(unlist(gene.sig.asso.pcs[row.sig.pcs.gene]))
        
        #extract the sig columns from the pc matrix with the common individuals
        # also limit the PCs to those less than nPCs
        sig.pc.gene <- pc.gene[com.ind.pc, sig.pcs.cols.gene[which (sig.pcs.cols.gene < nPCs)]]
        sig.pc.meth <- pc.meth[com.ind.pc, sig.pcs.cols.meth[which (sig.pcs.cols.meth < nPCs)]]
        
        
        #create matrix with the trios and the confounding variables
        final.mat <- cbind(trio.mat, sig.pc.gene, sig.pc.meth, age, race, sex)
        #extracting row with no NAS
        final.mat.1 <- rowSums(is.na(final.mat))
        complete_rows <- final.mat.1 == 0
        if(sum(complete_rows) > nObs) {
          final.mat.complete <- final.mat[complete_rows,]
          
          #check if a categorical variable has at least 2 levels
          if ((nlevels(as.factor(final.mat.complete$race)) >= 2)&
              (nlevels(as.factor(final.mat.complete[,1])) >=2)&
              (nlevels(as.factor(final.mat.complete$sex)) >= 2)) {
            
            # write to file if writeToFile is TRUE
            if (writeToFile) {
              write.table(final.mat.complete, file = file, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)
            }
          }
          else {
            # Remove the "race" column if it doesn't have at least 2 levels
            if (nlevels(as.factor(final.mat.complete$race)) < 2){
              final.mat.complete <- final.mat.complete[, -which(names(final.mat.complete) == "race")]
              
              
              # write to file if writeToFile is TRUE
              if (writeToFile) {
                write.table(final.mat.complete,file = file, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)
              }
            }
          }
          
          result <-  final.mat.complete
          
        }
      }
    }
  }
  #return the dataset
  colnames(result)[1] <- "V1"
  return(result)
  
}
data_BLCA=datamatrix(BLCA.meth, BLCA.gene, BLCA.cna, trios[1,], pc.meth, pc.gene, meth.sig.asso.pcs[[1]], gene.sig.asso.pcs[[1]],clinical.BLCA, meth.table, gene.table,age.col=5, race.col=26,sex.col=6)


#Remove race  column
t1 <- data_BLCA[, -which(names(data_BLCA) %in% c("race"))]
t1$sex <- ifelse(t1$sex == "Male", 0, 1)
t1[1:5,]
#unlist the data matrix data_BLCA
t2<-unlist(t1)
t2[1:10]
t3<-matrix(t2,byrow=FALSE,nrow=nrow(t1))
colnames(t3)<-colnames(t1)
#Move confounders up so they can be treated as genetic variants
t4 <- t3[,c(1,4:12, 2, 3)]


#For the 1st Trio
#Create an adjacency matrix which matches with the data_BLCA matrix
am_m1_BLCA<- matrix(c(0, 0,0,0,0,0,0,0,0,0,1,1,
                        0, 0,0,0,0,0,0,0,0,0,1,1,
                        0, 0,0,0,0,0,0,0,0,0,1,1,
                        0, 0,0,0,0,0,0,0,0,0,1,1,
                        0, 0,0,0,0,0,0,0,0,0,1,1,
                        0, 0,0,0,0,0,0,0,0,0,1,1,
                        0, 0,0,0,0,0,0,0,0,0,1,1,
                        0, 0,0,0,0,0,0,0,0,0,1,1,
                        0, 0,0,0,0,0,0,0,0,0,1,1,
                        0, 0,0,0,0,0,0,0,0,0,1,1,
                        0, 0,0,0,0,0,0,0,0,0,0,1,
                        0, 0,0,0,0,0,0,0,0,0,1,0),
                      byrow =TRUE,
                      nrow = 12)

mh_m1_BLCA<- mhEdge(data=t4,
                     adjMatrix =am_m1_BLCA,
                     prior = c(0.05,
                               0.05,
                               0.9),
                     nCPh = 0,
                     nGV = 10,
                     pmr = TRUE,
                     burnIn = 0.2,
                     iterations = 1000,
                     thinTo = 200,
                     progress = FALSE)

summary(mh_m1_BLCA)