library(data.table,lib.loc="/mnt/ceph/fern5249/Rpackages")
library(na.tools,lib.loc="/mnt/ceph/fern5249/Rpackages")
library(MRGN,lib.loc="/mnt/ceph/fern5249/Rpackages")
library(MRTrios,lib.loc="/mnt/ceph/fern5249/Rpackages")
library(baycn,lib.loc="/mnt/ceph/fern5249/Rpackages")

meth<- as.data.frame(fread("/mnt/ceph/fern5249/GDCdata/TCGA-BLCA/Analysis/split.names.BLCA.meth.logit.txt"))
dim(meth)
# Apply partial matching to column names
colnames(meth)[5:ncol(meth)] <- sapply(strsplit(colnames(meth)[5:ncol(meth)], "-"), function(parts) paste(parts[1:3], collapse="-"))
meth[1, 1:10]


#Gene Expression dataset
gene<- fread("/mnt/ceph/fern5249/GDCdata/TCGA-BLCA/blca_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt")
dim(gene)
# Apply partial matching to column names
colnames(gene)[3:ncol(gene)] <- sapply(strsplit(colnames(gene)[3:ncol(gene)], "-"), function(parts) paste(parts[1:3], collapse="-"))
gene[1, 1:10]


#CNA dataset
cna<- fread("/mnt/ceph/fern5249/GDCdata/TCGA-BLCA/blca_tcga_pan_can_atlas_2018/data_cna.txt")
dim(cna)
# Apply partial matching to column names
colnames(cna)[3:ncol(cna)] <- sapply(strsplit(colnames(cna)[3:ncol(cna)], "-"), function(parts) paste(parts[1:3], collapse="-"))
cna[1, 1:10]



#clinical data
clinical<-fread("/mnt/ceph/fern5249/GDCdata/TCGA-BLCA/Analysis/split_new_data_clinical_patient.txt")

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

#Read in the PC score matrix
pc.meth<- read.table("/mnt/ceph/fern5249/GDCdata/TCGA-BLCA/Analysis/PCA.meth.txt", header = TRUE)

pc.gene<- read.table("/mnt/ceph/fern5249/GDCdata/TCGA-BLCA/Analysis/PCA.gene.exp.txt",header=TRUE)
# Remove duplicated rows from pc.meth
result_meth <- removeDupsPCs(pc.meth, NULL)
pc.meth <- result_meth$pc
duplicated_row_names_meth <- result_meth$duplicated_row_names

# Remove duplicated rows from pc.gene
result_gene <- removeDupsPCs(pc.gene, NULL)
pc.gene <- result_gene$pc
duplicated_row_names_gene <- result_gene$duplicated_row_names

# Remove corresponding columns from BLCA.cna
BLCA.cna <- BLCA.cna[, !colnames(BLCA.cna) %in% duplicated_row_names_meth, with = FALSE]

# Remove corresponding columns from BLCA.gene
BLCA.gene <- BLCA.gene[, !colnames(BLCA.gene) %in% duplicated_row_names_gene, with = FALSE]

#reading in the Trios data
trios <- data.frame(fread("/mnt/ceph/fern5249/GDCdata/TCGA-BLCA/Analysis/trio.final.protein.coding.txt"))

#read in the indices table
meth.table<- fread("/mnt/ceph/fern5249/GDCdata/TCGA-BLCA/Analysis/meth.table.txt", drop = 1)
gene.table<- fread("/mnt/ceph/fern5249/GDCdata/TCGA-BLCA/Analysis/gene.exp.table.txt", drop = 1)


#read in the sig pcs data
meth.sig.asso.pcs<- readRDS("/mnt/ceph/fern5249/GDCdata/TCGA-BLCA/Analysis/meth.sig.asso.pcs.RData")
gene.sig.asso.pcs<- readRDS("/mnt/ceph/fern5249/GDCdata/TCGA-BLCA/Analysis/gene.exp.sig.asso.pcs.RData")


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
              # Remove the "sex" column if it doesn't have at least 2 levels
      if (nlevels(as.factor(final.mat.complete$sex)) < 2){
        final.mat.complete <- final.mat.complete[, -which(names(final.mat.complete) == "sex")]
      }
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

#data=datamatrix(meth, gene, cna, trios[1,], pc.meth, pc.gene, meth.sig.asso.pcs[[1]], gene.sig.asso.pcs[[1]],clinical, meth.table, gene.table,age.col=5, race.col=26,sex.col=6)


baycn_summary_results <- function(data, trios, p, writeToFile = FALSE, file) {
  results <- data.frame(Trio = integer(), Model = character(), stringsAsFactors = FALSE)  # Initialize an empty dataframe to store results
    # write to file if writeToFile is TRUE
    if (writeToFile) {
      results.colnames <- c("Trio", "Inferred.Model")
      write.table(t(results.colnames), file = file, sep = "\t", row.names = FALSE, col.names = FALSE, append = FALSE, quote = FALSE)
    }
  # Begin the loop for rows in trios 
  for (i in 1:nrow(trios)) {
    cat("Trio", i, ":\n")
    
    data <- datamatrix(meth, gene, cna, trios[i,], pc.meth, pc.gene, meth.sig.asso.pcs[[1]], gene.sig.asso.pcs[[1]], clinical, meth.table, gene.table, age.col = 5, race.col = 26, sex.col = 6)
      
    
    if (is.null(data)) {
      next  # Skip to the next trio if data is NULL
    }
    
    # Remove race column
    t1 <- data[, -which(names(data) %in% c("race"))]
    t1$sex <- ifelse(t1$sex == "Male", 0, 1)
    
    # Unlist the data matrix 
    t2 <- unlist(t1)
    t3 <- matrix(t2, byrow = FALSE, nrow = nrow(t1))
    colnames(t3) <- colnames(t1)
    
    # Move confounders up so they can be treated as genetic variants
    t4 <- t3[, c(1, 4:ncol(t3), 2, 3)]
    
    # Create an adjacency matrix
    am_m1 <- matrix(0, nrow = ncol(t4), ncol = ncol(t4))
    am_m1[, (ncol(t4) - 1):ncol(t4)] <- 1
    am_m1[(ncol(t4)), (ncol(t4))] <- 0
    am_m1[(ncol(t4) - 1), (ncol(t4) - 1)] <- 0
    
    # Run the mhEdge function
    mh_m1 <- mhEdge(data = t4,
                    adjMatrix = am_m1,
                    prior = c(0.05, 0.05, 0.9),
                    nCPh = 0,
                    nGV = ncol(t4) - 2,
                    pmr = TRUE,
                    burnIn = 0.2,
                    iterations = 1000,
                    thinTo = 200,
                    progress = FALSE)
    posterior_probs <- mh_m1@posteriorES
    
    # Extract the relevant rows for posterior probabilities
    edges <- posterior_probs[c(1, 2, nrow(posterior_probs)), ]
    
    if (!is.null(edges)) { # Check if edges are not NULL
      if (edges[1, "zero"] == max(edges[1, -1])) {#Edge between V1-->V2
        if (edges[2, "zero"] == max(edges[2, -1])) {#Edge between V1-->V3
          if (sum(edges[3, c("zero", "one")]) > 0.5) {# If there's a edge between V2-V3
            model_type <- "M4"
          } else {
            model_type <- "M3" #no Edge V2-V3
          }
        } else if (edges[2, "two"] == max(edges[2, -1])) {#no Edge between V1-V3
          if (sum(edges[3, c("zero", "one")]) > 0.5) {# If there's a edge between V2-V3
            if (edges[3, "zero"]-edges[3,"one"]>p) {#Edge V2-->V3
              model_type <- "M1.1"
            } else if (edges[3, "one"]-edges[3, "zero"]>p) {#Edge V3-->V2
              model_type <- "M2.1"
            }else if (abs(edges[3, "zero"] - edges[3, "one"]) < p) {#Edge V2<-->V3
              model_type<-"Other"
            }
          } else {
            model_type <- "M0.1" #no Edge V2-V3
          }
        }
      } else if (edges[1, "two"] == max(edges[1, -1])) {# no edge between V1-V2
        if (edges[2, "zero"] == max(edges[2, -1])) {
          if (sum(edges[3, c("zero", "one")]) > 0.5) {
            if (edges[3, "zero"]-edges[3, "one"]>p) {
              model_type <- "M2.2"
            } else if (edges[3, "one"]-edges[3,"zero"]>p) {
              model_type <- "M1.2"
            }else if (abs(edges[3, "zero"] - edges[3, "one"]) < p) {
              model_type<-"Other"
            }
          } else {
            model_type <- "M0.2"
          }
        } else if (edges[2, "two"] == max(edges[2, -1])) {#no Edge between V1-V3
          model_type <- "Other"
        }
      }
    }
    # Store results only if edges are not NULL
    results <- rbind(results, data.frame(Trio= i, Inferred.Model = model_type))
    
      # Write to file if writeToFile is TRUE
      if (writeToFile) {
        write.table(data.frame(Trio= i, Inferred.Model = model_type), file = file, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)
      }
   }

  return(results)
}


baycn.results <- baycn_summary_results(data, trios[30001:60000, ], p=0.2, writeToFile =TRUE, file= "/mnt/ceph/fern5249/GDCdata/TCGA-BLCA/Analysis1/baycn2.txt")


###########################################################################################3

data1<-fread("/mnt/ceph/fern5249/GDCdata/TCGA-BLCA/analyze.trios.BLCA.split1.txt")

data2<-fread("/mnt/ceph/fern5249/GDCdata/TCGA-BLCA/baycn11.txt")


# Function to create a confusion matrix
create_confusion_matrix <- function(data1, data2) {
  # Extract trios and inferred models from both datasets
  trios_data1 <- data1$Index
  models_data1 <- data1$Inferred.Model
  
  trios_data2 <- data2$Index
  models_data2 <- data2$Inferred.Model
  
  # List of all possible models
  all_models <- c("M0.1", "M0.2", "M1.1", "M1.2", "M2.1", "M2.2", "M3", "M4", "Other")
  
  # Create an empty confusion matrix
  confusion_matrix <- matrix(0, nrow = length(all_models), ncol = length(all_models), 
                             dimnames = list(all_models, all_models))
  
  # Iterate through each trio
  for (i in 1:length(trios_data1)) {
    trio <- trios_data1[i]
    model_data1 <- models_data1[i]
    model_data2 <- models_data2[i]
    
    # Update the confusion matrix based on the models
    confusion_matrix[model_data1, model_data2] <- confusion_matrix[model_data1, model_data2] + 1
  }
  
  return(confusion_matrix)
}

# Usage example:
# Assuming data1 and data2 are the datasets containing trios and inferred models
confusion_matrix <- create_confusion_matrix(data1, data2)
print(confusion_matrix)
