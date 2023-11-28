##.....Trio_gene_type_sep.....##
library(data.table,lib.loc="/mnt/ceph/oluw5072/Rpackages/MRGN_R/")

#read the datasets
trios <- data.frame(fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/Trios.final.txt"))

gene.info <- data.frame(fread("/mnt/ceph/oluw5072/GDCdata/mart_export_merged_lncRNA_fixed.txt"))

#find the unique genes in the trios
gene.names <- unique(trios[,1])

#create a text file for new trios
write.table(t(colnames(trios)), file = "/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/trio.final.protein.coding.txt", sep = "\t", row.names = FALSE,
            col.names = FALSE,quote=FALSE)
#initialize count
count = 0

#loop to go over each gene
for(i in 1:length(gene.names)){
  
  #find if the gene is in the bio mart data
  row <- which(gene.info$Gene.name == gene.names[i])
  row
  
  #check the length
  if(length(row) > 0){
    
    #check if the gene type is either protein coding or lncRNA
    if(gene.info$Gene.type[row[1]] == "protein_coding" | gene.info$Gene.type[row[1]] == "lncRNA"){
      
      
      print(i)
      
      #increase the count for each gene of the two types
      count = count + 1
      
      #find the rows in trios which has this gene
      row.in.trios <- which(trios[,1] == gene.names[i])
      
      #save the values in those rows to a new variables
      final <- trios[row.in.trios,]
      
      #write final to the new trios file
      write.table(final, file = "/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/trio.final.protein.coding.txt", sep = "\t", row.names = FALSE,
                  col.names = FALSE, append = TRUE, quote=FALSE)
    }
  }
  
}
