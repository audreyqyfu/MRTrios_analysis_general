library(data.table,lib.loc="/mnt/ceph/oluw5072/Rpackages/MRGN_R/")
library(na.tools,lib.loc="/mnt/ceph/oluw5072/Rpackages/MRGN_R/")

##read in the analyzeTrios 
LUAD1.analyzeTrios <- fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/analyze.trios.LUAD_1.txt", header=TRUE)

LUAD2.analyzeTrios <- fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/analyze.trios.LUAD_2.txt", header=TRUE)

LUAD3.analyzeTrios <- fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/analyze.trios.LUAD_3.txt", header=TRUE)

## combine the data files
LUAD.analyzeTrios <- rbind(LUAD1.analyzeTrios, LUAD2.analyzeTrios, LUAD3.analyzeTrios)


#Save the combined data
write.table(LUAD.analyzeTrios, file = "/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/LUAD.analyze_trios.txt", sep = "\t", row.names = FALSE,
            col.names = TRUE)

#Load the combined_data
LUAD<-data.frame(fread("/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/LUAD.analyze_trios.txt"))

#Check the list of unique Inferred models
model_classes<-sapply(LUAD$Inferred.Model,list)
unique_classes<-unique(model_classes)
unique_classes

#Check the number of trios for each Inferred model
model_counts<-table(LUAD$Inferred.Model)
model_counts

#Check the distribution of Inferred models using percentages
model_percentages<-(model_counts/sum(model_counts))*100
model_percentages

#Save the Inferred model percentages
write.table(model_percentages, file = "/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/distribution_of_inferred_models.txt", sep = "\t", row.names = FALSE,
            col.names = TRUE)

## read in the distribution data in desktop
Model_LUAD <-read_excel("Documents/Research/distribution_of_Inferred_Model_LUAD.xlsx")

##plot barplot
ggplot(Model_LUAD, aes(x=Var1, y=Freq))+ geom_bar(stat = "identity", fill = "lightblue") + labs(title = "Distribution of Inferred Models_LUAD", x="Inferred Model" , y ="Percentages") + geom_text(aes(label=round(Freq,1)), vjust =-0.5, size =3)


##################################################
#Check the list of Total PC.Count
Total.PC.Count_classes<-sapply(LUAD$Total.PC.Count,list)
unique_Total.PC.Count_classes<-unique(Total.PC.Count_classes)
unique_Total.PC.Count_classes


#Check the number of trios for Total.PC.Count
Total.PC.Count_counts<-table(LUAD$Total.PC.Count)
Total.PC.Count_counts

#Check the distribution of counts using percentages
Total.PC.Count_percentages<-(Total.PC.Count_counts/sum(Total.PC.Count_counts))*100
Total.PC.Count_percentages

#Save the percentages
write.table(Total.PC.Count_percentages, file = "/mnt/ceph/oluw5072/GDCdata/TCGA-LUAD/Analysis/distribution_of_Total.PC.Count.txt", sep = "\t", row.names = FALSE,
            col.names = TRUE)
## read in the distribution data in desktop
Total.PC.Count_LUAD <-read_excel("Documents/Research/distribution_of_Total.PC.Count_LUAD.xlsx")

## plot barplot
ggplot(Total.PC.Count_percentages, aes(x=Var1, y=Freq))+ geom_bar(stat = "identity", fill = "lightblue") + labs(title = "Distribution of Inferred Models_LUAD", x="Inferred Model" , y ="Percentages") + geom_text(aes(label=round(Freq,1)), vjust =-0.5, size =3)

