# MRTrios_analysis_general

The scripts should be used in the following order:

•	DownloadHarmonizedData.R : this file downloads the harmonized methylation data from TCGA, combine it with additional data (gotten from BRCA), and save the processed data for further analysis.

•	DataProcessing.R: this file performs logit transformation of the methylation data.

•	mainTrioMatch.R: this file generates the trio data matrix by integrating the CNA, methylation, and gene expression data, with each trio in a separate line and each line containing the row numbers of the probe or gene in the input data.

•	trio.gene.type.sep.R: this file takes the trios data matrix as input and keeps only "protein coding genes" or "lncRNAs".

•	main.findPCs.R: this file calculates the principal component (PC) score matrix for methylation and separately for gene expression data, and identifies PCs that are significantly associated each trio.

•	Formatting.PatientID.R: this file splits the patient ID to match each other.

•	main.analyzeTrios.R: this file performs causal network inference for trios and their associated confounders (i.e., PCs) to infer the causal models.

•	trioInferenceResults.Plots.R :this file summarizes the trio inference results and generating the bar plots. 
