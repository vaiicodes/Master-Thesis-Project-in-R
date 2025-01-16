
#### AUC/ROC and Violin Plots for Cetuximab Response in PDX Cohort #####

#Set options to not automatically convert strings to factors
options(stringsAsFactors = FALSE)

#Remove all objects from the current R environment
rm(list=ls())
gc() #Free up memory

#Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(harmony)
library(readxl)
library(plyr)

#Set the working directory to the specified path
setwd("D:/Master Thesis/cetuximab/xeno_cetuximab")

# Install and load the 'affy' package, used for Affymetrix microarray data analysis
install.packages("affy")
library(affy)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("affy")

#Define the folder path where CEL files are stored
folder_path <- "D:/Master Thesis/cetuximab/GSE84713_RAW"

#List all CEL files in the current working directory
cel_files <- list.files(pattern="\\.CEL$")

#Read each CEL file into a list, with header and tab-separated values
cel_data <- lapply(cel_files, read.table, header=TRUE, sep="\t")

View(cel_data) #View the imported CEL data

#Read Affymetrix CEL files from the specified path
cel_data <- ReadAffy(celfile.path="D:/Master Thesis/cetuximab/GSE84713_RAW")

#Get and display sample names from the CEL data
sampleNames(cel_data)

#Get the number of samples in the CEL data
N=length(cel_data)

#Perform Robust Multi-array Average (RMA) normalization on the CEL data
eset.rma<-rma(cel_data)

#Extract expression values from the normalized data
cancer_exprs<- exprs(eset.rma)

#Get probe IDs from the row names of the expression matrix
probeid<-rownames(cancer_exprs)

#Add probe IDs as the first column in the expression matrix
cancer_exprs<-cbind(probeid, cancer_exprs)

#Save the expression data in a text file
write.table(cancer_exprs, file="cancer_exprs.txt", sep="\t", quote=F,row.names=F)
View(cancer_exprs)

#Set row names as probe IDs and remove the first column
row.names(cancer_exprs) = cancer_exprs[,1]
cancer_exprs = cancer_exprs[,-1]

#Shorten column names to first 10 characters
colnames(cancer_exprs) <- substr(colnames(cancer_exprs), 1,10)

#Install and load GEOquery package to access GEO database data
BiocManager::install("GEOquery")
library("GEOquery")

#Download the GPL570 platform file and extract gene names
GPL570 <- getGEO("GPL570",destdir=".")
genename <- Table(GPL570)[,c(1,11)]

#Convert expression data to a data frame and add probe IDs
exprSet <- as.data.frame(cancer_exprs)
exprSet$ID <- rownames(exprSet)

#Merge expression data with gene names based ofn probe IDs
express<-merge(x = genename,y=exprSet,by="ID")
express <- express[,-1] #remove the ID column after merging

#Set gene symbols as row names in the expression data
row.names(express) = express[,1]

#Save the expression data with gene symbols in a text file
write.table(express,file="GeneSymbol.txt",sep="\t")

#Read the gene symbol file as input data
input <- read.table("D:/Master Thesis/cetuximab/xeno_cetuximab/GeneSymbol.txt")
view(input)

#Define a list of genes of interest for the analysis
genes_of_interest <- c("COL17A1", "FOSL1", "CXCL8", "INHBA", "ITGA3", "ITGA5", "ITGB4", "KLF6", "LAMA3", "LAMB3", "LAMC2", "MT2A", "ODC1", "PHLDA1", "PLEK2", "RGS2", "SERPINE1", "SPHK1")

#Subset data to select only the desired genes and their expression values
selected_genes <- subset(input, Gene.Symbol %in% genes_of_interest)
View(selected_genes)

#Save the selected genes to a csv file
write.csv(selected_genes, file="gene_input.csv")

#Calculate the mean expression for each gene symbol
final_genes <- selected_genes %>%
  group_by(Gene.Symbol) %>%
  summarise_all(mean)

View(final_genes)

#Write the final gene expression data to a csv file
write.csv(final_genes, file="final_input.csv")

#Read the final input csv file and transpose it for further analysis
data <- read.csv("final_input.csv")
data <- data[, -1]
data <- t(data)

#Set the second row as column names
colnames(data) <- as.character(data[1, ])

#Remove the first row
data <- data[-1, ]

#Convert to a data frame and add a row ID column
data <- as.data.frame(data)
data <- data %>% rownames_to_column(var = "ID")

View(data)

#Extract and view the expression data for a specific gene, COL17A
COL17A1 <- data[, c("ID", "COL17A1")]
view(COL17A1)

#Read a manually prepared Cetuximab Responder data
data2 <- read_xlsx("regression_input.xlsx")
View(data2)

#Merge the gene expression data with cetuximab responder data
merged <- merge(data, data2, by = "ID")
View(merged)


# Convert all columns in the merged dataset to numeric
merged <- mutate_all(merged, as.numeric)

#Save the merged data in a text file
write.table(merged, "merged.txt", sep = "\t", row.names = FALSE)


#############################################
######AUC CURVE IN LOOP FOR ALL GENES########

#Load required library
library(pROC)

# Read the merged data
rt <- read.table("merged.txt", sep = "\t", header = TRUE, check.names = FALSE, row.names = 1)

# List of gene names
gene_names <- c("COL17A1", "FOSL1", "CXCL8", "INHBA", "ITGA3", "ITGA5", "ITGB4", "KLF6", "LAMA3", 
                "LAMB3", "LAMC2", "MT2A", "ODC1", "PHLDA1", "PLEK2", "RGS2", "SERPINE1", "SPHK1")

# Loop through each gene to plot the ROC Curve
for (gene in gene_names) {
  # Calculate ROC curve
  roc <- roc(rt$Cetuximab_response, rt[, gene], smooth = FALSE)
  
  # Plot ROC curve
  plot(roc, print.auc = TRUE, print.auc.x = 0.5, print.auc.y = 0.5, 
       auc.polygon = TRUE, auc.polygon.col = "#99CCFF", 
       grid = FALSE, legacy.axes = TRUE, main = gene)
}


############################################
########### VIOLIN PLOTS LOOP ##############

# Load necessary libraries
library(ggplot2)
library(ggpubr)
library(dplyr)

#Read merged data for plotting
merged <- read.table("merged.txt", header = TRUE, sep = "\t")

# Gene names list
gene_names <- c("COL17A1", "FOSL1", "CXCL8", "INHBA", "ITGA3", "ITGA5", "ITGB4", "KLF6", "LAMA3", 
                "LAMB3", "LAMC2", "MT2A", "ODC1", "PHLDA1", "PLEK2", "RGS2", "SERPINE1", "SPHK1")

# Loop through each gene and create violin plot
for (gene in gene_names) {
  # Filter out groups with fewer than two data points
  filtered_merged <- merged %>%
    group_by(Cetuximab_response) %>%
    filter(n() >= 2)
  
  # Create plot only if there are sufficient data points
  if (nrow(filtered_merged) > 0) {
    p <- ggplot(filtered_merged, aes(x = factor(Cetuximab_response), y = !!sym(gene), fill = factor(Cetuximab_response))) +
      geom_violin(alpha = 0.25, position = position_dodge(width = .75), size = 1, color = "black") + 
      geom_point(position = position_jitter(width = 0.2), alpha = 0.5) +
      geom_boxplot(notch = FALSE, outlier.size = -1, color = "black", lwd = 1.2, alpha = 0.7) + 
      labs(x = "Cetuximab Response",
           y = paste(gene, "Expression")) +
      scale_fill_manual(values = c("blue", "orange"), labels = c("non-responder", "responder")) +  # Custom fill colors
      theme_pubr()
    
    # Add p-value
    p <- p + stat_compare_means(comparisons = list(c("0", "1")), label = "p.format", method = "t.test")
    
    # Print plot
    print(p)
  } else {
    cat("Insufficient data points for", gene, "\n")
  }
}

# Close any open graphics devices
dev.off()

