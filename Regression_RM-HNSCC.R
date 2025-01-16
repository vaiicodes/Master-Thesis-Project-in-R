
options(stringsAsFactors = FALSE)

rm(list=ls())

gc()

abc <- read.table("C:/Users/HP/Downloads/GSE65021_non-normalized.txt.gz", fill=TRUE)
View(abc)

#make second row as colnames
colnames(abc) <- as.character(abc[1, ])


# Remove the first row
abc <- abc[-1, ]

#remove all Pval columns
abc <- abc[, !grepl("Pval", names(abc))]

colnames(abc) <- substr(colnames(abc), 1, 12)

colnames(abc) <- gsub("7964427070_A", "GSM1585940", colnames(abc))
colnames(abc) <- gsub("7964427070_B", "GSM1585941", colnames(abc))
colnames(abc) <- gsub("7964427070_C", "GSM1585942", colnames(abc))
colnames(abc) <- gsub("7964427070_D", "GSM1585943", colnames(abc))
colnames(abc) <- gsub("7964427070_E", "GSM1585944", colnames(abc))
colnames(abc) <- gsub("7964427070_F", "GSM1585945", colnames(abc))
colnames(abc) <- gsub("7964427070_I", "GSM1585946", colnames(abc))
colnames(abc) <- gsub("7964427070_K", "GSM1585947", colnames(abc))
colnames(abc) <- gsub("7964427070_L", "GSM1585948", colnames(abc))


colnames(abc) <- gsub("7964427073_A", "GSM1585949", colnames(abc))
colnames(abc) <- gsub("7964427073_B", "GSM1585950", colnames(abc))
colnames(abc) <- gsub("7964427073_C", "GSM1585951", colnames(abc))
colnames(abc) <- gsub("7964427073_D", "GSM1585952", colnames(abc))
colnames(abc) <- gsub("7964427073_E", "GSM1585953", colnames(abc))
colnames(abc) <- gsub("7964427073_F", "GSM1585954", colnames(abc))
colnames(abc) <- gsub("7964427073_G", "GSM1585955", colnames(abc))
colnames(abc) <- gsub("7964427073_H", "GSM1585956", colnames(abc))
colnames(abc) <- gsub("7964427073_I", "GSM1585957", colnames(abc))
colnames(abc) <- gsub("7964427073_K", "GSM1585958", colnames(abc))
colnames(abc) <- gsub("7964427073_L", "GSM1585959", colnames(abc))

colnames(abc) <- gsub("7964427082_A", "GSM1585960", colnames(abc))
colnames(abc) <- gsub("7964427082_B", "GSM1585961", colnames(abc))
colnames(abc) <- gsub("7964427082_E", "GSM1585962", colnames(abc))
colnames(abc) <- gsub("7964427082_F", "GSM1585963", colnames(abc))
colnames(abc) <- gsub("7964427082_G", "GSM1585964", colnames(abc))
colnames(abc) <- gsub("7964427082_H", "GSM1585965", colnames(abc))
colnames(abc) <- gsub("7964427082_I", "GSM1585966", colnames(abc))
colnames(abc) <- gsub("7964427082_J", "GSM1585967", colnames(abc))
colnames(abc) <- gsub("7964427082_K", "GSM1585968", colnames(abc))
colnames(abc) <- gsub("7964427082_L", "GSM1585969", colnames(abc))


colnames(abc) <- gsub("7964427089_A", "GSM1585970", colnames(abc))
colnames(abc) <- gsub("7964427089_B", "GSM1585971", colnames(abc))
colnames(abc) <- gsub("7964427089_E", "GSM1585972", colnames(abc))
colnames(abc) <- gsub("7964427089_F", "GSM1585973", colnames(abc))
colnames(abc) <- gsub("7964427089_G", "GSM1585974", colnames(abc))
colnames(abc) <- gsub("7964427089_H", "GSM1585975", colnames(abc))
colnames(abc) <- gsub("7964427089_I", "GSM1585976", colnames(abc))
colnames(abc) <- gsub("7964427089_J", "GSM1585977", colnames(abc))
colnames(abc) <- gsub("7964427089_K", "GSM1585978", colnames(abc))
colnames(abc) <- gsub("7964427089_L", "GSM1585979", colnames(abc))


probe_matrix <- read.table("D:/Master Thesis/cetuximab/probeMatrix.txt")
View(probe_matrix)

#make second row as colnames
colnames(probe_matrix) <- as.character(probe_matrix[1, ])


# Remove the first row
probe_matrix <- probe_matrix[-1, ]

library("GEOquery")
GPL14951 <- getGEO("GPL14951",destdir=".")
genename <- Table(GPL14951)[,c(1,11)]
expSet <- as.data.frame(probe_matrix)
expSet$ID <- rownames(expSet)
regress<-merge(x = genename,y=expSet,by="ID")
view(regress)

gene_matrix <- read.table("D:/Master Thesis/cetuximab/geneMatrix.txt")
View(gene_matrix)

colnames(gene_matrix) <- as.character(gene_matrix[1, ])


# Remove the first row
gene_matrix <- gene_matrix[-1, ]


specific_genes <- gene_matrix %>% 
  filter(geneNames %in% c("COL17A1", "FOSL1", "IL8", "INHBA", "ITGA3", "ITGA5", "ITGB4", "KLF6", "LAMA3", "LAMB3", "LAMC2", "MT2A", "ODC1", "PHLDA1", "PLEK2", "RGS2", "SERPINE1", "SPHK1"))

view(specific_genes)

# Assuming your data frame is called specific_genes
#colnames(specific_genes)[colnames(specific_genes) == "geneNames"] <- "ID"

specific_genes <- t(specific_genes)


colnames(specific_genes) <- as.character(specific_genes[1, ])
specific_genes <- specific_genes[-1, ]

view(specific_genes)

specific_genes <- as.data.frame(specific_genes)
specific_genes <- specific_genes %>% rownames_to_column(var = "ID")



SPHK1 <- specific_genes[, c("ID", "SPHK1")]
SPHK1 <- SPHK1[order(as.numeric(SPHK1[,"SPHK1"])), ]
view(SPHK1)

# Assuming your data frame is called ITGB4
SPHK1$SPHK1[1:20] <- "low"
# Assuming your data frame is called ITGB4
SPHK1$SPHK1[21:40] <- "high"
view(SPHK1)

gene_list <- list(COL17A1, FOSL1, IL8, INHBA, ITGA3, ITGA5, ITGB4, KLF6, LAMA3, LAMB3, LAMC2, MT2A, ODC1, PHLDA1, PLEK2, RGS2, SERPINE1, SPHK1)

# Merging the data frames by the "ID" column
genes <- Reduce(function(x, y) merge(x, y, by = "ID"), gene_list)
view(genes)

write.csv(genes, file = "genes.csv", row.names= TRUE)

main <- read.csv("genes.csv")
View(main)

install.packages("Seurat")
install.packages("gtools")
install.packages("patchwork")
install.packages("ggpubr")
install.packages("harmony")
install.packages("plyr")

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
library(openxlsx)
library(tibble)

setwd("D:/Master Thesis/Vaidehi_Thesis/cetuximab/Bossi_Cetuximab/Data")

data <- read.csv("genes.csv")
View(data)

data <- data[, -1]

data2 <- read.table("coxInput.txt")
View(data2)

# Make the second row as column names
colnames(data2) <- as.character(data2[1, ])

#remove the second row
data2 <- data2[-1, ]

#change the column name "id" to "ID"
colnames(data2)[colnames(data2) == "id"] <- "ID"

#only select the following columns 
data_final <- data2[, c("ID", "Age","Gender","Stage","Grade","Radiotherapy","PFS")]
View(data_final)

#merge data and data_final

merged <- merge(data, data_final, by = "ID")
merged <- na.omit(merged)
View(merged)

#for regression model
merged$PFS <- revalue(merged$PFS, c("Short"=1))
merged$PFS <- revalue(merged$PFS, c("Long"=0))

#for violin plots
merged$PFS <- revalue(merged$PFS, c("Short"=0))
merged$PFS <- revalue(merged$PFS, c("Long"=1))
merged$PFS <- as.numeric(merged$PFS)
View(merged)

write.table(merged, "cetuximab.txt", sep = "\t", row.names = FALSE)


##################################
#####REGRESSION MODEL#############
##################################

merged$PFS <- as.numeric(merged$PFS) #Convert PFS to numeric for model fitting

#Load additional libraries for survival and visualization
#install.packages("forestmodel")
library(forestmodel) 
library(ggplot2)
library(dplyr) 
library(survival) #For survival analysis and model fitting

merged <- read.table("cetuximab.txt", header = TRUE)

#Fit a logistic regression model with PFS as the dependent variable
linear_model <- glm(PFS~ITGB4+Age+Gender+Stage+Grade+Radiotherapy, 
                    family="binomial", data=merged)

#Generate a forest plot to visualize model coefficients and confidence intervals
forest_model(linear_model)


######for all genes using loop######

# Assuming your gene list is a character vector named 'gene_names'
gene_names <- c("COL17A1", "FOSL1", "IL8", "INHBA", "ITGA3", "ITGA5", "ITGB4", "KLF6", "LAMA3", "LAMB3", "LAMC2", "MT2A", "ODC1", "PHLDA1", "PLEK2", "RGS2", "SERPINE1", "SPHK1")

# Iterate through each gene in the gene list
for (gene_name in gene_names) {
  
  # Exclude the current gene from the formula
  gene_formula <- paste("PFS ~ Age + Gender + Stage + Grade + Radiotherapy +", gene_name)
  
  # Create formula dynamically
  formula <- as.formula(gene_formula)
  
  # Fit linear model for the current gene
  linear_model <- glm(formula, family = "binomial", data = merged)
  
  # Perform additional tasks if needed
  
  # Display forest model results
  forest_model_result <- forest_model(linear_model)
  
  # Plot the forest model result
  plot(forest_model_result)
}

??forest_model


##################################
###########VIOLIN PLOTS###########
##################################

# Load necessary libraries
library(ggplot2)
library(ggpubr)
library(dplyr)

# Example data preparation (replace 'merged' and column names accordingly)
merged <- read.table("cetuximab.txt", header = TRUE, sep = "\t")

# Gene names list
gene_names <- c("COL17A1", "FOSL1", "IL8", "INHBA", "ITGA3", "ITGA5", "ITGB4", "KLF6", "LAMA3", 
                "LAMB3", "LAMC2", "MT2A", "ODC1", "PHLDA1", "PLEK2", "RGS2", "SERPINE1", "SPHK1")

# Loop through each gene and create violin plot
for (gene in gene_names) {
  # Filter out groups with fewer than two data points
  filtered_merged <- merged %>%
    group_by(PFS) %>%
    filter(n() >= 2)
  
  # Create plot only if there are sufficient data points
  if (nrow(filtered_merged) > 0) {
    p <- ggplot(filtered_merged, aes(x = factor(PFS), y = !!sym(gene), fill = factor(PFS))) +
      geom_violin(alpha = 0.25, position = position_dodge(width = .75), size = 1, color = "black") + 
      geom_point(position = position_jitter(width = 0.2), alpha = 0.5) +
      geom_boxplot(notch = FALSE, outlier.size = -1, color = "black", lwd = 1.2, alpha = 0.7) + 
      labs(x = "PFS",
           y = paste(gene, "Expression")) +
      scale_fill_manual(values = c("blue", "#EC3338"), labels = c("Short PFS", "Long PFS")) +  # Custom fill colors
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

#############################################
######AUC CURVE IN LOOP FOR ALL GENES########
#############################################

# Load required library
library(pROC)

# Read data
rt <- read.table("cetuximab.txt", sep = "\t", header = TRUE, check.names = FALSE, row.names = 1)

# List of gene names
gene_names <- c("COL17A1", "FOSL1", "IL8", "INHBA", "ITGA3", "ITGA5", "ITGB4", "KLF6", "LAMA3", 
                "LAMB3", "LAMC2", "MT2A", "ODC1", "PHLDA1", "PLEK2", "RGS2", "SERPINE1", "SPHK1")

# Loop through each gene name
for (gene in gene_names) {
  # Calculate ROC curve
  roc <- roc(rt$PFS, rt[, gene], smooth = FALSE)
  
  # Plot ROC curve
  plot(roc, print.auc = TRUE, print.auc.x = 0.5, print.auc.y = 0.5, 
       auc.polygon = TRUE, auc.polygon.col = "#99CCFF", 
       grid = FALSE, legacy.axes = TRUE, main = gene)
}


