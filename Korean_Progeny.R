#Load the required package if not already installed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#Install the "progeny" package from Bioconductor
BiocManager::install("progeny")

#Disable converting strings to factors
options(stringsAsFactors = FALSE)

#Clear the R workspace, removing all objects
rm(list=ls())


#Load the Seurat object from a saved RData file
load("scRNAEPI.RData")

#Load the seurat library
library(Seurat)
#Create a t-SNE plot with labels
DimPlot(scRNA, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()



#Compute the Progeny activity scores and add them to our Seurat object
#as a new assay called Progeny. 
library(progeny)
scRNA <- progeny(scRNA, scale=FALSE, organism="Human", top=500, perm=1, 
                 return_assay = TRUE)

## We can now directly apply Seurat functions in our Progeny scores. 
## For instance, we scale the pathway activity scores. 
scRNA <- Seurat::ScaleData(scRNA, assay = "progeny") 

## We transform Progeny scores into a data frame to better handling the results
library(tidyr)
library(dplyr)
library(tidyverse) 

progeny_scores_df <- 
  as.data.frame(t(GetAssayData(scRNA, slot = "scale.data", 
                               assay = "progeny"))) 


progeny_scores_df

#Load the gene expression data
load("data_raw.RData")

#Get the matrix of expression of ITGA5, INHBA, and SERPINE1, and transfer it to dataframe
#We create a data frame with the specification of the genes that interest 
#each gene to match with the Progeny scores.
Genes = dataan[c("CXCL8", "FOSL1", "KLF6", "COL17A1", "INHBA", "ITGA3", "ITGA5", "ITGB4", "LAMA3", 
                "LAMB3", "LAMC2", "MT2A", "PLEK2", "PHLDA1", "SERPINE1", "SPHK1","ODC1","RGS2"),]
Genes = t(Genes)
#Genes
Genes = data.frame(Genes)
rownames_to_column(Genes,"Cell") 
View(Genes)

#Find common samples between Progeny scores and gene expression data
sameSample=intersect(row.names(progeny_scores_df), row.names(Genes))	
progeny_scores_df=progeny_scores_df[sameSample,,drop=F]		
Genes=Genes[sameSample,1:(ncol(Genes)),drop=F]		


#Perform correlation analysis
outTab=data.frame()		#create empty data frame to store the results
for(i in colnames(progeny_scores_df)){	#iterate over progeny pathways
  for(j in colnames(Genes)){	            #iterate over selected genes
    x=as.numeric(progeny_scores_df[,i]) #extract progeny scores from current pathway(i)
    y=as.numeric(Genes[,j])             #and gene(j)
    corT=cor.test(x,y,method="spearman") #Perform Spearmann correlation test
    cor=corT$estimate      #extract correlation coefficient, p-value 
    pvalue=corT$p.value    #and significant text based on p-value
    text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
    outTab=rbind(outTab,cbind(Gene=j, Pathway=i, cor, text, pvalue))
  }	
}		

#Create a correlation heatmap
outTab$cor=as.numeric(outTab$cor)		#ensure that the 'cor' column in 'outTab' is numeric
pdf(file="COR_18GENES.pdf", width=8, height=12)	#create a pdf file to save the heatmap
ggplot(outTab, aes(Gene, Pathway)) + 		#create a ggplot object using 'outTab' for plotting
  geom_tile(aes(fill = cor), colour = "grey", size = 1)+	
  scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") + 	
  geom_text(aes(label=text),col ="black",size = 3) +	
  theme_minimal() +    #apply minimal theme, removing background elements	
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),	
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),   
        axis.text.y = element_text(size = 10, face = "bold")) +   
  labs(fill =paste0("***  p<0.001","\n", "**  p<0.01","\n", " *  p<0.05","\n", "\n","Correlation")) + 
  scale_x_discrete(position = "bottom") 	
dev.off()		#close the pdf device and save the heatmap as cor.pdf


#####
