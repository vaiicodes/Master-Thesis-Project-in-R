
#Disable automatic conversion of strings to factors
options(stringsAsFactors = FALSE)

#Clear the R environment and free up memory
rm(list=ls())
gc()

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



#Set the working directory
#setwd("D:/Master Thesis/Korean Dataset")

#Load expression data from the Korean dataset
data <- read.table("GSE181919_UMI_counts.txt", 
                   header = T, 
                   row.names = 1)

#Convert to data frame, then to a matrix for sparse storage
datan = data.frame(data)
dataan <- as(as.matrix(datan),"dgCMatrix")

#Save the file as R Data for future use
save(dataan,file = 'data_raw.RData')

#Load the R data into the environment
load("data_raw.RData")

#Create a Seurat object from the raw data matrix
scRNA = CreateSeuratObject(dataan, 
                           project="HNSC",   
                           min.cells = 3,  #Only include genes expressed in at least 3 cells
                           min.features = 300) #Only include cells with at least 300 detected features (genes)

class(scRNA)   
view(scRNA)   

#Save the Seurat object in a R data format
save(scRNA,file = 'scRNA_raw.RData')

#Load metadata for Korean Dataset cells
celldata=read.table(file = "GSE181919_Barcode_metadata.txt",
                    sep="\t",
                    header = T)

#Format the cell identifiers for consistency
rownames(celldata) <- gsub("-",".",celldata$X)    
celldata <-celldata[,-1]
head(celldata)
View(celldata)

#Filter rows where cell.type is "epithelial"
epithelial_cells <- celldata[celldata$cell.type == "Epithelial.cells", ]

#Count the total number of sample.id values (including duplicates)
total_sample_ids <- nrow(epithelial_cells)

#Display the count
total_sample_ids

#Add metadata to the Seurat object
scRNA=AddMetaData(scRNA, metadata = celldata)
head(scRNA)
dim(scRNA)
View(scRNA)

#Save the seurat object in R data format
save(scRNA,file = 'scRNA0.RData')


#Load the seurat object in the environment
load("scRNA0.RData")
dim(scRNA)

#Calculate mitochondrial and ribosomal gene expression percentages
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
scRNA[["percent.rb"]] <- PercentageFeatureSet(scRNA, pattern = "^RP[LS]")

#Calculate the percentage of hemoglobin genes, which may indicate low-quality or stressed cells
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(scRNA))
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes) 

#Count the number of unique cluster identities in the Seurat object
col.num <- length(levels(scRNA@active.ident)) 

#Visualize quality control metrics using violin plots
violin <- VlnPlot(scRNA,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb","percent.HB"), 
                  cols =rainbow(col.num),  
                  
                  pt.size = 0.01, 
                  ncol = 4) +    
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 

violin
#ggsave("vlnplot_before_qc.pdf", plot = violin, width = 12, height = 6) 
ggsave("vlnplot_after_qc.png", plot = violin, width = 12, height = 6)  


#Use scatter plots to inspect relationships between QC metrics
plot1=FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2=FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2
plot3=FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.rb")
plot3
plot4=FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.HB")
plot4
pearplot <- CombinePlots(plots = list(plot1, plot2, plot3,plot4), nrow=2, legend="none") 
pearplot


#Subset data to remove low-quality cells based on specific QC thresholds
scRNA <- subset(scRNA,                                
                subset = nFeature_RNA > 200 & 
                  nFeature_RNA < 8000 & 
                  percent.mt < 10 & 
                  percent.HB < 3 &
                  nCount_RNA < 100000)
scRNA

#Normalize data for downstream analysis
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)

#Identify highly variable genes (HVGs)
scRNA <- FindVariableFeatures(scRNA, 
                              selection.method = "vst", 
                              nfeatures = 2000,
                              mean.cutoff=c(0.0125,3),
                              dispersion.cutoff =c(1.5,Inf))

top20 <- head(VariableFeatures(scRNA), 20)

#Plot top 20 variable features for inspection
plot1 <- VariableFeaturePlot(scRNA)
plot1
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE, size=3.0)
plot2
feat_20 <- CombinePlots(plots = list(plot1, plot2),legend="bottom")
feat_20
ggsave(filename = 'KoreanTop20.pdf',plot = feat_20,he=10,wi=15)


#Scale the data to make it comparable across features
scale.genes <-  rownames(scRNA)
scRNA <- ScaleData(scRNA, features = scale.genes)



###########################################################
###########################################################


#Run PCA for dimensional reduction and plot the elbow plot for optimal dimensions
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA)) 
dimplot1 <- DimPlot(scRNA, reduction = "pca", group.by = "sample.id") 
??DimPlot

elbowplot1 <- ElbowPlot(scRNA, ndims=50, reduction="pca") 
sc_pca <- dimplot1+elbowplot1
sc_pca
ggsave(filename = 'KoreanSC_pca.pdf',plot = sc_pca,he=10,wi=15)

#Optional visualization of PCA loading
VizDimLoadings(scRNA, dims = 1:2, nfeatures = 20, reduction = "pca")

#Visualize PCA dimensions using a heatmap
par(mar = c(3, 3, 1, 1))
DimHeatmap(scRNA, dims = 1:20, cells = 500, balanced = TRUE)


#Run Harmony for batch effect correction
system.time({scRNA <- RunHarmony(scRNA, group.by.vars = "sample.id")})
#system.time({scRNA <- RunHarmony(scRNA, group.by.vars = "tissue.type")})


#Run t-SNE and UMAP for final visualization, correcting for batch effects
scRNA <- scRNA %>% 
  RunTSNE(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.2) %>% 
  identity()



#Plot t-SNE results to visualize clusters
p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "seurat_clusters", label = TRUE)
p1

p2 <- DimPlot(scRNA, reduction = "tsne", group.by = "cell.type", label = TRUE)
p2


#Save the final Seurat object
save(scRNA,file = 'scRNA1.RData')

#Load the saved Seurat object for further analysis
load('scRNA1.RData')
dim(scRNA)

#Subset the data to include only HPV-negative cells
scRNA=subset(scRNA,hpv %in% 'HPV-')   
dim(scRNA)
# Further subset for "Epithelial cell" and "Malignant" based on the cell.type column
scRNA_Epi_Mal <- subset(scRNA, cell.type %in% c("Epithelial.cells", "Malignant.cells"))

dim(scRNA_Epi_Mal)

#Display tissue type breakdown in the epithelial cell subset
table(scRNA_Epi_Mal@meta.data$tissue.type)
View(scRNA)

#Save epithelial cell subset for future use
save(scRNA_Epi_Mal, file = "scRNAEPI.Rdata")

#Load the saved epithelial cell subset for analysis
load("scRNAEPI.Rdata")

dim(scRNA_Epi_Mal)

#Normalize the Korean epithelial cell dataset
scRNA <- NormalizeData(scRNA_Epi_Mal, normalization.method = "LogNormalize", scale.factor = 10000)

#Identify variable features in the Korean epithelial cell dataset
scRNA <- FindVariableFeatures(scRNA_Epi_Mal, 
                                     selection.method = "vst", 
                                     nfeatures = 2000,
                                     mean.cutoff=c(0.0125,3),
                                     dispersion.cutoff =c(1.5,Inf))
scale.genes <-  rownames(scRNA_Korean)

#Scale data in the Korean epithelial cell dataset
scRNA <- ScaleData(scRNA_Epi_Mal, features = scale.genes)

#Display the epithelial cell subset
View(scRNA_Epi)

#Generate t-SNE plot for epithelial cells by tissue type
sc_tsne_picked = DimPlot(scRNA_Epi_Mal,
                         cols = c("#0000FF","#DC143C","pink","orange"),
                         group.by='tissue.type',        
                         reduction="tsne",
                         
                         pt.size = 0.2,
                         #    label = "T", 
                         label.size = 5) +
  ggtitle("tissue_type", )  +     
  theme(plot.title = element_text(hjust = 0.5,size=15))+  
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        axis.ticks = element_blank()
  )

sc_tsne_picked

#Generate FeaturePlots to show gene expression in specific tissues

#Primary Cancer Tissue (CA)
scRNA_CA = subset(scRNA_Epi_Mal, tissue.type %in% 'CA')
Expression_CA <- FeaturePlot(scRNA_CA, features = c("CXCL8", "FOSL1", "KLF6", "COL17A1", "INHBA", "ITGA3", "ITGA5", "ITGB4", "LAMA3", 
                                                    "LAMB3", "LAMC2", "MT2A", "PLEK2", "PHLDA1", "SERPINE1", "SPHK1")
                             , min.cutoff = "q9", reduction="tsne")

Expression_CA

#NORMAL TISSUE (NL)
scRNA_NL = subset(scRNA_Epi_Mal, tissue.type %in% 'NL')

Expression_NL <- FeaturePlot(scRNA_NL, features = c("CXCL8", "FOSL1", "KLF6", "COL17A1", "INHBA", "ITGA3", "ITGA5", "ITGB4", "LAMA3", 
                                                    "LAMB3", "LAMC2", "MT2A", "PLEK2", "PHLDA1", "SERPINE1", "SPHK1"), min.cutoff = "q9", reduction="tsne")
Expression_NL

#LEUKOPLAKIA TISSUE (LP)
scRNA_LP = subset(scRNA_Epi_Mal, tissue.type %in% 'LP')

Expression_LP <- FeaturePlot(scRNA_LP, features = c("CXCL8", "FOSL1", "KLF6", "COL17A1", "INHBA", "ITGA3", "ITGA5", "ITGB4", "LAMA3", 
                                                    "LAMB3", "LAMC2", "MT2A", "PLEK2", "PHLDA1", "SERPINE1", "SPHK1"), min.cutoff = "q9", reduction="tsne")
Expression_LP

#LYMPH NODE TISSUE (LN)
scRNA_LN = subset(scRNA_Epi_Mal, tissue.type %in% 'LN')

Expression_LN <- FeaturePlot(scRNA_LN, features = c("CXCL8", "FOSL1", "KLF6", "COL17A1", "INHBA", "ITGA3", "ITGA5", "ITGB4", "LAMA3", 
                                                    "LAMB3", "LAMC2", "MT2A", "PLEK2", "PHLDA1", "SERPINE1", "SPHK1"), min.cutoff = "q9", reduction="tsne")
Expression_LN

#Generate DotPlot to show expression of selected genes by tissue type
genes_of_interest <- c("CXCL8", "FOSL1", "KLF6", "COL17A1", "INHBA", "ITGA3", "ITGA5", "ITGB4", "LAMA3", "LAMB3", "LAMC2", "MT2A", "PLEK2", "PHLDA1", "SERPINE1", "SPHK1")

#Subset data to include specific tissue types (CA, LP, NL) for a focused DotPlot
scRNA_Korean = subset(scRNA_Epi_Mal, tissue.type %in% c("CA","LP","NL"))
scRNA_Korean

#Visualize a dot plot
scRNA_Korean$tissue.type <- factor(scRNA_Korean$tissue.type, levels = c("NL", "LP", "CA"))
KoreanDotPlot <- DotPlot(scRNA_Korean, 
                          features = genes_of_interest, 
                          cols = c("blue", "red"), 
                          dot.scale = 8, 
                          group.by = "tissue.type") + coord_flip()
KoreanDotPlot
