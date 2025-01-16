

#Load required libraries
library(Seurat)
library(dplyr)
library(SeuratData)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(R.utils)
library(harmony)

#Decompress the data file
gunzip("GSE188737_hnscc_seu.RData.gz", remove = FALSE)

#Load the scRNA-seq data object
load("GSE188737_hnscc_seu.RData")
View(hnscc_seu)



#Simplify metadata by selecting only relevant columns
a = hnscc_seu@meta.data[,1:3]
a

b = hnscc_seu@meta.data[,-(4:11)]
b

hnscc_seu@meta.data = hnscc_seu@meta.data[,-(4:11)]

scRNA = hnscc_seu

View(scRNA)


#Load additional metadata file for the cells
gunzip("GSE188737_metacells.csv.gz")

celldata = read.csv("GSE188737_metacells.csv")
View(celldata)

#Set row names of the metadata to match cell IDs
rownames(celldata) <- celldata$X    #找人改错过
celldata <-celldata[,-1]
head(celldata)

#Add the additional metadata to the Seurat object
scRNA=AddMetaData(scRNA, metadata = celldata)
view(scRNA)

#Save the seurat object for further use
save(scRNA,file = 'scRNA0.RData')


#Load the seurat object in the environment
load("scRNA0.RData")
dim(scRNA)


#Quality Control
#Calculate mitochondrial and ribosomal content percentages
scRNA_Epi[["percent.mt"]] <- PercentageFeatureSet(scRNA_Epi, pattern = "^MT-")
scRNA_Epi[["percent.rb"]] <- PercentageFeatureSet(scRNA_Epi, pattern = "^RP[LS]")

#Define haemoglobin genes and calculate their percentage in cells
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(scRNA_Epi))
scRNA_Epi[["percent.HB"]]<-PercentageFeatureSet(scRNA_Epi, features=HB.genes) 


#Visualize quality control metrics 
col.num <- length(levels(scRNA_Epi@active.ident))   



violin <- VlnPlot(scRNA_Epi,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb","percent.HB"), 
                  cols =rainbow(col.num),  
                  
                  pt.size = 0.01, 
                  ncol = 4) +    
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 


violin
#Save the violin plot
ggsave("vlnplotQCSingapore.png", plot = violin, width = 12, height = 6)  

#Use scatter plots to inspect relationships between QC metrics
plot1=FeatureScatter(scRNA_Epi, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2=FeatureScatter(scRNA_Epi, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2
plot3=FeatureScatter(scRNA_Epi, feature1 = "nCount_RNA", feature2 = "percent.rb")
plot3
plot4=FeatureScatter(scRNA_Epi, feature1 = "nCount_RNA", feature2 = "percent.HB")
plot4
pearplot <- CombinePlots(plots = list(plot1, plot2, plot3,plot4), nrow=2, legend="none") 
pearplot


#Filter cells based on quality control thresholds
scRNA_Epi <- subset(scRNA_Epi,                                  
                subset = nFeature_RNA > 200 & 
                  nFeature_RNA < 8000 & 
                  percent.mt < 10 & 
                  percent.HB < 3 &
                  nCount_RNA < 100000)
scRNA_Epi



###########################################################
##############Normalization and Scaling####################
###########################################################

#Normalize data for downstream analysis

#Identify highly variable genes (HVGs)
scRNA_Epi <- NormalizeData(scRNA_Epi, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA_Epi <- FindVariableFeatures(scRNA_Epi, 
                              selection.method = "vst", 
                              nfeatures = 2000,
                              mean.cutoff=c(0.0125,3),
                              dispersion.cutoff =c(1.5,Inf))

#Visualize the top 20 variable genes
top20 <- head(VariableFeatures(scRNA), 20)
plot1 <- VariableFeaturePlot(scRNA)
plot1
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE, size=3.0)
plot2
feat_20 <- CombinePlots(plots = list(plot1, plot2),legend="bottom")
feat_20
ggsave(filename = 'SingaporeTop20.pdf',plot = feat_20,he=10,wi=15)



#Scale the Data to make it comparable across features
scale.genes <-  rownames(scRNA_Epi)
scRNA_Epi <- ScaleData(scRNA_Epi, features = scale.genes)



###########################################################
##############PCA & Dimentionality Reduction###############
###########################################################


#Run PCA for dimensional reduction and plot the elbow plot for optimal dimensions
scRNA_Epi <- RunPCA(scRNA_Epi, features = VariableFeatures(scRNA_Epi)) 
dimplot1 <- DimPlot(scRNA_Epi, reduction = "pca", group.by = "P_Mid") 

elbowplot1 <- ElbowPlot(scRNA_Epi, ndims=50, reduction="pca") 
sc_pca <- dimplot1+elbowplot1
sc_pca
ggsave(filename = 'SingaporeSC_pca.pdf',plot = sc_pca,he=10,wi=15)

#Optional visualization of PCA loading
VizDimLoadings(scRNA_Epi, dims = 1:2, nfeatures = 20, reduction = "pca")

#Visualize PCA dimensions using a heatmap
par(mar = c(3, 3, 1, 1))
DimHeatmap(scRNA, dims = 1:20, cells = 500, balanced = TRUE)


#Run Harmony for batch effect correction
system.time({scRNA_Epi <- RunHarmony(scRNA_Epi, group.by.vars = "P_Mid")})
#system.time({scRNA <- RunHarmony(scRNA, group.by.vars = "tissue.type")})


#Run t-SNE for final visualization, correcting for batch effects
scRNA_Epi <- scRNA_Epi %>% 
  RunTSNE(reduction = "harmony", dims = 1:50) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:50) %>% 
  FindClusters(resolution = 0.2) %>% 
  identity()

#Plot t-SNE results to visualize clusters

p1 <- DimPlot(scRNA_Epi, reduction = "tsne", group.by = "P_Mid", label = TRUE)
p1


#Save the final Seurat object
save(scRNA,file = 'scRNA1.RData')

load("scRNA1.RData")
dim(scRNA)
view(scRNA)

#Subset the data to include only Epithelial cells
scRNA_Epi <- subset(scRNA, subset = cell_type == "Epithelial")

view(scRNA_Epi)
dim(scRNA_Epi)
summary(scRNA_Epi)

#Save epithelial cell subset for future use
save(scRNA_Epi,file = 'scRNAEpi.RData')


#Generate FeaturePlots to show gene expression in specific tissues

#Primary Cancer Tissue (CA)
scRNA_CA = subset(scRNA_Epi, tissue.type %in% "CA")
dim(scRNA_CA)
Expression_CA <- FeaturePlot(scRNA_CA, features = c("CXCL8", "FOSL1", "KLF6", "COL17A1", "INHBA", "ITGA3", "ITGA5", "ITGB4", "LAMA3", 
                                                    "LAMB3", "LAMC2", "MT2A", "PLEK2", "PHLDA1", "SERPINE1", "SPHK1","ODC1","RGS2"), min.cutoff = "q9", reduction="tsne")
Expression_CA

#LYMPH NODE TISSUE (LN)
scRNA_LN = subset(scRNA_Epi, tissue.type %in% "LN")
dim(scRNA_LN)

Expression_LN <- FeaturePlot(scRNA_LN, features = c("CXCL8", "FOSL1", "KLF6", "COL17A1", "INHBA", "ITGA3", "ITGA5", "ITGB4", "LAMA3", 
                                                  "LAMB3", "LAMC2", "MT2A", "PLEK2", "PHLDA1", "SERPINE1", "SPHK1","ODC1","RGS2"), min.cutoff = "q9", reduction="tsne")
Expression_LN


genes_of_interest <- c("CXCL8", "FOSL1", "KLF6", "COL17A1", "INHBA", "ITGA3", "ITGA5", "ITGB4", "LAMA3", "LAMB3", "LAMC2", "MT2A", "PLEK2", "PHLDA1", "SERPINE1", "SPHK1", "ODC1", "RGS2")

#Generate DotPlot to show expression of selected genes by tissue type
scRNA_Epi$P_Mid <- factor(scRNA_Epi$P_Mid, levels = c("P", "M"))
SingaporeDotPlot2 <- DotPlot(scRNA_Epi, 
                          features = genes_of_interest, 
                          cols = c("blue", "red"), 
                          dot.scale = 8,col.min = -1,
                          col.max = 2) + coord_flip()
SingaporeDotPlot2




