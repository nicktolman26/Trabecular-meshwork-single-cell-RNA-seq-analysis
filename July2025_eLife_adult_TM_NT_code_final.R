library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library("biomaRt")
library(plyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(data.table)
library(scCustomize)
library(viridis)
library(data.table)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggpubr)
library(reticulate)
use_condaenv("r-reticulate")
library(clusterProfiler)
library(org.Mm.eg.db)
library(utils)
library(VennDiagram)
library(circlize)
library(nichenetr)

############ Load in data and create seurat objects - B6 background data only

p60_b6_john_wt1 <- Read10X(data.dir = "~/Desktop/Lmx1b/Development/Single_Cell_Sequencing/SN073_WT/SN073_cellranger_count_outs/filtered_feature_bc_matrix")
p60_b6_john_wt1 <- CreateSeuratObject(counts = p60_b6_john_wt1, project = "p60_b6_john_wt1", min.cells = 3, min.features = 200)

p60_b6_john_wt2 <- Read10X(data.dir = "~/Desktop/Lmx1b/Development/Single_Cell_Sequencing/SN081_WT/SN081_cellranger_count_outs/filtered_feature_bc_matrix")
p60_b6_john_wt2 <- CreateSeuratObject(counts = p60_b6_john_wt2, project = "p60_b6_john_wt2", min.cells = 3, min.features = 200)

p60_129_john_wt1 <- Read10X(data.dir = "~/Desktop/Lmx1b/Development/Single_Cell_Sequencing/SN083_WT/SN083_cellranger_count_outs/filtered_feature_bc_matrix")
p60_129_john_wt1 <- CreateSeuratObject(counts = p60_129_john_wt1, project = "p60_129_john_wt1", min.cells = 3, min.features = 200)

p90_b6_stamer_wt1 <- Read10X(data.dir = "~/Desktop/Lmx1b/Development/Single_Cell_Sequencing/Stamer_lab_3mo_B6_scseq")
p90_b6_stamer_wt1 <- CreateSeuratObject(counts = p90_b6_stamer_wt1, project = "p90_b6_stamer_wt1", min.cells = 3, min.features = 200)

############# Remove cells with bad mitochondrial DNA % or likely doublet and visualize by UMAP all samples

p60_b6_john_wt1[["percent.mt"]] <- PercentageFeatureSet(p60_b6_john_wt1, pattern = "^mt-")
VlnPlot(p60_b6_john_wt1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p60_b6_john_wt1 <- subset( p60_b6_john_wt1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
p60_b6_john_wt1 <- NormalizeData( p60_b6_john_wt1, normalization.method = "LogNormalize", scale.factor = 10000)
p60_b6_john_wt1 <- FindVariableFeatures( p60_b6_john_wt1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames( p60_b6_john_wt1)
p60_b6_john_wt1 <- ScaleData( p60_b6_john_wt1, features = all.genes)
p60_b6_john_wt1 <- RunPCA( p60_b6_john_wt1, features = VariableFeatures(object = p60_b6_john_wt1))
p60_b6_john_wt1<- FindNeighbors(p60_b6_john_wt1, dims = 1:10)
p60_b6_john_wt1 <- FindClusters(p60_b6_john_wt1, resolution = 0.1)
p60_b6_john_wt1<- RunUMAP( p60_b6_john_wt1, dims = 1:10)
DimPlot(p60_b6_john_wt1, reduction = "umap")
table(Idents( p60_b6_john_wt1))
p60_b6_john_wt1 <- AddMetaData(object = p60_b6_john_wt1, metadata = "B6", col.name = "Strain")
p60_b6_john_wt1 <- AddMetaData(object = p60_b6_john_wt1, metadata = "John", col.name = "Lab")
p60_b6_john_wt1 <- AddMetaData(object = p60_b6_john_wt1, metadata = "WT", col.name = "Genotype")

p60_b6_john_wt2[["percent.mt"]] <- PercentageFeatureSet(p60_b6_john_wt2, pattern = "^mt-")
VlnPlot(p60_b6_john_wt2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p60_b6_john_wt2 <- subset( p60_b6_john_wt2, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
p60_b6_john_wt2 <- NormalizeData( p60_b6_john_wt2, normalization.method = "LogNormalize", scale.factor = 10000)
p60_b6_john_wt2 <- FindVariableFeatures( p60_b6_john_wt2, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames( p60_b6_john_wt2)
p60_b6_john_wt2 <- ScaleData( p60_b6_john_wt2, features = all.genes)
p60_b6_john_wt2 <- RunPCA( p60_b6_john_wt2, features = VariableFeatures(object = p60_b6_john_wt2))
p60_b6_john_wt2<- FindNeighbors(p60_b6_john_wt2, dims = 1:10)
p60_b6_john_wt2 <- FindClusters(p60_b6_john_wt2, resolution = 0.1)
p60_b6_john_wt2<- RunUMAP( p60_b6_john_wt2, dims = 1:10)
DimPlot(p60_b6_john_wt2, reduction = "umap")
table(Idents( p60_b6_john_wt2))
p60_b6_john_wt2 <- AddMetaData(object = p60_b6_john_wt2, metadata = "B6", col.name = "Strain")
p60_b6_john_wt2 <- AddMetaData(object = p60_b6_john_wt2, metadata = "John", col.name = "Lab")
p60_b6_john_wt2 <- AddMetaData(object = p60_b6_john_wt2, metadata = "WT", col.name = "Genotype")

p60_129_john_wt1[["percent.mt"]] <- PercentageFeatureSet(p60_129_john_wt1, pattern = "^mt-")
VlnPlot(p60_129_john_wt1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p60_129_john_wt1 <- subset( p60_129_john_wt1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
p60_129_john_wt1 <- NormalizeData( p60_129_john_wt1, normalization.method = "LogNormalize", scale.factor = 10000)
p60_129_john_wt1 <- FindVariableFeatures( p60_129_john_wt1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames( p60_129_john_wt1)
p60_129_john_wt1 <- ScaleData( p60_129_john_wt1, features = all.genes)
p60_129_john_wt1 <- RunPCA( p60_129_john_wt1, features = VariableFeatures(object = p60_129_john_wt1))
p60_129_john_wt1<- FindNeighbors(p60_129_john_wt1, dims = 1:10)
p60_129_john_wt1 <- FindClusters(p60_129_john_wt1, resolution = 0.1)
p60_129_john_wt1<- RunUMAP( p60_129_john_wt1, dims = 1:10)
DimPlot(p60_129_john_wt1, reduction = "umap")
table(Idents( p60_129_john_wt1))
p60_129_john_wt1 <- AddMetaData(object = p60_129_john_wt1, metadata = "129", col.name = "Strain")
p60_129_john_wt1 <- AddMetaData(object = p60_129_john_wt1, metadata = "John", col.name = "Lab")

p90_b6_stamer_wt1[["percent.mt"]] <- PercentageFeatureSet(p90_b6_stamer_wt1, pattern = "^mt-")
VlnPlot(p90_b6_stamer_wt1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p90_b6_stamer_wt1 <- subset( p90_b6_stamer_wt1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
p90_b6_stamer_wt1 <- NormalizeData( p90_b6_stamer_wt1, normalization.method = "LogNormalize", scale.factor = 10000)
p90_b6_stamer_wt1 <- FindVariableFeatures( p90_b6_stamer_wt1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames( p90_b6_stamer_wt1)
p90_b6_stamer_wt1 <- ScaleData( p90_b6_stamer_wt1, features = all.genes)
p90_b6_stamer_wt1 <- RunPCA( p90_b6_stamer_wt1, features = VariableFeatures(object = p90_b6_stamer_wt1))
p90_b6_stamer_wt1<- FindNeighbors(p90_b6_stamer_wt1, dims = 1:10)
p90_b6_stamer_wt1 <- FindClusters(p90_b6_stamer_wt1, resolution = 0.1)
p90_b6_stamer_wt1<- RunUMAP( p90_b6_stamer_wt1, dims = 1:10)
DimPlot(p90_b6_stamer_wt1, reduction = "umap")
table(Idents( p90_b6_stamer_wt1))
p90_b6_stamer_wt1 <- AddMetaData(object = p90_b6_stamer_wt1, metadata = "Stamer", col.name = "Lab")

##### Integrate John Lab single cell seq P60 with 2 B6 and 1 129 runs
anchors_adult_wt_john <- FindIntegrationAnchors(object.list = list(p60_b6_john_wt1, p60_b6_john_wt2, p60_129_john_wt1), dims = 1:10)
combined_adult_wt_john <- IntegrateData(anchorset = anchors_adult_wt_john, dims = 1:10)
DefaultAssay(combined_adult_wt_john) <- "integrated"
combined_adult_wt_john <- ScaleData(combined_adult_wt_john, verbose = FALSE)
combined_adult_wt_john <- RunPCA(combined_adult_wt_john, npcs = 30, verbose = FALSE)
combined_adult_wt_john <- RunUMAP(combined_adult_wt_john, reduction = "pca", dims = 1:10, metric = 'correlation')
combined_adult_wt_john <- FindNeighbors(combined_adult_wt_john, reduction = "pca", dims = 1:10)
combined_adult_wt_john <- FindClusters(combined_adult_wt_john, resolution = 0.1)
combined_adult_wt_john$seurat_clusters <- Idents(combined_adult_wt_john)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
p1 <- DimPlot(combined_adult_wt_john, reduction = "umap", group.by = "orig.ident", cols = brewer.pal(n=4,name="Set2"))
p2 <- DimPlot(combined_adult_wt_john, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
nb.cols <- 10
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
Idents(combined_adult_wt_john) <- "seurat_clusters"
DimPlot(combined_adult_wt_john, reduction = "umap", cols = mycolors)
new.cluster.ids <- c("TM_containing","Epithelial","Epithelial","TM_containing","Iris_and_CB","Endothelial","Immune","Epithelial","Neuron")
names(new.cluster.ids) <- levels(combined_adult_wt_john)
combined_adult_wt_john <- RenameIdents(combined_adult_wt_john, new.cluster.ids)
combined_adult_wt_john$CellType <- Idents(combined_adult_wt_john)
combined_adult_wt_john$CellType <- factor(combined_adult_wt_john$CellType, levels = c("TM_containing","Epithelial","Iris_and_CB","Endothelial","Immune","Neuron"))
mycolors <- c("#8dd3c7","#80b1d3","#fdb462","#b3de69","#bebada","#fb8072")
Idents(combined_adult_wt_john) <- "CellType"
DimPlot(combined_adult_wt_john, reduction = "umap", cols = mycolors, label = TRUE)

######## Can run this way and achieve same cluster resolution as above
######## Not used in paper though
######## Do not run for subsequent analysis
DefaultAssay(combined_adult_wt_john) <- "integrated"
combined_adult_wt_john <- ScaleData(combined_adult_wt_john, verbose = FALSE)
combined_adult_wt_john <- RunPCA(combined_adult_wt_john, npcs = 30, verbose = FALSE)
combined_adult_wt_john <- RunUMAP(combined_adult_wt_john, reduction = "pca", dims = 1:10, metric = 'correlation')
combined_adult_wt_john <- FindNeighbors(combined_adult_wt_john, reduction = "pca", dims = 1:10)
combined_adult_wt_john <- FindClusters(combined_adult_wt_john, resolution = 0.03)
combined_adult_wt_john$seurat_clusters <- Idents(combined_adult_wt_john)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
p1 <- DimPlot(combined_adult_wt_john, reduction = "umap", group.by = "orig.ident", cols = brewer.pal(n=4,name="Set2"))
p2 <- DimPlot(combined_adult_wt_john, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
nb.cols <- 10
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
Idents(combined_adult_wt_john) <- "seurat_clusters"
DimPlot(combined_adult_wt_john, reduction = "umap", cols = mycolors)
new.cluster.ids <- c("Epithelial","TM_containing","Endothelial","Iris_and_CB","Immune","Neuron")
names(new.cluster.ids) <- levels(combined_adult_wt_john)
combined_adult_wt_john <- RenameIdents(combined_adult_wt_john, new.cluster.ids)
combined_adult_wt_john$CellType <- Idents(combined_adult_wt_john)
combined_adult_wt_john$CellType <- factor(combined_adult_wt_john$CellType, levels = c("TM_containing","Epithelial","Iris_and_CB","Endothelial","Immune","Neuron"))
mycolors <- c("#8dd3c7","#80b1d3","#fdb462","#b3de69","#bebada","#fb8072")
Idents(combined_adult_wt_john) <- "CellType"
DimPlot(combined_adult_wt_john, reduction = "umap", cols = mycolors, label = TRUE)
#########
######### Pick up here
DefaultAssay(combined_adult_wt_john) <- "RNA"
Idents(combined_adult_wt_john) <- "CellType"
gene_list_plot <- c("Rho","Cd52","Egfl7","Tyrp1","Krt14","Pitx2","Tfap2b","Myoc","Acta2")
DotPlot(object = combined_adult_wt_john, features = gene_list_plot, cols = c("yellow","darkorchid3"))

DefaultAssay(combined_adult_wt_john) <- "RNA"
combined_adult_wt_john <- ScaleData(combined_adult_wt_john, verbose = FALSE)
Idents(combined_adult_wt_john) <- "CellType"
markers_combined_adult_wt_john <- FindAllMarkers(combined_adult_wt_john, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <- markers_combined_adult_wt_john %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(combined_adult_wt_john, features = top5$gene, group.colors = mycolors) +scale_fill_viridis()

####### Save all cells combined object
saveRDS(combined_adult_wt_john, file = "~/Desktop/Lmx1b/Development/Single_Cell_Sequencing/RDS_files/11172023_combined_adult_wt_john.RDS")

####### load all cells combined object
combined_adult_wt_john <- readRDS("~/Desktop/Lmx1b/Development/Single_Cell_Sequencing/RDS_files/11172023_combined_adult_wt_john.RDS")

############ Making POM-derived object
Idents(combined_adult_wt_john) <- "seurat_clusters"
combined_adult_wt_john_pom <- subset(combined_adult_wt_john, idents = c("0","3"))
DefaultAssay(combined_adult_wt_john_pom) <- "RNA"
combined_adult_wt_john_pom <- NormalizeData(combined_adult_wt_john_pom, normalization.method = "LogNormalize", scale.factor = 10000)
combined_adult_wt_john_pom <- FindVariableFeatures(combined_adult_wt_john_pom, nfeatures = 2000)
combined_adult_wt_john_pom <- ScaleData(combined_adult_wt_john_pom)
combined_adult_wt_john_pom <- RunPCA(combined_adult_wt_john_pom, features = VariableFeatures(object = combined_adult_wt_john_pom))
combined_adult_wt_john_pom <- FindNeighbors(combined_adult_wt_john_pom, dims = 1:10)
combined_adult_wt_john_pom <- FindClusters(combined_adult_wt_john_pom, resolution = 1)
set.seed(500)
combined_adult_wt_john_pom <- RunUMAP(combined_adult_wt_john_pom, dims = 1:10)
nb.cols <- 18
mycolors <- colorRampPalette(brewer.pal(9, "Paired"))(nb.cols)
DimPlot(combined_adult_wt_john_pom, reduction = "umap", cols = mycolors, raster = FALSE, label = TRUE)

#Because of the low cell numbers of ciliary muscle cells, they cluster with pericytes at 0.3 resolution. 
#Therefore, we clustered at 1 resolution and combined clusters to match 0.3 resolution except for the ciliary muscle
#Here is the clustering at 0.3 resolution
#Do not run for subsequent analysis
DefaultAssay(combined_adult_wt_john_pom) <- "RNA"
combined_adult_wt_john_pom <- NormalizeData(combined_adult_wt_john_pom, normalization.method = "LogNormalize", scale.factor = 10000)
combined_adult_wt_john_pom <- FindVariableFeatures(combined_adult_wt_john_pom, nfeatures = 2000)
combined_adult_wt_john_pom <- ScaleData(combined_adult_wt_john_pom)
combined_adult_wt_john_pom <- RunPCA(combined_adult_wt_john_pom, features = VariableFeatures(object = combined_adult_wt_john_pom))
combined_adult_wt_john_pom <- FindNeighbors(combined_adult_wt_john_pom, dims = 1:10)
combined_adult_wt_john_pom <- FindClusters(combined_adult_wt_john_pom, resolution = 0.3)
set.seed(500)
combined_adult_wt_john_pom <- RunUMAP(combined_adult_wt_john_pom, dims = 1:10)
nb.cols <- 17
mycolors <- colorRampPalette(brewer.pal(9, "Paired"))(nb.cols)
DimPlot(combined_adult_wt_john_pom, reduction = "umap", cols = mycolors, raster = FALSE, label = TRUE)
###########

########### Pick back up here
Idents(combined_adult_wt_john_pom) <- "seurat_clusters"
new.cluster.ids <- c("TM3","TM1","TM2","IS","IS","TM1","TM1","TM2","IS","TM3","TM2","TM3","Kera","SF","Peri","Schw","CE","CM")
names(new.cluster.ids) <- levels(combined_adult_wt_john_pom)
combined_adult_wt_john_pom <- RenameIdents(combined_adult_wt_john_pom, new.cluster.ids)
combined_adult_wt_john_pom$CellType <- Idents(combined_adult_wt_john_pom)
combined_adult_wt_john_pom$CellType <- factor(combined_adult_wt_john_pom$CellType, levels = c('TM1','TM2','TM3','IS','SF','Kera','CE','CM','Peri','Schw'))

Idents(combined_adult_wt_john_pom) <- "CellType"
mycolors <- c("#a6cee3","#e31a1c","#33a02c","#fb9a99","khaki3","#b2df8a","#cab2d6","#ff7f00","#1f78b4","#6a3d9a")
DimPlot(combined_adult_wt_john_pom, reduction = "umap", cols = mycolors)
DimPlot(combined_adult_wt_john_pom, reduction = "umap", cols = mycolors, split.by = "Strain")

################# Save file if needed
saveRDS(combined_adult_wt_john_pom, file = "~/Desktop/Lmx1b/Development/Single_Cell_Sequencing/RDS_files/11172023_combined_adult_wt_john_pom.RDS")
#################

################# Load file if needed
combined_adult_wt_john_pom <- readRDS("~/Desktop/Lmx1b/Development/Single_Cell_Sequencing/RDS_files/11172023_combined_adult_wt_john_pom.RDS")
#################

# Dendrogram

Idents(combined_adult_wt_john_pom) <- "CellType"
combined_adult_wt_john_pom_reduced <- subset(combined_adult_wt_john_pom, idents = c('TM1','TM2','TM3','IS','SF','Kera','CE','CM'))

Idents(combined_adult_wt_john_pom_reduced) <- "CellType"
markers_combined_adult_wt_john_pom_reduced <- FindAllMarkers(combined_adult_wt_john_pom_reduced, min.pct = 0.1, logfc.threshold = 0.25)
top500 <- markers_combined_adult_wt_john_pom_reduced %>% group_by(cluster) %>% top_n(n = 500, wt = avg_log2FC)

combined_adult_wt_john_pom_tree = BuildClusterTree(
  combined_adult_wt_john_pom_reduced,
  features = top500$gene,
  slot = "data",
)

PlotClusterTree(combined_adult_wt_john_pom_tree, direction = "downwards")

############ Heatmap of markergenes in each POM cluster

####### Select marker genes
Idents(combined_adult_wt_john_pom) <- "CellType"
gene_list_plot <- c("Chil1","Tnmd","Cldn10","Myoc","Fmod","Pgf","Inmt","Cygb","Crym","Gpc3","Lypd1","Tfap2b","Gpr137b","Tagln","Acta2")
DotPlot(object = combined_adult_wt_john_pom, features = gene_list_plot, idents = c('TM1','TM2','TM3'),cols = c("yellow","darkorchid3"))

# TM cells only
Idents(combined_adult_wt_john_pom) <- "CellType"
combined_adult_wt_john_tm <- subset(combined_adult_wt_john_pom, idents = c('TM1','TM2','TM3'))

Idents(combined_adult_wt_john_tm) <- "CellType"
markers_combined_adult_wt_john_tm <- FindAllMarkers(combined_adult_wt_john_tm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <- markers_combined_adult_wt_john_tm %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top10 <- markers_combined_adult_wt_john_tm %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top20 <- markers_combined_adult_wt_john_tm %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
mycolors <- c("#a6cee3","#e31a1c","#33a02c")
my_levels <- c("TM1","TM2","TM3")
levels(combined_adult_wt_john_tm) <- my_levels
DoHeatmap(combined_adult_wt_john_tm, features = top20$gene, group.colors = mycolors) +scale_fill_viridis()

# All TM containing cluster
Idents(combined_adult_wt_john_pom) <- "CellType"
markers_combined_adult_wt_john_pom <- FindAllMarkers(combined_adult_wt_john_pom, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top20 <- markers_combined_adult_wt_john_pom %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top10 <- markers_combined_adult_wt_john_pom %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top5 <- markers_combined_adult_wt_john_pom %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top100 <- markers_combined_adult_wt_john_pom %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
mycolors <- c("#a6cee3","#e31a1c","#33a02c","#fb9a99","khaki3","#b2df8a","#cab2d6","#ff7f00","#1f78b4","#6a3d9a")
my_levels <- c('TM1','TM2','TM3','IS','SF','Kera','CE','CM','Peri','Schw')
levels(combined_adult_wt_john_pom) <- my_levels
DoHeatmap(combined_adult_wt_john_pom, features = top5$gene, group.colors = mycolors) +scale_fill_viridis()

# Marker comparison across strains
Idents(combined_adult_wt_john_pom) <- "Strain"
combined_adult_wt_john_pom_b6 <- subset(combined_adult_wt_john_pom, idents = c("B6"))

Idents(combined_adult_wt_john_pom) <- "Strain"
combined_adult_wt_john_pom_129 <- subset(combined_adult_wt_john_pom, idents = c("129"))

Idents(combined_adult_wt_john_pom_b6) <- "CellType"
markers_combined_adult_wt_john_pom_b6 <- FindAllMarkers(combined_adult_wt_john_pom_b6, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top20 <- markers_combined_adult_wt_john_pom_b6 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top5 <- markers_combined_adult_wt_john_pom_b6 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top10 <- markers_combined_adult_wt_john_pom_b6 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
my_levels <- c('TM1','TM2','TM3','IS','SF','Kera','CE','CM','Peri','Schw')
levels(combined_adult_wt_john_pom_b6) <- my_levels
DoHeatmap(combined_adult_wt_john_pom_b6, features = top10$gene, group.colors = mycolors) +scale_fill_viridis()

Idents(combined_adult_wt_john_pom_129) <- "CellType"
markers_combined_adult_wt_john_pom_129 <- FindAllMarkers(combined_adult_wt_john_pom_129, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top20 <- markers_combined_adult_wt_john_pom_129 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top5 <- markers_combined_adult_wt_john_pom_129 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top10 <- markers_combined_adult_wt_john_pom_129 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
my_levels <- c('TM1','TM2','TM3','IS','SF','Kera','CE','CM','Peri','Schw')
levels(combined_adult_wt_john_pom_129) <- my_levels
DoHeatmap(combined_adult_wt_john_pom_129, features = top10$gene, group.colors = mycolors) +scale_fill_viridis()

#### TM cell DE analysis

Idents(combined_adult_wt_john_pom) <- "CellType"
combined_adult_wt_john_tm <- subset(combined_adult_wt_john_pom, idents = c("TM1","TM2","TM3"))

Idents(combined_adult_wt_john_tm) <- "CellType"
new.cluster.ids <- c("TM1","TM","TM")
names(new.cluster.ids) <- levels(combined_adult_wt_john_tm)
combined_adult_wt_john_tm <- RenameIdents(combined_adult_wt_john_tm, new.cluster.ids)
combined_adult_wt_john_tm$CellTypeTM1 <- Idents(combined_adult_wt_john_tm)

Idents(combined_adult_wt_john_tm) <- "CellType"
new.cluster.ids <- c("TM","TM2","TM")
names(new.cluster.ids) <- levels(combined_adult_wt_john_tm)
combined_adult_wt_john_tm <- RenameIdents(combined_adult_wt_john_tm, new.cluster.ids)
combined_adult_wt_john_tm$CellTypeTM2 <- Idents(combined_adult_wt_john_tm)

Idents(combined_adult_wt_john_tm) <- "CellType"
new.cluster.ids <- c("TM","TM","TM3")
names(new.cluster.ids) <- levels(combined_adult_wt_john_tm)
combined_adult_wt_john_tm <- RenameIdents(combined_adult_wt_john_tm, new.cluster.ids)
combined_adult_wt_john_tm$CellTypeTM3 <- Idents(combined_adult_wt_john_tm)

Idents(combined_adult_wt_john_tm) <- "CellTypeTM1"
combined_adult_wt_john_tm1_vs_tm <- FindMarkers(combined_adult_wt_john_tm, ident.1 = "TM1", ident.2 = "TM")
combined_adult_wt_john_tm1_vs_tm
write.csv(combined_adult_wt_john_tm1_vs_tm, 'combined_adult_wt_john_tm1_vs_tm.csv')

Idents(combined_adult_wt_john_tm) <- "CellTypeTM2"
combined_adult_wt_john_tm2_vs_tm <- FindMarkers(combined_adult_wt_john_tm, ident.1 = "TM2", ident.2 = "TM")
combined_adult_wt_john_tm2_vs_tm
write.csv(combined_adult_wt_john_tm2_vs_tm, 'combined_adult_wt_john_tm2_vs_tm.csv')

Idents(combined_adult_wt_john_tm) <- "CellTypeTM3"
combined_adult_wt_john_tm3_vs_tm <- FindMarkers(combined_adult_wt_john_tm, ident.1 = "TM3", ident.2 = "TM")
combined_adult_wt_john_tm3_vs_tm
write.csv(combined_adult_wt_john_tm3_vs_tm, 'combined_adult_wt_john_tm3_vs_tm.csv')

#TM vs other POM
Idents(combined_adult_wt_john_pom) <- "CellType"
new.cluster.ids <- c('TM','TM','TM','other','other','other','other','other','other','other')
names(new.cluster.ids) <- levels(combined_adult_wt_john_pom)
combined_adult_wt_john_pom <- RenameIdents(combined_adult_wt_john_pom, new.cluster.ids)
combined_adult_wt_john_pom$CellTypeTM <- Idents(combined_adult_wt_john_pom)

Idents(combined_adult_wt_john_pom) <- "CellTypeTM"
combined_adult_wt_john_tm_vs_pom <- FindMarkers(combined_adult_wt_john_pom, ident.1 = "TM", ident.2 = "other")
combined_adult_wt_john_tm_vs_pom
write.csv(combined_adult_wt_john_tm_vs_pom, 'combined_adult_wt_john_tm_vs_pom.csv')

#### TM vs all other cells
Idents(combined_adult_wt_john) <- "CellType"
combined_adult_wt_john_no_pom <- subset(combined_adult_wt_john, idents = c("Epithelial","Iris_and_CB","Endothelial","Immune","Neuron"))

Idents(combined_adult_wt_john_no_pom) <- "CellType"
new.cluster.ids <- c("other","other","other","other","other")
names(new.cluster.ids) <- levels(combined_adult_wt_john_no_pom)
combined_adult_wt_john_no_pom <- RenameIdents(combined_adult_wt_john_no_pom, new.cluster.ids)
combined_adult_wt_john_no_pom$CellTypeTM <- Idents(combined_adult_wt_john_no_pom)

combined_adult_wt_john_tm_vs_all <- merge(combined_adult_wt_john_no_pom,combined_adult_wt_john_pom)
combined_adult_wt_john_tm_vs_all$CellTypeTM <- as.factor(combined_adult_wt_john_tm_vs_all$CellTypeTM)

Idents(combined_adult_wt_john_tm_vs_all) <- "CellTypeTM"
combined_adult_wt_john_tm_vs_all_deg <- FindMarkers(combined_adult_wt_john_tm_vs_all, ident.1 = "TM", ident.2 = "other")
combined_adult_wt_john_tm_vs_all_deg
write.csv(combined_adult_wt_john_tm_vs_all_deg, 'combined_adult_wt_john_tm_vs_all_deg.csv')

######### Specific TM cell comparisons
setwd("~/Desktop")

# All DEG excel files were modified
TM1_en <- fread("combined_adult_wt_john_tm_TM1_markers_markergenes_de.csv")

# DE pathways against all POM derived genes
TM1_en_pathway <- bitr(TM1_en$SYMBOL[which(TM1_en$abs_avg_log2FC > 0)], fromType = "SYMBOL",
                         toType = c("ENSEMBL", "ENTREZID"),
                         OrgDb = org.Mm.eg.db)
universe.df <- bitr(rownames(combined_adult_wt_john_tm), fromType = "SYMBOL",
                    toType = c("ENSEMBL", "ENTREZID"),
                    OrgDb = org.Mm.eg.db)
ego <- enrichGO(gene          = TM1_en_pathway$ENTREZID,
                universe      = universe.df$ENTREZID,
                OrgDb         = org.Mm.eg.db,
                ont           = "CC", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.01)

head(summary(ego))
TM1_en_pathway_cc <- as.data.frame(summary(ego))
write.csv(TM1_en_pathway_cc, "TM1_en_pathway_cc.csv")

ego <- enrichGO(gene          = TM1_en_pathway$ENTREZID,
                universe      = universe.df$ENTREZID,
                OrgDb         = org.Mm.eg.db,
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.01)

head(summary(ego), n=10)
TM1_en_pathway_bp <- as.data.frame(summary(ego))
write.csv(TM1_en_pathway_bp, "TM1_en_pathway_bp.csv")

ego <- enrichGO(gene          = TM1_en_pathway$ENTREZID,
                universe      = universe.df$ENTREZID,
                OrgDb         = org.Mm.eg.db,
                ont           = "MF", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.01)

head(summary(ego), n=30)
TM1_en_pathway_mf <- as.data.frame(summary(ego))
write.csv(TM1_en_pathway_mf, "TM1_en_pathway_mf.csv")

##### GSEA analysis
original_gene_list <- TM1_en$avg_log2FC
names(original_gene_list) <- TM1_en$SYMBOL
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

gse <- gseGO(geneList=gene_list, 
             ont ="CC", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Mm.eg.db, 
             pAdjustMethod = "BH")

dotplot(gse, showCategory=20, split=".sign") + facet_grid(.~.sign)

TM1_en_cc <- as.data.frame(summary(gse))
write.csv(TM1_en_cc, "TM1_en_cc.csv")

gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Mm.eg.db, 
             pAdjustMethod = "BH")

dotplot(gse, showCategory=20, split=".sign") + facet_grid(.~.sign)

TM1_en_bp <- as.data.frame(summary(gse))
write.csv(TM1_en_bp, "TM1_en_bp.csv")

gse <- gseGO(geneList=gene_list, 
             ont ="MF", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Mm.eg.db, 
             pAdjustMethod = "BH")

dotplot(gse, showCategory=20, split=".sign") + facet_grid(.~.sign)

TM1_en_mf <- as.data.frame(summary(gse))
write.csv(TM1_en_mf, "TM1_en_mf.csv")

######### Make barplot of results

setwd("~/Desktop")
data <- read.csv("combined_adult_wt_john_pom_TM_pathway_MF_all_tm.csv")
data <- read.csv("05302024_v265d_tm_wt_het_significant_pathways_with_gsea.csv")

p<-ggplot(data=data, aes(x=Group, y=Genes_in_pathway)) +
  geom_bar(aes(fill = log_p.adjust), stat = "identity")

#Adjust pathway IDs in scale_x_discrete for each individual analysis
#Representative pathways shown

p + scale_x_discrete(limits=c("TM3_CC_GO:0005746","TM3_CC_GO:0098803","TM3_CC_GO:0045271","TM3_CC_GO:0030964","TM3_CC_GO:0005747","TM3_CC_GO:0022626","TM3_CC_GO:0043209","TM3_CC_GO:0062023","TM3_CC_GO:0030312","TM3_CC_GO:0031012")) +
  scale_y_continuous(expand=c(0,0), limits=c(0,30)) +
  scale_color_gradient(
    limits = c(2,10),
    low = "blue",
    high = "red",
    space = "Lab",
    na.value = "red",
    guide = "colourbar",
    aesthetics = c("color", "fill")
  ) +
  coord_flip()

######## Plot TM marker location data

setwd("~/Desktop")

#### Data was generated in Imaris and entered into excel sheets
data <- read.csv("06172024_TM_cluster_low_resolution_marker_area_data.csv")

p<-ggplot(data=data, aes(x=Group, y=Size)) +
  geom_bar(aes(fill = pct_outer), stat = "identity")

#### Adjust order and identity of markers as desired
p + scale_x_discrete(limits=c("Chil1_Antibody","Myoc_Antibody","Myoc_InSitu","Crym_Antibody","Edn3_InSitu","Inmt_InSitu","aSMA_Antibody","aSMA_InSitu","Lypd1_InSitu","Tfap2b_Antibody","TM1_com","TM2_com","TM3_com","Test1","Test2")) +
  scale_y_continuous(expand=c(0,0), limits=c(0,10)) + 
  scale_fill_gradientn(colours = c("blue", "gold", "red"),
                       values = scales::rescale(c(0, 9, 18, 27, 35))) +
  coord_flip()

######### Analyze Stamer lab B6 data seprately

Idents(p90_b6_stamer_wt1) <- "seurat_clusters"
new.cluster.ids <- c("Fibro_1","Epi","Fibro_2","Immune_1","Immune_2","Endo","Neuron","Iris_CB_1","Iris_CB_2")
names(new.cluster.ids) <- levels(p90_b6_stamer_wt1)
p90_b6_stamer_wt1 <- RenameIdents(p90_b6_stamer_wt1, new.cluster.ids)
p90_b6_stamer_wt1$CellType <- Idents(p90_b6_stamer_wt1)
p90_b6_stamer_wt1$CellType <- factor(p90_b6_stamer_wt1$CellType, levels = c("Fibro_1","Fibro_2","Epi","Iris_CB_1","Iris_CB_2","Endo","Immune_1","Immune_2","Neuron"))

Idents(p90_b6_stamer_wt1) <- "CellType"
new.cluster.ids <- c("TM_containing","TM_containing","Epithelial","Iris_and_CB","Iris_and_CB","Endothelial","Immune","Immune","Iris_and_CB")
names(new.cluster.ids) <- levels(p90_b6_stamer_wt1)
p90_b6_stamer_wt1 <- RenameIdents(p90_b6_stamer_wt1, new.cluster.ids)
p90_b6_stamer_wt1$CellType <- Idents(p90_b6_stamer_wt1)
p90_b6_stamer_wt1$CellType <- factor(p90_b6_stamer_wt1$CellType, levels = c("TM_containing","Epithelial","Iris_and_CB","Endothelial","Immune"))

Idents(p90_b6_stamer_wt1) <- "CellType"
mycolors <- c("#8dd3c7","#80b1d3","#fdb462","#b3de69","#bebada")
DimPlot(p90_b6_stamer_wt1, reduction = "umap", cols = mycolors, label = TRUE)

Idents(p90_b6_stamer_wt1) <- "CellType"
markers_p90_b6_stamer_wt1 <- FindAllMarkers(p90_b6_stamer_wt1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <- markers_p90_b6_stamer_wt1 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top10 <- markers_p90_b6_stamer_wt1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
my_levels <- c("TM_containing","Epithelial","Iris_and_CB","Endothelial","Immune")
levels(p90_b6_stamer_wt1) <- my_levels
DoHeatmap(p90_b6_stamer_wt1, features = top5$gene, group.colors = mycolors) +NoLegend() +scale_fill_viridis()

####### Save all cells combined object
saveRDS(p90_b6_stamer_wt1, file = "~/Desktop/Lmx1b/Development/Single_Cell_Sequencing/RDS_files/12052023_p90_b6_stamer_all_cells.RDS")

####### load all cells combined object
p90_b6_stamer_wt1 <- readRDS("~/Desktop/Lmx1b/Development/Single_Cell_Sequencing/RDS_files/12052023_p90_b6_stamer_all_cells.RDS")

### Generate TM containing subset
Idents(p90_b6_stamer_wt1) <- "seurat_clusters"
p90_b6_stamer_pom <- subset(p90_b6_stamer_wt1, idents = c(0,2))
DefaultAssay(p90_b6_stamer_pom) <- "RNA"
p90_b6_stamer_pom <- NormalizeData(p90_b6_stamer_pom, normalization.method = "LogNormalize", scale.factor = 10000)
p90_b6_stamer_pom <- FindVariableFeatures(p90_b6_stamer_pom, nfeatures = 2000)
p90_b6_stamer_pom <- ScaleData(p90_b6_stamer_pom)
p90_b6_stamer_pom <- RunPCA(p90_b6_stamer_pom, features = VariableFeatures(object = p90_b6_stamer_pom))
p90_b6_stamer_pom <- FindNeighbors(p90_b6_stamer_pom, dims = 1:10)
p90_b6_stamer_pom <- FindClusters(p90_b6_stamer_pom, resolution = 0.3)
set.seed(500)
p90_b6_stamer_pom <- RunUMAP(p90_b6_stamer_pom, dims = 1:10)
nb.cols <- 30
mycolors <- colorRampPalette(brewer.pal(9, "Paired"))(nb.cols)
DimPlot(p90_b6_stamer_pom, reduction = "umap", cols = mycolors, raster = FALSE, label = TRUE)

Idents(p90_b6_stamer_pom) <- "seurat_clusters"
new.cluster.ids <- c("TM2","TM1","CM","TM3","Kera","CE","SF","Peri","IS","Schw")
names(new.cluster.ids) <- levels(p90_b6_stamer_pom)
p90_b6_stamer_pom <- RenameIdents(p90_b6_stamer_pom, new.cluster.ids)
p90_b6_stamer_pom$CellType <- Idents(p90_b6_stamer_pom)
p90_b6_stamer_pom$CellType <- factor(p90_b6_stamer_pom$CellType, levels = c('TM1','TM2','TM3','IS','SF','Kera','CE','CM','Peri','Schw'))
mycolors <- c("#a6cee3","#e31a1c","#33a02c","#fb9a99","khaki3","#b2df8a","#cab2d6","#ff7f00","#1f78b4","#6a3d9a")
Idents(p90_b6_stamer_pom) <- "CellType"
DimPlot(p90_b6_stamer_pom, reduction = "umap", cols = mycolors, raster = FALSE, label = TRUE)

####### Save all cells combined object
saveRDS(p90_b6_stamer_pom, file = "~/Desktop/Lmx1b/Development/Single_Cell_Sequencing/RDS_files/12052023_p90_b6_stamer_pom.RDS")

####### load all cells combined object
p90_b6_stamer_pom <- readRDS("~/Desktop/Lmx1b/Development/Single_Cell_Sequencing/RDS_files/12052023_p90_b6_stamer_pom.RDS")

Idents(p90_b6_stamer_pom) <- "CellType"
markers_p90_b6_stamer_pom <- FindAllMarkers(p90_b6_stamer_pom, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <- markers_p90_b6_stamer_pom %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top10 <- markers_p90_b6_stamer_pom %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
my_levels <- c('TM1','TM2','TM3','IS','SF','Kera','CE','CM','Peri','Schw')
levels(p90_b6_stamer_pom) <- my_levels
DoHeatmap(p90_b6_stamer_pom, features = top5$gene, group.colors = mycolors) +NoLegend() +scale_fill_viridis()

Idents(p90_b6_stamer_pom) <- "CellType"
p90_b6_stamer_tm <- subset(p90_b6_stamer_pom, idents = c('TM1','TM2','TM3'))

Idents(p90_b6_stamer_tm) <- "CellType"
markers_p90_b6_stamer_tm <- FindAllMarkers(p90_b6_stamer_tm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top20 <- markers_p90_b6_stamer_tm %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
my_levels <- c('TM1','TM2','TM3','IS','SF','Kera','CE','CM','Peri','Schw')
levels(p90_b6_stamer_tm) <- my_levels
DoHeatmap(p90_b6_stamer_tm, features = top20$gene, group.colors = mycolors) +NoLegend() +scale_fill_viridis()

Idents(p90_b6_stamer_pom) <- "CellType"
gene_list_plot <- c("Chil1","Tnmd","Cldn10","Myoc","Fmod","Pgf","Inmt","Cygb","Crym","Gpc3","Lypd1","Tfap2b","Gpr137b","Tagln","Acta2")
DotPlot(object = p90_b6_stamer_pom, features = gene_list_plot, idents = c('TM1','TM2','TM3'),cols = c("yellow","darkorchid3"))

######### Integrate John Lab B6 data with Stamer lab B6 data

anchors_adult_wt_b6_john_stamer <- FindIntegrationAnchors(object.list = list(p60_b6_john_wt1, p60_b6_john_wt2, p90_b6_stamer_wt1), dims = 1:10)
combined_adult_wt_b6_john_stamer <- IntegrateData(anchorset = anchors_adult_wt_b6_john_stamer, dims = 1:10)
DefaultAssay(combined_adult_wt_b6_john_stamer) <- "integrated"
combined_adult_wt_b6_john_stamer <- ScaleData(combined_adult_wt_b6_john_stamer, verbose = FALSE)
combined_adult_wt_b6_john_stamer <- RunPCA(combined_adult_wt_b6_john_stamer, npcs = 30, verbose = FALSE)
combined_adult_wt_b6_john_stamer <- RunUMAP(combined_adult_wt_b6_john_stamer, reduction = "pca", dims = 1:10, metric = 'correlation')
combined_adult_wt_b6_john_stamer <- FindNeighbors(combined_adult_wt_b6_john_stamer, reduction = "pca", dims = 1:10)
combined_adult_wt_b6_john_stamer <- FindClusters(combined_adult_wt_b6_john_stamer, resolution = 0.15)
combined_adult_wt_b6_john_stamer$seurat_clusters <- Idents(combined_adult_wt_b6_john_stamer)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
p1 <- DimPlot(combined_adult_wt_b6_john_stamer, reduction = "umap", group.by = "orig.ident", cols = brewer.pal(n=4,name="Set2"))
p2 <- DimPlot(combined_adult_wt_b6_john_stamer, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
nb.cols <- 11
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
mycolors <- c("#8dd3c7","#80b1d3","#ffffb3","#fccde5","#b3de69","#bebada","#fdb462","#fdb462","#fccde5","#fb8072","#d9d9d9")
DimPlot(combined_adult_wt_b6_john_stamer, reduction = "umap", cols = mycolors, label = TRUE)
DimPlot(combined_adult_wt_b6_john_stamer, reduction = "umap", split.by = "Lab", cols = mycolors, label = TRUE)

####### Save all cells combined object
saveRDS(combined_adult_wt_b6_john_stamer, file = "~/Desktop/Lmx1b/Development/Single_Cell_Sequencing/RDS_files/11202023_combined_adult_wt_b6_john_stamer.RDS")

####### load all cells combined object
combined_adult_wt_b6_john_stamer <- readRDS("~/Desktop/Lmx1b/Development/Single_Cell_Sequencing/RDS_files/11202023_combined_adult_wt_b6_john_stamer.RDS")

DefaultAssay(combined_adult_wt_b6_john_stamer) <- "RNA"
Idents(combined_adult_wt_b6_john_stamer) <- "seurat_clusters"
combined_adult_wt_b6_john_stamer_pom <- subset(combined_adult_wt_b6_john_stamer, idents = c(0,2))
DefaultAssay(combined_adult_wt_b6_john_stamer_pom) <- "RNA"
combined_adult_wt_b6_john_stamer_pom <- NormalizeData(combined_adult_wt_b6_john_stamer_pom, normalization.method = "LogNormalize", scale.factor = 10000)
combined_adult_wt_b6_john_stamer_pom <- FindVariableFeatures(combined_adult_wt_b6_john_stamer_pom, nfeatures = 2000)
combined_adult_wt_b6_john_stamer_pom <- ScaleData(combined_adult_wt_b6_john_stamer_pom)
combined_adult_wt_b6_john_stamer_pom <- RunPCA(combined_adult_wt_b6_john_stamer_pom, features = VariableFeatures(object = combined_adult_wt_b6_john_stamer_pom))
combined_adult_wt_b6_john_stamer_pom <- FindNeighbors(combined_adult_wt_b6_john_stamer_pom, dims = 1:5)
combined_adult_wt_b6_john_stamer_pom <- FindClusters(combined_adult_wt_b6_john_stamer_pom, resolution = 0.5)
set.seed(500)
combined_adult_wt_b6_john_stamer_pom <- RunUMAP(combined_adult_wt_b6_john_stamer_pom, dims = 1:5)
nb.cols <- 30
mycolors <- colorRampPalette(brewer.pal(9, "Paired"))(nb.cols)
DimPlot(combined_adult_wt_b6_john_stamer_pom, reduction = "umap", split.by = "Lab", cols = mycolors, raster = FALSE, label = TRUE)
DimPlot(combined_adult_wt_b6_john_stamer_pom, reduction = "umap", cols = mycolors, raster = FALSE, label = TRUE)

new.cluster.ids <- c("TM1","TM3","TM2","TM2","CM","TM1","IS","TM2","Kera","TM3","TM1","Kera","SF","Schw","Peri","Kera","CE")
names(new.cluster.ids) <- levels(combined_adult_wt_b6_john_stamer_pom)
combined_adult_wt_b6_john_stamer_pom <- RenameIdents(combined_adult_wt_b6_john_stamer_pom, new.cluster.ids)
combined_adult_wt_b6_john_stamer_pom$CellType <- Idents(combined_adult_wt_b6_john_stamer_pom)
combined_adult_wt_b6_john_stamer_pom$CellType <- factor(combined_adult_wt_b6_john_stamer_pom$CellType, levels = c('TM1','TM2','TM3','IS','SF','Kera','CE','CM','Peri','Schw'))

Idents(combined_adult_wt_b6_john_stamer_pom) <- "CellType"
mycolors <- c("#a6cee3","#e31a1c","#33a02c","#fb9a99","khaki3","#b2df8a","#cab2d6","#ff7f00","#1f78b4","#6a3d9a")
DimPlot(combined_adult_wt_b6_john_stamer_pom, reduction = "umap", cols = mycolors)
DimPlot(combined_adult_wt_b6_john_stamer_pom, reduction = "umap", cols = mycolors, split.by = "Lab")

####### Save all cells combined object
saveRDS(combined_adult_wt_b6_john_stamer_pom, file = "~/Desktop/Lmx1b/Development/Single_Cell_Sequencing/RDS_files/11202023_combined_adult_wt_b6_john_stamer_pom.RDS")

####### load all cells combined object
combined_adult_wt_b6_john_stamer_pom <- readRDS("~/Desktop/Lmx1b/Development/Single_Cell_Sequencing/RDS_files/11202023_combined_adult_wt_b6_john_stamer_pom.RDS")

####### Multiome data - get additional details from Taibo

####### load all cells combined object
P60_wt_single_nuc_atac <- readRDS("~/Desktop/Lmx1b/Development/Single_Cell_Sequencing/RDS_files/12042023_P60_wt_single_nuc_atac.RDS")

#Subset POM-derived Pitx2+ cluster
Idents(P60_wt_single_nuc_atac) <- "CellType"
P60_wt_single_nuc_atac_pom <- subset(P60_wt_single_nuc_atac, idents = c("Fibro_1","Fibro_2"))

#Cluster all RNA cells
DefaultAssay(P60_wt_single_nuc_atac_pom) <- "RNA"
P60_wt_single_nuc_atac_pom <- NormalizeData(P60_wt_single_nuc_atac_pom, normalization.method = "LogNormalize", scale.factor = 10000)
P60_wt_single_nuc_atac_pom <- FindVariableFeatures(P60_wt_single_nuc_atac_pom, nfeatures = 2000)
P60_wt_single_nuc_atac_pom <- ScaleData(P60_wt_single_nuc_atac_pom)
P60_wt_single_nuc_atac_pom <- RunPCA(P60_wt_single_nuc_atac_pom, features = VariableFeatures(object = P60_wt_single_nuc_atac_pom))
P60_wt_single_nuc_atac_pom <- FindNeighbors(P60_wt_single_nuc_atac_pom, dims = 1:10)
P60_wt_single_nuc_atac_pom <- FindClusters(P60_wt_single_nuc_atac_pom, resolution = 0.4)
#P60_wt_single_nuc_atac_pom <- FindClusters(P60_wt_single_nuc_atac_pom, resolution = 0.12) 
set.seed(500)
P60_wt_single_nuc_atac_pom <- RunUMAP(P60_wt_single_nuc_atac_pom, dims = 1:10)
nb.cols <- 12
mycolors <- colorRampPalette(brewer.pal(9, "Paired"))(nb.cols)
DimPlot(P60_wt_single_nuc_atac_pom, reduction = "umap", cols = mycolors, raster = FALSE, label = TRUE)

Idents(P60_wt_single_nuc_atac_pom) <- "seurat_clusters"
new.cluster.ids <- c("TM3","TM2","IS","LE","TM1","TM2","CE","CM","SF","Peri","Kera")
names(new.cluster.ids) <- levels(P60_wt_single_nuc_atac_pom)
P60_wt_single_nuc_atac_pom <- RenameIdents(P60_wt_single_nuc_atac_pom, new.cluster.ids)
P60_wt_single_nuc_atac_pom$CellType <- Idents(P60_wt_single_nuc_atac_pom)
P60_wt_single_nuc_atac_pom$CellType <- factor(P60_wt_single_nuc_atac_pom$CellType, levels = c('TM1','TM2','TM3','IS','SF','Kera','CE',"LE",'CM','Peri'))

Idents(P60_wt_single_nuc_atac_pom) <- "CellType"
mycolors <- c("#a6cee3","#e31a1c","#33a02c","#fb9a99","khaki3","#b2df8a","#cab2d6","violetred4","#ff7f00","#1f78b4")
DimPlot(P60_wt_single_nuc_atac_pom, reduction = "umap", cols = mycolors, label = "TRUE")

Idents(P60_wt_single_nuc_atac_pom) <- "CellType"
P60_wt_single_nuc_atac_tm <- subset(P60_wt_single_nuc_atac_pom, idents = c("TM1","TM2","TM3"))

Idents(P60_wt_single_nuc_atac_tm) <- "CellType"
markers_P60_wt_single_nuc_atac_tm <- FindAllMarkers(P60_wt_single_nuc_atac_tm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top20 <- markers_P60_wt_single_nuc_atac_tm %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
my_levels <- c('TM1','TM2','TM3')
levels(P60_wt_single_nuc_atac_tm) <- my_levels
DoHeatmap(P60_wt_single_nuc_atac_tm, features = top20$gene, group.colors = mycolors) +NoLegend() +scale_fill_viridis()

####### Save all cells combined object
saveRDS(P60_wt_single_nuc_atac_pom, file = "~/Desktop/Lmx1b/Development/Single_Cell_Sequencing/RDS_files/12072023_P60_wt_single_nuc_atac_pom.RDS")

####### load all cells combined object
P60_wt_single_nuc_atac_pom <- readRDS("~/Desktop/Lmx1b/Development/Single_Cell_Sequencing/RDS_files/12072023_P60_wt_single_nuc_atac_pom.RDS")

########### Making endothelial object

Idents(combined_adult_wt_john) <- "CellType"
combined_adult_wt_john_sc <- subset(combined_adult_wt_john, idents = c('Endothelial'))

DefaultAssay(combined_adult_wt_john_sc) <- "RNA"
combined_adult_wt_john_sc <- NormalizeData(combined_adult_wt_john_sc, normalization.method = "LogNormalize", scale.factor = 10000)
combined_adult_wt_john_sc <- FindVariableFeatures(combined_adult_wt_john_sc, nfeatures = 2000)
combined_adult_wt_john_sc <- ScaleData(combined_adult_wt_john_sc)
combined_adult_wt_john_sc <- RunPCA(combined_adult_wt_john_sc, features = VariableFeatures(object = combined_adult_wt_john_sc))
combined_adult_wt_john_sc <- FindNeighbors(combined_adult_wt_john_sc, dims = 1:10)
combined_adult_wt_john_sc <- FindClusters(combined_adult_wt_john_sc, resolution = 0.5)
set.seed(500)
combined_adult_wt_john_sc <- RunUMAP(combined_adult_wt_john_sc, dims = 1:10)
nb.cols <- 6
mycolors <- colorRampPalette(brewer.pal(9, "Paired"))(nb.cols)
DimPlot(combined_adult_wt_john_sc, reduction = "umap", cols = mycolors, raster = FALSE, label = TRUE)

new.cluster.ids <- c("SEC","BEC","BEC","BEC","CC","LEC")
names(new.cluster.ids) <- levels(combined_adult_wt_john_sc)
combined_adult_wt_john_sc <- RenameIdents(combined_adult_wt_john_sc, new.cluster.ids)
combined_adult_wt_john_sc$CellType <- Idents(combined_adult_wt_john_sc)

######## Disease score enrichment

combined_adult_wt_john_pom$CellType1 <- combined_adult_wt_john_pom$CellType
combined_adult_wt_john_sc$CellType1 <- combined_adult_wt_john_sc$CellType
combined_adult_wt_john$CellType1 <- combined_adult_wt_john$CellType

Idents(combined_adult_wt_john) <- "CellType1"
combined_adult_wt_john_no_pom_sc <- subset(combined_adult_wt_john, idents = c("Epithelial","Iris_and_CB","Immune","Neuron"))
combined_adult_wt_john_disease_analysis_1 <- merge(combined_adult_wt_john_no_pom_sc,combined_adult_wt_john_pom)
combined_adult_wt_john_disease_analysis <- merge(combined_adult_wt_john_disease_analysis_1,combined_adult_wt_john_sc)
combined_adult_wt_john_disease_analysis$CellType1 <- as.factor(combined_adult_wt_john_disease_analysis$CellType1)
combined_adult_wt_john_disease_analysis$CellType1 <- factor(combined_adult_wt_john_disease_analysis$CellType1, levels = c("TM1","TM2","TM3","CE","Kera","SF","CM","IS","Peri","Schw","SEC","BEC","CC","LEC","Epithelial","Iris_and_CB","Immune","Neuron"))

iop_gwas_genes <- read.csv("05212024_iop_gwas_gene_list.csv", header = TRUE)
Idents(combined_adult_wt_john_disease_analysis) <- "CellType1"
DefaultAssay(combined_adult_wt_john_disease_analysis) <- "RNA"

######
disease_score_iop_short <- AddModuleScore(combined_adult_wt_john_disease_analysis, features = list(iop_gwas_genes$IOP.loci.short), name = "IOP_Disease_score_short")
disease_score_poag_short <- AddModuleScore(combined_adult_wt_john_disease_analysis, features = list(iop_gwas_genes$POAG.loci.short), name = "POAG_Disease_score_short")

### Combine into one object
gwas_combined <- disease_score_iop_short
gwas_combined$POAG_Disease_score_short1 <- disease_score_poag_short$POAG_Disease_score_short1

Idents(gwas_combined) <- "CellType1"
DefaultAssay(gwas_combined) <- "RNA"
DotPlot(object = gwas_combined, features = c('IOP_Disease_score_short1','POAG_Disease_score_short1'),cols = c("yellow","darkorchid3"))

##### extracellular matrix organization BP
pathway_genes <- read.csv("extracellular matrix organization BP.csv", header = TRUE)
disease_score_pathway <- AddModuleScore(combined_adult_wt_john_disease_analysis, features = list(pathway_genes$All_pathway), name = "All_pathway")

Idents(disease_score_pathway) <- "CellType"
DotPlot(object = disease_score_pathway, features = 'All_pathway1')

disease_score_pathway <- AddModuleScore(combined_adult_wt_john_disease_analysis, features = list(pathway_genes$GWAS_overlap), name = "GWAS_overlap")

Idents(disease_score_pathway) <- "CellType"
DotPlot(object = disease_score_pathway, features = 'GWAS_overlap1',cols = c("yellow","darkorchid3"))

##### actin binding MF
pathway_genes <- read.csv("actin binding MF.csv", header = TRUE)
disease_score_pathway <- AddModuleScore(combined_adult_wt_john_disease_analysis, features = list(pathway_genes$All_pathway), name = "All_pathway")

Idents(disease_score_pathway) <- "CellType1"
DotPlot(object = disease_score_pathway, features = 'All_pathway1')

disease_score_pathway <- AddModuleScore(combined_adult_wt_john_disease_analysis, features = list(pathway_genes$GWAS_overlap), name = "GWAS_overlap")

Idents(disease_score_pathway) <- "CellType"
DotPlot(object = disease_score_pathway, features = 'GWAS_overlap1')

##### cell-substrate adhesion BP
pathway_genes <- read.csv("cell-substrate adhesion BP.csv", header = TRUE)
disease_score_pathway <- AddModuleScore(combined_adult_wt_john_disease_analysis, features = list(pathway_genes$All_pathway), name = "All_pathway")

Idents(disease_score_pathway) <- "CellType"
DotPlot(object = disease_score_pathway, features = 'All_pathway1')

disease_score_pathway <- AddModuleScore(combined_adult_wt_john_disease_analysis, features = list(pathway_genes$GWAS_overlap), name = "GWAS_overlap")

Idents(disease_score_pathway) <- "CellType"
DotPlot(object = disease_score_pathway, features = 'GWAS_overlap1',cols = c("yellow","darkorchid3"))

##### regulation of Rho protein signal transduction BP
pathway_genes <- read.csv("regulation of Rho protein signal transduction BP.csv", header = TRUE)
disease_score_pathway <- AddModuleScore(combined_adult_wt_john_disease_analysis, features = list(pathway_genes$All_pathway), name = "All_pathway")

Idents(disease_score_pathway) <- "CellType"
DotPlot(object = disease_score_pathway, features = 'All_pathway1')

disease_score_pathway <- AddModuleScore(combined_adult_wt_john_disease_analysis, features = list(pathway_genes$GWAS_overlap), name = "GWAS_overlap")

Idents(disease_score_pathway) <- "CellType"
DotPlot(object = disease_score_pathway, features = 'GWAS_overlap1',cols = c("yellow","darkorchid3"))

##### cellular response to hormone stimulus BP
pathway_genes <- read.csv("cellular response to hormone stimulus BP.csv", header = TRUE)
disease_score_pathway <- AddModuleScore(combined_adult_wt_john_disease_analysis, features = list(pathway_genes$All_pathway), name = "All_pathway")

Idents(disease_score_pathway) <- "CellType"
DotPlot(object = disease_score_pathway, features = 'All_pathway1')

disease_score_pathway <- AddModuleScore(combined_adult_wt_john_disease_analysis, features = list(pathway_genes$GWAS_overlap), name = "GWAS_overlap")

Idents(disease_score_pathway) <- "CellType"
DotPlot(object = disease_score_pathway, features = 'GWAS_overlap1')

##### regulation of cellular response to growth factor stimulus BP
pathway_genes <- read.csv("regulation of cellular response to growth factor stimulus BP.csv", header = TRUE)
disease_score_pathway <- AddModuleScore(combined_adult_wt_john_disease_analysis, features = list(pathway_genes$All_pathway), name = "All_pathway")

Idents(disease_score_pathway) <- "CellType"
DotPlot(object = disease_score_pathway, features = 'All_pathway1')

disease_score_pathway <- AddModuleScore(combined_adult_wt_john_disease_analysis, features = list(pathway_genes$GWAS_overlap), name = "GWAS_overlap")

Idents(disease_score_pathway) <- "CellType"
DotPlot(object = disease_score_pathway, features = 'GWAS_overlap1')

##### vascular endothelial growth factor receptor signaling pathway BP
pathway_genes <- read.csv("vascular endothelial growth factor receptor signaling pathway BP.csv", header = TRUE)
disease_score_pathway <- AddModuleScore(combined_adult_wt_john_disease_analysis, features = list(pathway_genes$All_pathway), name = "All_pathway")

Idents(disease_score_pathway) <- "CellType"
DotPlot(object = disease_score_pathway, features = 'All_pathway1')

disease_score_pathway <- AddModuleScore(combined_adult_wt_john_disease_analysis, features = list(pathway_genes$GWAS_overlap), name = "GWAS_overlap")

Idents(disease_score_pathway) <- "CellType"
DotPlot(object = disease_score_pathway, features = 'GWAS_overlap1',cols = c("yellow","darkorchid3"))

##### Combine into same dotplot

pathway_genes <- read.csv("vascular endothelial growth factor receptor signaling pathway BP.csv", header = TRUE)
disease_score_pathway <- AddModuleScore(combined_adult_wt_john_disease_analysis, features = list(pathway_genes$All_pathway), name = "vegf_bp")

pathway_genes <- read.csv("regulation of Rho protein signal transduction BP.csv", header = TRUE)
disease_score_pathway1 <- AddModuleScore(combined_adult_wt_john_disease_analysis, features = list(pathway_genes$All_pathway), name = "rho_bp")

pathway_genes <- read.csv("cell-substrate adhesion BP.csv", header = TRUE)
disease_score_pathway2 <- AddModuleScore(combined_adult_wt_john_disease_analysis, features = list(pathway_genes$All_pathway), name = "cell_sub_bp")

pathway_genes <- read.csv("extracellular matrix organization BP.csv", header = TRUE)
disease_score_pathway3 <- AddModuleScore(combined_adult_wt_john_disease_analysis, features = list(pathway_genes$All_pathway), name = "ecm_bp")

pathway_genes <- read.csv("actin binding MF.csv", header = TRUE)
disease_score_pathway4 <- AddModuleScore(combined_adult_wt_john_disease_analysis, features = list(pathway_genes$All_pathway), name = "actin_mf")

pathway_genes <- read.csv("regulation of cellular response to growth factor stimulus BP.csv", header = TRUE)
disease_score_pathway5 <- AddModuleScore(combined_adult_wt_john_disease_analysis, features = list(pathway_genes$All_pathway), name = "gf_stim_bp")

pathway_genes <- read.csv("cellular response to hormone stimulus BP.csv", header = TRUE)
disease_score_pathway6 <- AddModuleScore(combined_adult_wt_john_disease_analysis, features = list(pathway_genes$All_pathway), name = "hor_stim_bp")

disease_score_pathway$rho_bp1 <- disease_score_pathway1$rho_bp1
disease_score_pathway$cell_sub_bp1 <- disease_score_pathway2$cell_sub_bp1
disease_score_pathway$ecm_bp1 <- disease_score_pathway3$ecm_bp1
disease_score_pathway$actin_mf1 <- disease_score_pathway4$actin_mf1
disease_score_pathway$gf_stim_bp1 <- disease_score_pathway5$gf_stim_bp1
disease_score_pathway$hor_stim_bp1 <- disease_score_pathway6$hor_stim_bp1
disease_score_pathway$gtpase_bp1 <- disease_score_pathway7$gtpase_bp1

Idents(disease_score_pathway) <- "CellType1"
DotPlot(object = disease_score_pathway, features = c('ecm_bp1','cell_sub_bp1','gf_stim_bp1','hor_stim_bp1','vegf_bp1','rho_bp1','actin_mf1'),cols = c("yellow","darkorchid3"))

#Mitochondrial analysis
lipid.mito.metab <- AddModuleScore(combined_adult_wt_john_disease_analysis, features = list(iop_gwas_genes$lipid.mito.metab), name = "lipid_mito_score")

Idents(lipid.mito.metab) <- "CellType1"
DotPlot(object = lipid.mito.metab, features = 'lipid_mito_score1',cols = c("yellow","darkorchid3"))

##
carb.mito.metab <- AddModuleScore(combined_adult_wt_john_disease_analysis, features = list(iop_gwas_genes$carb.mito.metab), name = "carb_mito_score")

Idents(carb.mito.metab) <- "CellType1"
DotPlot(object = carb.mito.metab, features = 'carb_mito_score1',cols = c("yellow","darkorchid3"))

### Combine into one object
mito_combined <- lipid.mito.metab
mito_combined$carb_mito_score1 <- carb.mito.metab$carb_mito_score1

DotPlot(object = mito_combined, features = c('lipid_mito_score1','carb_mito_score1'),cols = c("yellow","darkorchid3"))

######## Nichenet
######## This is one example with SECs as receiver cell 
######## For analyses with different sender and receiver cells, switch cell types in this section and rerun
######## Seurat object include all cell types here to get more accurate representation of enriched genes in receiver type

ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))

### Convert to mouse symbols
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols() 
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols() 

ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

dim(ligand_target_matrix)

lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))

weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
weighted_networks_lr = weighted_networks_lr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()

#### receiver
#### SEC
receiver = "SEC"
expressed_genes_receiver <- get_expressed_genes(receiver, combined_adult_wt_john_disease_analysis, pct = 0.1)

all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)

potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()

#### sender
sender_celltypes <- c("TM1", "TM2", "TM3")

# Use lapply to get the expressed genes of every sender cell type separately here
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, combined_adult_wt_john_disease_analysis, 0.1)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

length(expressed_genes_sender)
length(potential_ligands)
length(potential_ligands_focused)

####### Select enriched geneset for receiver cell type
# SEC
Idents(combined_adult_wt_john_disease_analysis) <- "CellType"
combined_adult_wt_john_disease_analysis_sec <- FindMarkers(combined_adult_wt_john_disease_analysis, ident.1 = "SEC", min.pct = 0.1, only.pos = TRUE)
combined_adult_wt_john_disease_analysis_sec

geneset_oi <- combined_adult_wt_john_disease_analysis_sec %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25)
geneset_oi <- rownames(geneset_oi)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

######## Background geneset
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

length(background_expressed_genes)
length(geneset_oi)

######## Select top ligands
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(dplyr::desc(aupr_corrected)))
ligand_activities

p_hist_lig_activity <- ggplot(ligand_activities, aes(x=aupr_corrected)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(30, aupr_corrected) %>% pull(aupr_corrected))),
             color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()

p_hist_lig_activity

best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)

vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% dplyr::select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)

(make_heatmap_ggplot(vis_ligand_aupr,
                     "Prioritized ligands", "Ligand activity", 
                     legend_title = "AUPR", color = "darkorange") + 
    theme(axis.text.x.top = element_blank()))  

###### Target gene inference
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()

nrow(active_ligand_target_links_df)
head(active_ligand_target_links_df)

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33) 

nrow(active_ligand_target_links)
head(active_ligand_target_links)

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

###### Make dotplot of ligand receptor interactions
make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                    color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

p_dotplot <- DotPlot(subset(combined_adult_wt_john_disease_analysis, CellType %in% sender_celltypes),
                     features = rev(best_upstream_ligands), cols = "RdYlBu") + 
  coord_flip() +
  scale_y_discrete(position = "right")

p_dotplot

###### Circos plots build up

######## The seurat object can't have more clusters than you are including in the plot, subset here-
Idents(combined_adult_wt_john_disease_analysis) <- "CellType"
combined_adult_wt_john_nichenet <- subset(combined_adult_wt_john_disease_analysis, idents = c("TM1","TM2","TM3","SEC"))

#### Define set of potential ligands 
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

########## NicheNet ligand activity analysis 
ligand_activities = predict_ligand_activities(geneset = geneset_oi, 
                                              background_expressed_genes = background_expressed_genes, 
                                              ligand_target_matrix = ligand_target_matrix, 
                                              potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(dplyr::desc(pearson)))

#############################
## Circos plot
#############################

top_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()

DefaultAssay(combined_adult_wt_john_nichenet) <- "RNA"
avg_expression_ligands = AverageExpression(combined_adult_wt_john_nichenet, features = top_ligands)

## assign ligands to sender cells
sender_ligand_assignment = avg_expression_ligands$RNA %>% apply(1, function(ligand_expression){
  ligand_expression > (ligand_expression %>% mean() + ligand_expression %>% sd())
}) %>% t()
sender_ligand_assignment = sender_ligand_assignment %>% apply(2, function(x){x[x == TRUE]}) %>% purrr::keep(function(x){length(x) > 0})
names(sender_ligand_assignment)

###############################

tmp <- list("TM1" = sender_ligand_assignment$TM1,
            "TM2" = sender_ligand_assignment$TM2,
            "TM3" = sender_ligand_assignment$TM3)
all_assigned_ligands = tmp %>% lapply(function(x){names(x)}) %>% unlist()

unique_ligands = all_assigned_ligands %>% table() %>% .[. == 1] %>% names()

ligand_type_indication_df = tibble(
  ligand_type = c(rep("TM1", times = tmp$TM1 %>% length()),
                  rep("TM2", times = tmp$TM2 %>% length()),
                  rep("TM3", times = tmp$TM3 %>% length())),
  ligand = c(names(tmp$TM1), names(tmp$TM2), names(tmp$TM3)))

################### Define ligand-target links of interest 

active_ligand_target_links_df = top_ligands %>% 
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi, 
         ligand_target_matrix = ligand_target_matrix, 
         n = 200) %>% 
  bind_rows() %>% 
  drop_na()

active_ligand_target_links_df = active_ligand_target_links_df %>% mutate(target_type = "TM_SEC") %>% inner_join(ligand_type_indication_df) 

### Adjust TM cell assignments based on RNA examination
### Skip this step for code example
write.csv(active_ligand_target_links_df, "TM_assignments.csv")
active_ligand_target_links_df <- read.csv("TM_assignments.csv", header = TRUE)

### Start here again
cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% 
  quantile(0.75)

active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)

ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())

circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)

### prepare the circos visualization 
### give each segment of ligands and targets a specific color and order 

grid_col_ligand =c("TM1" = "#a6cee3",
                   "TM2" = "#e31a1c",
                   "TM3" = "#33a02c"
)

grid_col_target =c("TM_SEC" = "pink")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

circos_links = circos_links %>% dplyr::mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% dplyr::inner_join(grid_col_tbl_ligand) %>% dplyr::inner_join(grid_col_tbl_target)
links_circle = circos_links %>% dplyr::select(ligand,target, weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(target,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$target)

grid_col =c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparant ~ ligand-target potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 

target_order = circos_links$target %>% unique()
ligand_order = c("Crlf1","Tgfb2","Fn1","Tgfb3","Vegfa","Edn3","Angpt1","Vegfc","Cxcl12") %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,target_order)

### define gaps between different segments

width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_target = 15
width_same_cell_same_target_type = 0.5

gaps = c(
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "TM1") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "TM2") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "TM3") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_type == "TM_SEC") %>% distinct(target) %>% nrow() -1)),
  width_ligand_target
)

#### Rendering 

circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = 0, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #

legend("topright",
       legend = c("TM1", "TM2", "TM3",  "SEC"),
       fill = c("#a6cee3","#e31a1c","#33a02c","pink"),       
       border = "black")

########## dotplot of ligands in TM cell subtypes

Idents(combined_adult_wt_john_tm) <- "CellType"
gene_list_plot <- c("Crlf1","Tgfb2","Fn1","Tgfb3","Vegfa","Edn3","Angpt1","Vegfc","Cxcl12")
DotPlot(object = combined_adult_wt_john_tm, features = gene_list_plot,cols = c("yellow","darkorchid3"))

############ B6 V265D P60 analysis

p60_het1 <- Read10X(data.dir = "~/Desktop/Lmx1b/Development/Single_Cell_Sequencing/SN074_Het/SN074_cellranger_count_outs/filtered_feature_bc_matrix")
p60_het1 <- CreateSeuratObject(counts = p60_het1, project = "p60_het1", min.cells = 3, min.features = 200)

p60_het2 <- Read10X(data.dir = "~/Desktop/Lmx1b/Development/Single_Cell_Sequencing/SN082_Het/SN082_cellranger_count_outs/filtered_feature_bc_matrix")
p60_het2 <- CreateSeuratObject(counts = p60_het2, project = "p60_het2", min.cells = 3, min.features = 200)

p60_het1[["percent.mt"]] <- PercentageFeatureSet(p60_het1, pattern = "^mt-")
VlnPlot(p60_het1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p60_het1 <- subset( p60_het1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
p60_het1 <- NormalizeData( p60_het1, normalization.method = "LogNormalize", scale.factor = 10000)
p60_het1 <- FindVariableFeatures( p60_het1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames( p60_het1)
p60_het1 <- ScaleData( p60_het1, features = all.genes)
p60_het1 <- RunPCA( p60_het1, features = VariableFeatures(object = p60_het1))
p60_het1<- FindNeighbors(p60_het1, dims = 1:10)
p60_het1 <- FindClusters(p60_het1, resolution = 0.1)
p60_het1<- RunUMAP( p60_het1, dims = 1:10)
DimPlot(p60_het1, reduction = "umap")
table(Idents( p60_het1))
p60_het1 <- AddMetaData(object = p60_het1, metadata = "B6", col.name = "Strain")
p60_het1 <- AddMetaData(object = p60_het1, metadata = "het", col.name = "Genotype")

p60_het2[["percent.mt"]] <- PercentageFeatureSet(p60_het2, pattern = "^mt-")
VlnPlot(p60_het2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p60_het2 <- subset( p60_het2, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
p60_het2 <- NormalizeData( p60_het2, normalization.method = "LogNormalize", scale.factor = 10000)
p60_het2 <- FindVariableFeatures( p60_het2, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames( p60_het2)
p60_het2 <- ScaleData( p60_het2, features = all.genes)
p60_het2 <- RunPCA( p60_het2, features = VariableFeatures(object = p60_het2))
p60_het2<- FindNeighbors(p60_het2, dims = 1:10)
p60_het2 <- FindClusters(p60_het2, resolution = 0.1)
p60_het2<- RunUMAP( p60_het2, dims = 1:10)
DimPlot(p60_het2, reduction = "umap")
table(Idents( p60_het2))
p60_het2 <- AddMetaData(object = p60_het2, metadata = "B6", col.name = "Strain")
p60_het2 <- AddMetaData(object = p60_het2, metadata = "het", col.name = "Genotype")

####### Integrate with WT B6 data

anchors_adult_wt_het_p60_john <- FindIntegrationAnchors(object.list = list(p60_b6_john_wt1, p60_b6_john_wt2, p60_het1, p60_het2), dims = 1:10)
combined_adult_wt_het_p60_john <- IntegrateData(anchorset = anchors_adult_wt_het_p60_john, dims = 1:10)
DefaultAssay(combined_adult_wt_het_p60_john) <- "integrated"
combined_adult_wt_het_p60_john <- ScaleData(combined_adult_wt_het_p60_john, verbose = FALSE)
combined_adult_wt_het_p60_john <- RunPCA(combined_adult_wt_het_p60_john, npcs = 30, verbose = FALSE)
combined_adult_wt_het_p60_john <- RunUMAP(combined_adult_wt_het_p60_john, reduction = "pca", dims = 1:10, metric = 'correlation')
combined_adult_wt_het_p60_john <- FindNeighbors(combined_adult_wt_het_p60_john, reduction = "pca", dims = 1:10)
combined_adult_wt_het_p60_john <- FindClusters(combined_adult_wt_het_p60_john, resolution = 0.1)
combined_adult_wt_het_p60_john$seurat_clusters <- Idents(combined_adult_wt_het_p60_john)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
p1 <- DimPlot(combined_adult_wt_het_p60_john, reduction = "umap", group.by = "orig.ident", cols = brewer.pal(n=4,name="Set2"))
p2 <- DimPlot(combined_adult_wt_het_p60_john, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
nb.cols <- 15
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
DimPlot(combined_adult_wt_het_p60_john, reduction = "umap", cols = mycolors, label = TRUE)
DimPlot(combined_adult_wt_het_p60_john, reduction = "umap", cols = mycolors, split.by = "orig.ident")

Idents(combined_adult_wt_het_p60_john) <- "seurat_clusters"
new.cluster.ids <- c("TM_containing","Epithelial","TM_containing","Epithelial","Endothelial","Immune","Iris_and_CB","Epithelial","Neuron")
names(new.cluster.ids) <- levels(combined_adult_wt_het_p60_john)
combined_adult_wt_het_p60_john <- RenameIdents(combined_adult_wt_het_p60_john, new.cluster.ids)
combined_adult_wt_het_p60_john$CellType <- Idents(combined_adult_wt_het_p60_john)
combined_adult_wt_het_p60_john$CellType <- factor(combined_adult_wt_het_p60_john$CellType, levels = c("TM_containing","Epithelial","Iris_and_CB","Endothelial","Immune","Neuron"))
mycolors <- c("#8dd3c7","#80b1d3","#fdb462","#b3de69","#bebada","#fb8072")
Idents(combined_adult_wt_het_p60_john) <- "CellType"
DimPlot(combined_adult_wt_het_p60_john, reduction = "umap", cols = mycolors, label = TRUE)

############ Making POM-derived object
DefaultAssay(combined_adult_wt_het_p60_john) <- "RNA"
Idents(combined_adult_wt_het_p60_john) <- "seurat_clusters"
combined_adult_wt_het_p60_john_pom <- subset(combined_adult_wt_het_p60_john, idents = c("0","2"))
DefaultAssay(combined_adult_wt_het_p60_john_pom) <- "RNA"
combined_adult_wt_het_p60_john_pom <- NormalizeData(combined_adult_wt_het_p60_john_pom, normalization.method = "LogNormalize", scale.factor = 10000)
combined_adult_wt_het_p60_john_pom <- FindVariableFeatures(combined_adult_wt_het_p60_john_pom, nfeatures = 2000)
combined_adult_wt_het_p60_john_pom <- ScaleData(combined_adult_wt_het_p60_john_pom)
combined_adult_wt_het_p60_john_pom <- RunPCA(combined_adult_wt_het_p60_john_pom, features = VariableFeatures(object = combined_adult_wt_het_p60_john_pom))
combined_adult_wt_het_p60_john_pom <- FindNeighbors(combined_adult_wt_het_p60_john_pom, dims = 1:10)
combined_adult_wt_het_p60_john_pom <- FindClusters(combined_adult_wt_het_p60_john_pom, resolution = 0.3)
set.seed(500)
combined_adult_wt_het_p60_john_pom <- RunUMAP(combined_adult_wt_het_p60_john_pom, dims = 1:10)
nb.cols <- 10
mycolors <- colorRampPalette(brewer.pal(9, "Paired"))(nb.cols)
DimPlot(combined_adult_wt_het_p60_john_pom, reduction = "umap", cols = mycolors, raster = FALSE, label = TRUE)
DimPlot(combined_adult_wt_het_p60_john_pom, reduction = "umap", split.by = "Genotype", cols = mycolors, raster = FALSE, label = TRUE)
DimPlot(combined_adult_wt_het_p60_john_pom, reduction = "umap", split.by = "orig.ident", cols = mycolors, label = TRUE)

Idents(combined_adult_wt_het_p60_john_pom) <- "seurat_clusters"
combined_adult_wt_het_p60_john_pom[["percent.mt"]] <- PercentageFeatureSet(combined_adult_wt_het_p60_john_pom, pattern = "^mt-")
VlnPlot(combined_adult_wt_het_p60_john_pom, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

### Remove doublets and dying cell clusters
Idents(combined_adult_wt_het_p60_john_pom) <- "seurat_clusters"
combined_adult_wt_het_p60_john_pom <- subset(combined_adult_wt_het_p60_john_pom, idents = c("0","1","3","4","5","6","7","8"))

DefaultAssay(combined_adult_wt_het_p60_john_pom) <- "RNA"
combined_adult_wt_het_p60_john_pom <- NormalizeData(combined_adult_wt_het_p60_john_pom, normalization.method = "LogNormalize", scale.factor = 10000)
combined_adult_wt_het_p60_john_pom <- FindVariableFeatures(combined_adult_wt_het_p60_john_pom, nfeatures = 2000)
combined_adult_wt_het_p60_john_pom <- ScaleData(combined_adult_wt_het_p60_john_pom)
combined_adult_wt_het_p60_john_pom <- RunPCA(combined_adult_wt_het_p60_john_pom, features = VariableFeatures(object = combined_adult_wt_het_p60_john_pom))
combined_adult_wt_het_p60_john_pom <- FindNeighbors(combined_adult_wt_het_p60_john_pom, dims = 1:20)
combined_adult_wt_het_p60_john_pom <- FindClusters(combined_adult_wt_het_p60_john_pom, resolution = 0.5)
set.seed(500)
combined_adult_wt_het_p60_john_pom <- RunUMAP(combined_adult_wt_het_p60_john_pom, dims = 1:20)
nb.cols <- 15
mycolors <- colorRampPalette(brewer.pal(9, "Paired"))(nb.cols)
DimPlot(combined_adult_wt_het_p60_john_pom, reduction = "umap", cols = mycolors, raster = FALSE, label = TRUE)

Idents(combined_adult_wt_het_p60_john_pom) <- "seurat_clusters"
combined_adult_wt_het_p60_john_pom[["percent.mt"]] <- PercentageFeatureSet(combined_adult_wt_het_p60_john_pom, pattern = "^mt-")
VlnPlot(combined_adult_wt_het_p60_john_pom, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

Idents(combined_adult_wt_het_p60_john_pom) <- "seurat_clusters"
combined_adult_wt_het_p60_john_pom <- subset(combined_adult_wt_het_p60_john_pom, idents = c("0","1","2","3","4","5","6","7","8","10","12"))

DefaultAssay(combined_adult_wt_het_p60_john_pom) <- "RNA"
combined_adult_wt_het_p60_john_pom <- NormalizeData(combined_adult_wt_het_p60_john_pom, normalization.method = "LogNormalize", scale.factor = 10000)
combined_adult_wt_het_p60_john_pom <- FindVariableFeatures(combined_adult_wt_het_p60_john_pom, nfeatures = 2000)
combined_adult_wt_het_p60_john_pom <- ScaleData(combined_adult_wt_het_p60_john_pom)
combined_adult_wt_het_p60_john_pom <- RunPCA(combined_adult_wt_het_p60_john_pom, features = VariableFeatures(object = combined_adult_wt_het_p60_john_pom))
combined_adult_wt_het_p60_john_pom <- FindNeighbors(combined_adult_wt_het_p60_john_pom, dims = 1:20)
combined_adult_wt_het_p60_john_pom <- FindClusters(combined_adult_wt_het_p60_john_pom, resolution = 0.2)
set.seed(500)
combined_adult_wt_het_p60_john_pom <- RunUMAP(combined_adult_wt_het_p60_john_pom, dims = 1:20)
nb.cols <- 15
mycolors <- colorRampPalette(brewer.pal(9, "Paired"))(nb.cols)
DimPlot(combined_adult_wt_het_p60_john_pom, reduction = "umap", cols = mycolors, raster = FALSE, label = TRUE)
DimPlot(combined_adult_wt_het_p60_john_pom, reduction = "umap", cols = mycolors, split.by = "Genotype", raster = FALSE, label = TRUE)

Idents(combined_adult_wt_het_p60_john_pom) <- "seurat_clusters"
combined_adult_wt_het_p60_john_pom[["percent.mt"]] <- PercentageFeatureSet(combined_adult_wt_het_p60_john_pom, pattern = "^mt-")
VlnPlot(combined_adult_wt_het_p60_john_pom, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

Idents(combined_adult_wt_het_p60_john_pom) <- "seurat_clusters"
combined_adult_wt_het_p60_john_pom <- subset(combined_adult_wt_het_p60_john_pom, idents = c("0","1","2","3","4","5","6","7","9"))

Idents(combined_adult_wt_het_p60_john_pom) <- "seurat_clusters"
new.cluster.ids <- c("TM2","TM1","TM3","IS","Kera","SF","Peri","Schw","CE")
names(new.cluster.ids) <- levels(combined_adult_wt_het_p60_john_pom)
combined_adult_wt_het_p60_john_pom <- RenameIdents(combined_adult_wt_het_p60_john_pom, new.cluster.ids)
combined_adult_wt_het_p60_john_pom$CellType <- Idents(combined_adult_wt_het_p60_john_pom)
combined_adult_wt_het_p60_john_pom$CellType <- factor(combined_adult_wt_het_p60_john_pom$CellType, levels = c('TM1','TM2','TM3','IS','SF','Kera','CE','Peri','Schw'))

combined_adult_wt_het_p60_john_pom$Genotype <- factor(combined_adult_wt_het_p60_john_pom$Genotype, levels = c('WT','het'))
Idents(combined_adult_wt_het_p60_john_pom) <- "Genotype"
new.cluster.ids <- c("ctrl","het")
names(new.cluster.ids) <- levels(combined_adult_wt_het_p60_john_pom)
combined_adult_wt_het_p60_john_pom <- RenameIdents(combined_adult_wt_het_p60_john_pom, new.cluster.ids)
combined_adult_wt_het_p60_john_pom$Genotype <- Idents(combined_adult_wt_het_p60_john_pom)
combined_adult_wt_het_p60_john_pom$Genotype <- factor(combined_adult_wt_het_p60_john_pom$Genotype, levels = c('ctrl','het'))

Idents(combined_adult_wt_het_p60_john_pom) <- "CellType"
mycolors <- c("#a6cee3","#e31a1c","#33a02c","#fb9a99","khaki3","#b2df8a","#cab2d6","#1f78b4","#6a3d9a")
DimPlot(combined_adult_wt_het_p60_john_pom, reduction = "umap", cols = mycolors)
DimPlot(combined_adult_wt_het_p60_john_pom, reduction = "umap", cols = mycolors, split.by = "Genotype")

Idents(combined_adult_wt_het_p60_john_pom) <- "CellType"
gene_list_plot <- c("Lmx1b")
DotPlot(object = combined_adult_wt_het_p60_john_pom, features = gene_list_plot,cols = c("yellow","darkorchid3"))

#Compare WT and het marker genes repectively
Idents(combined_adult_wt_het_p60_john_pom) <- "CellType"
combined_adult_wt_het_p60_john_tm <- subset(combined_adult_wt_het_p60_john_pom, idents = c("TM1","TM2","TM3"))

Idents(combined_adult_wt_het_p60_john_tm) <- "Genotype"
combined_adult_wt_het_p60_john_tm_wt <- subset(combined_adult_wt_het_p60_john_tm, idents = c("ctrl"))
combined_adult_wt_het_p60_john_tm_het <- subset(combined_adult_wt_het_p60_john_tm, idents = c("het"))

Idents(combined_adult_wt_het_p60_john_tm_wt) <- "CellType"
Idents(combined_adult_wt_het_p60_john_tm_het) <- "CellType"
markers_combined_adult_wt_het_p60_john_tm_wt <- FindAllMarkers(combined_adult_wt_het_p60_john_tm_wt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top20 <- markers_combined_adult_wt_het_p60_john_tm_wt %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
DoHeatmap(combined_adult_wt_het_p60_john_tm_wt, features = top20$gene, group.colors = mycolors) +NoLegend() +scale_fill_viridis()
DoHeatmap(combined_adult_wt_het_p60_john_tm_het, features = top20$gene, group.colors = mycolors) +NoLegend() +scale_fill_viridis()
DoHeatmap(combined_adult_wt_het_p60_john_tm, features = top20$gene, group.colors = mycolors) +NoLegend() +scale_fill_viridis()

####### Save all cells combined object
saveRDS(combined_adult_wt_het_p60_john_pom, file = "~/Desktop/Lmx1b/Development/Single_Cell_Sequencing/RDS_files/05302024_combined_v265d_pom.RDS")

####### load all cells combined object
combined_adult_wt_het_p60_john_pom <- readRDS("~/Desktop/Lmx1b/Development/Single_Cell_Sequencing/RDS_files/05302024_combined_v265d_pom.RDS")

######## V265D TM analysis
Idents(combined_adult_wt_het_p60_john_pom) <- "CellType"
combined_adult_wt_het_p60_john_pom$cell.geno <- paste(Idents(combined_adult_wt_het_p60_john_pom), combined_adult_wt_het_p60_john_pom$Genotype, sep = "_")

Idents(combined_adult_wt_het_p60_john_pom) <- "cell.geno"
combined_adult_wt_het_p60_john_pom_tm1_deg <- FindMarkers(combined_adult_wt_het_p60_john_pom, ident.1 = "TM1_ctrl", ident.2 = "TM1_het")
combined_adult_wt_het_p60_john_pom_tm1_deg
write.csv(combined_adult_wt_het_p60_john_pom_tm1_deg, 'combined_adult_wt_het_p60_john_pom_tm1_deg.csv')

combined_adult_wt_het_p60_john_pom_tm2_deg <- FindMarkers(combined_adult_wt_het_p60_john_pom, ident.1 = "TM2_ctrl", ident.2 = "TM2_het")
combined_adult_wt_het_p60_john_pom_tm2_deg
write.csv(combined_adult_wt_het_p60_john_pom_tm2_deg, 'combined_adult_wt_het_p60_john_pom_tm2_deg.csv')

combined_adult_wt_het_p60_john_pom_tm3_deg <- FindMarkers(combined_adult_wt_het_p60_john_pom, ident.1 = "TM3_ctrl", ident.2 = "TM3_het", logfc.threshold = 0.1,)
combined_adult_wt_het_p60_john_pom_tm3_deg
write.csv(combined_adult_wt_het_p60_john_pom_tm3_deg, 'combined_adult_wt_het_p60_john_pom_tm3_deg.csv')

Idents(combined_adult_wt_het_p60_john_pom_tm) <- "Genotype"
combined_adult_wt_het_p60_john_tm_deg <- FindMarkers(combined_adult_wt_het_p60_john_pom_tm, ident.1 = "ctrl", ident.2 = "het")
combined_adult_wt_het_p60_john_tm_deg
write.csv(combined_adult_wt_het_p60_john_tm_deg, 'combined_adult_wt_het_p60_john_tm_deg.csv')

##### Pathway analysis
setwd("~/Desktop")

# All DEG excel files were modified
# This is example for  TM1 cell subtype, repeated for all subtype comparisons
TM1_v265d <- fread("combined_adult_wt_het_p60_john_pom_tm1_deg_de.csv")

Idents(combined_adult_wt_het_p60_john_pom) <- "CellType"
combined_adult_wt_het_p60_john_pom_tm <- subset(combined_adult_wt_het_p60_john_pom, idents = c("TM1","TM2","TM3"))

TM1_v265d_pathway <- bitr(TM1_v265d$SYMBOL[which(TM1_v265d$abs_avg_log2FC > 0)], fromType = "SYMBOL",
                         toType = c("ENSEMBL", "ENTREZID"),
                         OrgDb = org.Mm.eg.db)
universe.df <- bitr(rownames(combined_adult_wt_het_p60_john_pom_tm), fromType = "SYMBOL",
                    toType = c("ENSEMBL", "ENTREZID"),
                    OrgDb = org.Mm.eg.db)
ego <- enrichGO(gene          = TM1_v265d_pathway$ENTREZID,
                universe      = universe.df$ENTREZID,
                OrgDb         = org.Mm.eg.db,
                ont           = "CC", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.1,
                qvalueCutoff  = 0.1)

head(summary(ego))
dotplot(ego, showCategory=30)
TM1_cc <- as.data.frame(summary(ego))
write.csv(TM1_cc, "TM1_cc.csv")

ego <- enrichGO(gene          = TM1_v265d_pathway$ENTREZID,
                universe      = universe.df$ENTREZID,
                OrgDb         = org.Mm.eg.db,
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

head(summary(ego))
dotplot(ego, showCategory=30)
TM1_bp <- as.data.frame(summary(ego))
write.csv(TM1_bp, "TM1_bp.csv")

ego <- enrichGO(gene          = TM1_v265d_pathway$ENTREZID,
                universe      = universe.df$ENTREZID,
                OrgDb         = org.Mm.eg.db,
                ont           = "MF", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

head(summary(ego))
dotplot(ego, showCategory=30)
TM1_mf <- as.data.frame(summary(ego))
write.csv(TM1_mf, "TM1_mf.csv")

##### GSEA analysis
original_gene_list <- TM1_v265d$avg_log2FC
names(original_gene_list) <- TM1_v265d$SYMBOL
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

gse <- gseGO(geneList=gene_list, 
             ont ="CC", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 1, 
             verbose = TRUE, 
             OrgDb = org.Mm.eg.db, 
             pAdjustMethod = "BH")

dotplot(gse, showCategory=20, split=".sign") + facet_grid(.~.sign)

TM1_v265d_gse_cc <- as.data.frame(summary(gse))
write.csv(TM1_v265d_gse_cc, "TM1_v265d_gse_cc.csv")

gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 1, 
             verbose = TRUE, 
             OrgDb = org.Mm.eg.db, 
             pAdjustMethod = "BH")

dotplot(gse, showCategory=20, split=".sign") + facet_grid(.~.sign)

TM1_v265d_gse_bp <- as.data.frame(summary(gse))
write.csv(TM1_v265d_gse_bp, "TM1_v265d_gse_bp.csv")

gse <- gseGO(geneList=gene_list, 
             ont ="MF", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 1, 
             verbose = TRUE, 
             OrgDb = org.Mm.eg.db, 
             pAdjustMethod = "BH")

dotplot(gse, showCategory=20, split=".sign") + facet_grid(.~.sign)

TM1_v265d_gse_mf <- as.data.frame(summary(gse))
write.csv(TM1_v265d_gse_mf, "TM1_v265d_gse_mf.csv")

######### Make barplot of results as on line 482
######### End of code
