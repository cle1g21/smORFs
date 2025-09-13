## convert h5ad to seurat rds

# load necessary libraries
library(Seurat)
library(reticulate)
library(ggplot2)
library(dplyr)

data <- readRDS("seurat_object.rds")
head(data@meta.data, 5)

# pct.mito outlier detection
mad(data@meta.data$pct.mito)
thresh <- median(data@meta.data$pct.mito) + 3 * mad(data@meta.data$pct.mito) # stats method for outliers in non-normally distributed in ScRNA-seq data

# Initial QC metrics
plot1 <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "pct.mito"))

# apply subsetting
data <- subset(data, subset = nFeature_RNA <= 6000 & pct.mito < thresh) # remove doublets and stressed/dying cells

# split by patient ID before normalisation to avoid batch effects
# SCTransform works well by normalising each datset separately and creating separate RPCAs and aligns the mutual nearest neighbours (anchors)
HNSC <- SplitObject(data, split.by = "donor_id")
colnames(data@meta.data)

# run SCTransform individually on each sample using lapply()
HNSC <- lapply(X = HNSC, FUN = SCTransform) # no need to regress out pct.mito as it is already removed in the subsetted data

# Running SCTransform on assay: RNAtted data
# Running SCTransform on layer: counts
# vst.flavor='v2' set. Using model with fixed slope and excluding poisson genes.
# Variance stabilizing transformation of count matrix of size 21597 by 9612
# Model formula is y ~ log_umi
# Get Negative Binomial regression parameters per gene
# Using 2000 genes, 5000 cells
# Error in getGlobalsAndPackages(expr, envir = envir, globals = globals) : 
#   The total size of the 19 globals exported for future expression ('FUN()') is 505.63 MiB.. This exceeds the maximum allowed size of 500.00 MiB (option 'future.globals.maxSize'). The three largest globals are 'FUN' (484.88 MiB of class 'function'), 'umi_bin' (19.21 MiB of class 'numeric') and 'data_step1' (1.46 MiB of class 'list')

# Initial QC metrics
plot1 <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "pct.mito"))
plot2 <- FeatureScatter(data, feature1 ="nCount_RNA", feature2 = "pct.mito") 
plot3 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# filter removing empty droplets, doublets, and stressed/dying cells
data <- subset(data, subset = nFeature_RNA > 200 & # keeps cells that express more than 200 genes
  nFeature_RNA < 2500 & # keeps cells with less than 6000 genes
  pct.mito < 20) # keeps cells with less than 20% mitochondrial genes

# Replot QC metrics - is the data tighter and cleaner?
plot1_after <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "pct.mito"))
plot2_after <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "pct.mito")
plot3_after <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# SCTransform replaces NormalizeData, ScaleData, FIndVariableFeatures
data <- SCTransform(data, vars.to.regress = "pct.mito")

saveRDS(data, file = "SCTransform.rds")

data <- readRDS("SCTransform.rds")

# PCA
data <- RunPCA(data, features = VariableFeatures(object = data))

# which genes drive variation between the PCs?
VizDimLoadings(data, dims = 1:2, reduction = "pca")
DimPlot(data, reduction = "pca") + NoLegend()

# picking PCs to use
DimHeatmap(data, dims = 1:15, cells = 500, balanced = TRUE)

x <- 1/mean(1/as.numeric(data[["pca"]]@stdev))
ElbowPlot(data, ndims = 50) + geom_hline(yintercept = x, color = "grey")
# 26 PCs chosen

# find clusters
data <- FindNeighbors(data, dims = 1:26)
data <- FindClusters(data, resolution = 0.5)

# visualise clusters via UMAP
data <- RunUMAP(data, dims = 1:26)
DimPlot(data, reduction = "umap")

table(data$seurat_clusters)

# gene marker for fibroblasts?
fibro_marker <- FindMarkers(data, ident.1 = 2)
head(fibro_marker, n = 5)

# gene marker for every cluster
all_markers <- FindAllMarkers(data, only.pos = TRUE)

all_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.5, p_val_adj < 0.05)

# expression of a gene across cell types
vlnPlot_LINC01614 <- VlnPlot(data, features = c("LINC01614", "PRRX1-ENSG00000116132"))
ggsave(filename = "/lyceum/cle1g21/smORFs/TCGA/HNSC/Outputs/Figures/gareth_vlnPlot_LINC01614.jpg", height = 7, width = 12, plot = vlnPlot_LINC01614, dpi = 300)

# expression of a gene across clusters
featurePlot_LINC01614 <- FeaturePlot(data, features = c("LINC01614", "PRRX1-ENSG00000116132"))
ggsave(filename = "/lyceum/cle1g21/smORFs/TCGA/HNSC/Outputs/Figures/gareth_featurePlot_LINC01614.jpg", height = 7, width = 12, plot = featurePlot_LINC01614, dpi = 300)

# expression of a gene across clusters
dotPlot_LINC01614 <- DotPlot(data, features = c("LINC01614", "PRRX1-ENSG00000116132"))
ggsave(filename = "/lyceum/cle1g21/smORFs/TCGA/HNSC/Outputs/Figures/gareth_dotPlot_LINC01614.jpg", height = 7, width = 12, plot = dotPlot_LINC01614, dpi = 300)

umap <- DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(filename = "/lyceum/cle1g21/smORFs/TCGA/HNSC/Outputs/Figures/gareth_umap.jpg", height = 7, width = 12, plot = umap, dpi = 300)


