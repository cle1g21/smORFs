library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

pbmc.data <- Read10X(data.dir = "/lyceum/cle1g21/smORFs/TCGA/BRCA/cellxgene/filtered_gene_bc_matrices 3/hg19")
# Read10X takes output from cellranger pipeline from 10X and returns a unique molecular identified (UMI) count matrix
rownames(pbmc.data) # genes/features
colnames(pbmc.data) # cells

# Use count matrix to create Seurat object
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
# project named for record keeping
# min.cells = keep only genes that are expressed in at least x amount of cells. Filters out rare genes and likely noise
# min.features = keep only cells that have at least x amount of genes expressed. Filters out low-quality cells with few detected genes that are dead

pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30] # looks at 30 cells for CD3D, TCL1A, and MS4A1 genes
# . printed where 0s
# Uses sparse.size where only non-zero values and their positions are kept to save storage when >90% of values are 0 in single cell 

# [[ ]] adds columns to the object metadata -> used to keep quality control (QC) stats 
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-") # percentage of mitochondrial genes (genes that start with MT-)

# Show metrics for the first 5 cells
head(pbmc@meta.data, 5)
# orig.ident = the sample or dataset the each cell came
# nCount_RNA = total number of RNA molecules (UMI counts) detected in that cell. High values could be high transcriptional activity or doublets (2 cells captured together)
# nFeature_RNA = how many genes/features in that cell. Low vlaues could indicate low quality cells or dying cells

# Visualise QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")

plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# Healthy cells show a strong positive correlation between the nCount_RNA and nFeature_RNA
# Doublets or poor-quality cells may have high nCount_RNA but low nFeature_RNA

# Quality control filtering to keep only good-quality cells 
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA <2500 & percent.mt <5)
# nFeature_RNA > 200 = keeps cells that express more than 200 genes -> removes empty droplets
# nFeatureRNA <2500 = keeps cells with less than 2500 genes -> removes doublets
# percent.mt < 5 = keeps cells with less than 5% mitochondrial genes -> removes stressed/dying cells

# Normalisation aftersubsetting for high quality genes only
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
# Scaled so each cell has the same total count. Log transform to reduce highly expressed genes and normalize variance
# Can use CLR or RC normalization.method. CLR used for cell surface protein data
# Can use SCTransform instead of LogNormalize -> replaces the need to run NormalizeData, FindVariableFeatures, and ScaleData

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# Identifies the 2000 most variable genes across all cells -> most informative for clustering and dimensionality reduction
# vst best for vUMI-based data. Can also have mean.var.plot (better for SMART-seq when low cell numbers or dense count data -> vst breaks when very small datasets (<100 cells) and dispersion (older method)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10) # VariableFeatures retrieves first 10 genes from FindVariableFeatures (highest biological variability)

# Plot the most variable features
plot1 <- VariableFeaturePlot(pbmc) # scatter plot of all genes in the dataset. The top 2000 variable genes are highlighted
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) # adds labels to plot1 of the 10 most variable genes

# Scaling is a pre-processing step prior to dimensional reduction techniques
# ScaleData() shifts expression of each gene so the mean expression across cells is 0. Scales the expression of each gene so the variance across cells is 1
# Important so highly-expressed genes do not dominate dimensionality reduction
# By default, only variable features are scaled. To scale all features, specific using the features argument
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
# The default regression model is linear, but can be changed via model.use = . Can use poisson or negbinom (raw count data), glmGamPoi (large datasets with many low-count genes (high overdispersion))
# The results are stored in pbmc[["RNA"]]$scale.data

# You can also use ScaleData() to remove unwanted variation for example, mitochondrial contamination via vars.to.regress parameter
# pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
# However, it is recommended to use SCTransform instead if you would like to use vars.to.regress

# Perform PCA on the scaled data
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# features = defaults to the FindVariableFeatures, but can change to specific subset eg. to focus on specific biological pathways (but these must be scaled before using ScaleData())
# Each PC separates major biological states or cells types. For example, positive PC1 is monocyte/macrophage genes and negative PC1 is genes related to T cells. Whereas negative PC2 is B cell genes and positive PC2 is genes related to NK/T cytotoxic cells
# npcs = changes number of PCs -> default is 50 

# You can view which genes drive variation in the PC and identify which PC separates specific cell types by what genes are expressed in that cell type using VizDimLoadings()
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")  + NoLegend()
# Lots of parameters for reduction = . Could use pca, umap, tsne. 
# Can group.by = "celltype" from colnames(pbmc@meta.data) but this column is not automatically included in the meta.data
# Can split.by = to check whether a condition changes cell type coposition or cluster structure

# To assess what number of PCs to use for each dataset, plot DimHeatmap() and ElbowPlot()
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
# Here, you are looking for genes known to be associated with technical variation. For example, mitochondrial genes, robosomal genes, high-abundance housekeeping genes, or genes globally expressed in all cells

# Elbow plot allows easier interpretation of number of PCs to use
# Here, a harmonic mean is calculated to plot a line on the elbow plot. This separates high-variance PCs from low-variance PCs (noise). PCs above the line are more informative and should be used
x <- 1/mean(1/as.numeric(pbmc[["pca"]]@stdev))
# Calculates the harmonic mean of the standard deviations from each PC. Harmonic mean is taken because it is less sensitive to outliers
# The line is plotted on the elbow plot using geom_hline()
ElbowPlot(pbmc, ndims = 50) + geom_hline(yintercept = x, color = "grey")
# ElbowPlot(pbmc) can be used alone and visually cut off where the curve flattens
# It is good practice to repeat downstream analyses with different numbers of PCs -> often the results do not differ dramatically
# Consider if you use 10 PCs you may miss rare cell types and advise to higher amount of PCs in doubt

# Use FindNeighbors() and FindClusters() to identify clusters
# FindNeighbors() uses the euclidean distance (straightest and shortest path between two points) to construct a shared nearest neighbor (SNN) graph based on their PCA-reduced coordinates
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5) # Uses FindNeighbors() SNN graph to identify clusters

# Visualise clusters using UMAP
# cells with very similar gene expression profiles co-localise
pbmc <- RunUMAP(pbmc, dims = 1:10)
# Default metric cosine (scaled data) but can use euclidean (PCA-reduced data), manhattan (when data is sparse)

DimPlot(pbmc, reduction = "umap")

# FindClusters classifies cells that are transcrptionally similar into clusters
table(pbmc$seurat_clusters)

# For this tutorial, cluster 2 is used. You should change this depending what cell type you want to distingish from all others

# Find Markers()finds genes differentially expressed in this cluster (cell type) compared to all other clusters
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
# Can add idnet.2 if you want to compare 1 cluster to another
# Default test.use is wilcox. roc useful for cell markers by generating binary scoring for how well a gene separates groups for marker detection. 0.5 (no separation) → 1 (perfect separation). Use LR to regress out eg. batch effects
# min.pct = minimum percentage of cells in a cluster that express a gene. Default is 0.25
head(cluster2.markers, n = 5)

# FindAllMarkers() finds gene markers for every cluster compared to the rest of the cells, typically these are only positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
    group_by(cluster) %>% # This groups the data by the cluster column. Good to keep for comparable results per cluster not globally if applying multiple filters
    dplyr::filter(avg_log2FC > 1) # filters rows to only keep those where avg_log2FC is greater than 1 (gene is expressed at least 2x more in the cluster compared to others)

# VlnPlot() shows expression probability distributions of 1 gene across clusters (cell types)
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

# FeaturePlot() visualises feature expression on a PCA plot
FeaturePlot(pbmc, c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
    "CD8A"))

# In the same way as VlnPlot and FeaturePlot, DotPlot() visualises how feature expression changes across different clusters

# Use canonical markers upregulated in each cluster to match the unbiased clustering to known cell types
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
    "NK", "DC", "Platelet") # order for each cluster

# names() sets the names of each cluster ID (level) to those in the above vector
names(new.cluster.ids) <- levels(pbmc)

# RenameIdents() replaces the current default cluster labels (e.g., "0", "1", "2", etc.) with custom names
pbmc <- RenameIdents(pbmc, new.cluster.ids)
# Idents(pbmc) controls how cells are grouped in visualisations like DimPlot() and VlnPlot()

# Visualise clusterrs using UMAP
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# Saving complete plot
plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + xlab("UMAP 1") + ylab("UMAP 2") +
    theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + 
    guides(colour = guide_legend(override.aes = list(size = 10)))

ggsave(filename = "/lyceum/cle1g21/smORFs/TCGA/BRCA/Outputs/pbmc3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)
