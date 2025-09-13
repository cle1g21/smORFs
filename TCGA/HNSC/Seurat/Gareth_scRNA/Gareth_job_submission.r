# load necessary libraries
cat("Loading libraries...\n"); flush.console()
library(Seurat)
library(ggplot2)
library(dplyr)
library(SingleR)
library(celldex)
library(SummarizedExperiment)

# # Read in Seurat object
cat("Creating seurat object...\n"); flush.console()
counts <- readRDS("counts_matrix.rds")
metadata <- readRDS("metadata_df.rds")

data <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata,
  min.features = 200,
  min.cells = 30
)
grep("LINC01614", rownames(data), value = TRUE)
# cat("Assays present:\n"); flush.console()
# print(Assays(data))

# pct.mito outlier detection
cat("Calculating mitochondrial threshold...\n"); flush.console()
mad(data@meta.data$pct.mito)
thresh <- median(data@meta.data$pct.mito) + 3 * mad(data@meta.data$pct.mito) # stats method for outliers in non-normally distributed in ScRNA-seq data

# # Initial QC metrics
# cat("Generating QC plots...\n"); flush.console()
# plot1 <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "pct.mito"), slot = "counts")
# plot2 <- FeatureScatter(data, feature1 ="nCount_RNA", feature2 = "pct.mito") 
# plot3 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# apply subsetting
cat("Filtering cells based on nFeature_RNA and pct.mito...\n"); flush.console()
data <- subset(data, subset = nFeature_RNA <= 6000 & pct.mito < thresh) # remove doublets and stressed/dying cells

# split by patient ID before normalisation to avoid batch effects
# SCTransform works well by normalising each datset separately and creating separate RPCAs and aligns the mutual nearest neighbours (anchors)
cat("Splitting object by donor_id...\n"); flush.console()
HNSC <- SplitObject(data, split.by = "donor_id")
cat("Split into", length(HNSC), "objects\n"); flush.console()

# run SCTransform individually on each sample using lapply() -> treated as separate datasets during integration to correct for batch effects, technical variation, or inter-sample differences
cat("Running SCTransform on each split object...\n"); flush.console()
HNSC <- lapply(X = HNSC, FUN = function(x) {
  cat("Processing", unique(x$donor_id), "...\n"); flush.console()
  SCTransform(x, min_cells = 30) # attempt to keep LINC01614 appearing in 41 cells
}) # no need to regress out pct.mito as it is already removed in the subsetted data

# Save intermediate result
cat("Saving SCTransformed objects...\n"); flush.console()
saveRDS(HNSC, file = "SCTransform_2.rds")

# HNSC <- readRDS("/lyceum/cle1g21/smORFs/SCTransform.rds")
print(grep("LINC01614", rownames(HNSC[[1]]), value = TRUE))

# Select integration features
cat("Selecting integration features...\n"); flush.console()
features <- SelectIntegrationFeatures(object.list = HNSC, nfeatures = 3000) # Selects the top 3000 highly variable genes that are shared across datasets -> used to align the datasets during integration

# Prep for SCT integration -> Required before running FindIntegrationAnchors()
cat("Preparing for integration...\n"); flush.console()
HNSC <- PrepSCTIntegration(object.list = HNSC, anchor.features = features)

# PCA before integration -> needed for anchor-finding using rpca
cat("Running PCA on each object before integration...\n"); flush.console()
HNSC <- lapply(X = HNSC, FUN = function(x) {
  cat("Running PCA for", unique(x$donor_id), "...\n"); flush.console()
  RunPCA(x, features = features)
})

# Find anchors -> rpca requires PCA to be run on each object beforehand
cat("Finding integration anchors...\n"); flush.console()
HNSCC.anchors <- FindIntegrationAnchors(object.list = HNSC, normalization.method = "SCT",
                                        anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 5)

# Integrate data
cat("Integrating data...\n"); flush.console()
HNSCC.integrated <- IntegrateData(anchorset = HNSCC.anchors, normalization.method = "SCT", dims = 1:30)

# Downstream analysis -> PCA for UMAP
cat("Running PCA on integrated object...\n"); flush.console()
HNSCC.integrated <- RunPCA(HNSCC.integrated, verbose = FALSE)
cat("Running UMAP...\n"); flush.console()
HNSCC.integrated <- RunUMAP(HNSCC.integrated, reduction = "pca", dims = 1:20)
cat("Finding neighbors...\n"); flush.console()
HNSCC.integrated <- FindNeighbors(object = HNSCC.integrated, dims = 1:20)
cat("Clustering cells...\n"); flush.console()
HNSCC.integrated <- FindClusters(object = HNSCC.integrated)

# Find markers
cat("Finding marker genes...\n"); flush.console()
Markers <- FindAllMarkers(HNSCC.integrated, min.pct = 0.25, logfc.threshold = 0.5, only.pos = T)
# FindMarkers object
saveRDS(Markers, "/lyceum/cle1g21/smORFs/TCGA/findMarkers.rds")

# Save final object
cat("Saving integrated Seurat object...\n"); flush.console()
saveRDS(HNSCC.integrated, file = "/lyceum/cle1g21/smORFs/TCGA/complete_seurat_2.rds")

HNSCC.integrated <- readRDS("/lyceum/cle1g21/smORFs/TCGA/complete_seurat_2.rds")
# grep("LINC01614", rownames(HNSCC.integrated), value = TRUE)

cat("Generating plots...\n"); flush.console()
# (Then paste the ggsave() lines)
# ggsave("/lyceum/cle1g21/smORFs/TCGA/vlnplot_qc.png", plot = plot1, width = 8, height = 4)
# ggsave("/lyceum/cle1g21/smORFs/TCGA/scatter_nCount_vs_pct.mito.png", plot = plot2, width = 6, height = 6)
# ggsave("/lyceum/cle1g21/smORFs/TCGA/scatter_nCount_vs_nFeature.png", plot = plot3, width = 6, height = 6)

umap <- DimPlot(HNSCC.integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave("/lyceum/cle1g21/smORFs/TCGA/umap_2.png", plot = umap, width = 8, height = 6)

cat("Plotting expression for LINC01614 and PRRX1...\n"); flush.console()

vlnPlot_LINC01614 <- VlnPlot(HNSCC.integrated, features = c("LINC01614", "PRRX1-ENSG00000116132"))
ggsave(filename = "/lyceum/cle1g21/smORFs/TCGA/gareth_vlnPlot_LINC01614_2.jpg", height = 7, width = 12, plot = vlnPlot_LINC01614, dpi = 300)

DefaultAssay(HNSCC.integrated) <- "RNA"
featurePlot_LINC01614 <- FeaturePlot(HNSCC.integrated, features = c("LINC01614"))
ggsave(filename = "/lyceum/cle1g21/smORFs/TCGA/gareth_featurePlot_LINC01614.jpg", height = 6, width = 8, plot = featurePlot_LINC01614, dpi = 300)
featurePlot_LINC02154 <- FeaturePlot(HNSCC.integrated, features = c("ENSG00000244675.2"))
ggsave(filename = "/lyceum/cle1g21/smORFs/TCGA/gareth_featurePlot_LINC02154.jpg", height = 6, width = 8, plot = featurePlot_LINC01614, dpi = 300)
grep("LINC02321", rownames(HNSCC.integrated[["RNA"]]), value = TRUE)


dotPlot_LINC01614 <- DotPlot(HNSCC.integrated, features = "LINC01614")
ggsave(filename = "/lyceum/cle1g21/smORFs/TCGA/gareth_dotPlot_LINC01614_2.jpg", height = 7, width = 12, plot = dotPlot_LINC01614, dpi = 300)

cat("All steps completed successfully.\n"); flush.console()

# ADDED experimenting for cluster labelling
HNSCC.integrated <- readRDS("/lyceum/cle1g21/smORFs/TCGA/BRCA/Seurat/Gareth_scRNA/complete_seurat_2.rds")
Markers <- readRDS("/lyceum/cle1g21/smORFs/TCGA/BRCA/Seurat/Gareth_scRNA/findMarkers.rds")

# labels cell types from .h5ad file in metadata$cell_type
metadata <- readRDS("/lyceum/cle1g21/smORFs/TCGA/metadata_df.rds")

# name clusters according to metadata$cell_type
Idents(HNSCC.integrated) <- HNSCC.integrated@meta.data$cell_type

dimplot <- DimPlot(HNSCC.integrated, reduction = "umap", label = FALSE, group.by = "cell_type")
ggsave("/lyceum/cle1g21/smORFs/TCGA/BRCA/Seurat/Gareth_scRNA/UMAP/umap_labelled_3.png", plot = dimplot, width = 12, height = 6)

# Load necessary library
library(patchwork)

# Combine side by side
combined_plot <- dimplot + featurePlot_LINC01614

# Save to a PNG file
ggsave("/lyceum/cle1g21/smORFs/TCGA/combined_plot.png", plot = combined_plot, width = 10, height = 20, dpi = 300)
