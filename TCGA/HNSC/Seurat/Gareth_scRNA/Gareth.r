## convert h5ad to seurat rds

# load necessary libraries
library(Seurat)
library(reticulate)
library(ggplot2)
library(dplyr)

# Set the Python environment for reticulate
#use_python("/local/software/conda/miniconda-py3-new/bin/python", required = TRUE)
use_condaenv("anndata_env", required = TRUE)
# Check the Python config
py_config()

ad <- import("anndata")

# Load the h5ad file
data <- ad$read_h5ad("/lyceum/cle1g21/smORFs/TCGA/BRCA/cellxgene/f72bc93c-b13e-4413-8fb9-0cd73f456b0d.h5ad")

# Extract metadata
# Convert obs metadata to an R data.frame
metadata <- py_to_r(data$obs)
# Extract barcodes
cell_ids <- rownames(py_to_r(data$obs))

# Extract AnnData count matrix as a R matrix and transpose
counts <- t(py_to_r(data$raw$X))

# Extract gene info (gene IDs and gene symbols) and assign to counts matrix
gene_meta <- py_to_r(data$raw$var)
# Need to transpose because in AnnData rows = cells and columns = genes, whereas in R you want rownames = genes and colnames = cells
gene_ids <- as.character(py_to_r(data$raw$var_names$to_list()))
gene_symbols <- as.character(gene_meta$feature_name)

# Replace missing or empty symbols with Ensembl ID
gene_symbols_clean <- ifelse(
  is.na(gene_symbols) | gene_symbols == "",
  gene_ids,
  gene_symbols
)

# Add suffix only to duplicates (not to all)
gene_symbols_unique <- make.unique(gene_symbols_clean)

rownames(counts) <- gene_symbols_unique # assigns gene symbols as rownames of counts matrix

colnames(counts) <- cell_ids # assigns cell barcodes as the column names of counts matrix so they match the rows (cells) in metadata

all(colnames(counts) == rownames(metadata))  # Should return TRUE. Seurat expects the metadata rownames to match the column names of the counts matrix

# Save them as RDS files
saveRDS(counts, file = "counts_matrix.rds")
saveRDS(metadata, file = "metadata_df.rds")

counts <- readRDS("counts_matrix.rds")
metadata <- readRDS("metadata_df.rds")
sum(counts["LINC01614", ] > 0)
colnames(metadata)
metadata$cell_type

seu <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata,
  min.features = 200,
  min.cells = 30
)

saveRDS(seu, file = "seurat_object.rds")

data <- readRDS("seurat_object.rds")

# store gene annotation (ensembl id and gene symbol)
gene_annotation <- data.frame(
  gene_name = rownames(data),
  ensembl_id = gene_ids[match(rownames(data), gene_names)],
  gene_symbol = gene_symbols[match(rownames(data), gene_names)],
  stringsAsFactors = FALSE
)

saveRDS(gene_annotation, "gene_annotation.rds")

head(data@meta.data, 5)

# Initial QC metrics
plot1 <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "pct.mito"))
plot2 <- FeatureScatter(data, feature1 ="nCount_RNA", feature2 = "pct.mito") 
plot3 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# filter removing empty droplets, doublets, and stressed/dying cells
data <- subset(data, subset = nFeature_RNA > 200 & # keeps cells that express more than 200 genes
  nFeature_RNA < 6000 & # keeps cells with less than 6000 genes
  pct.mito < 20) # keeps cells with less than 20% mitochondrial genes

# Replot QC metrics - is the data tighter and cleaner?
plot1_after <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "pct.mito"))
plot2_after <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "pct.mito")
plot3_after <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Normalise the subsetted data
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)

data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(data), 10) # finds top 10 most variable to label in plot

plot4 <- VariableFeaturePlot(data)
LabelPoints(plot = plot4, points = top10, repel = TRUE)

# scaling
all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes, vars.to.regress = c("nCount_RNA", "pct.mito"))

saveRDS(data, file = "current.rds")

data <- readRDS("current.rds")

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


