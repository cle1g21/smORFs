# do any covariates (stage, subtype, or grade) predict or correlate with gene expression of significant genes across 100% of cancers?
# created loop for interesting covariates for each of the 9 genes that are 4fold and >60% of cancers
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tools)

# ENSG00000230838 x stage
meta <- read.csv("/mainfs/ddnb/TCGA_data/full_TCGA-BRCA_meta.csv")
cts <- read.csv("/mainfs/ddnb/TCGA_data/full_TCGA-BRCA_matrix.csv")

View(meta)
colnames(cts)
rownames(cts)

# Setup cts ----
# Set the row names to Gene name
rownames(cts) <- cts[, "X"]

# Sort gene name in gene_ids object
gene_ids <- cts[, "X"]
# Remove gene name column, row names are gene names
cts_clean <- cts[, -1]

# makes column names in cts the same as X and barcode in meta
colnames(cts_clean) <- gsub("\\.", "-", colnames(cts_clean)) 

# keep only samples (rows) in mtadata whos sample id (X and barcode columns) is in column names of cts
meta <- meta[meta$X %in% colnames(cts_clean), ]
# keep only samples found in meta in the cts. DESeq2 assumes the same order of the cts and meta
cts_clean <- cts_clean[, meta$X]
# check they are aligned correctly before DESeq2
all(colnames(cts_clean) == meta$X)


# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = cts_clean,
                              colData = meta,
                              design = ~ 1)


# Transform counts for visualisation
# vst normalises and removes the mean–variance relationship - otherwsie the variance will depend on the mean. highly expressed genes dominate your analysis, and low-expression genes are lost in the noise
vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)  # log2-transformed counts

head(vsd_mat)
nrow(vsd_mat)

# genes_of_interest <- c(
#   "ENSG00000229970", "ENSG00000230838", "ENSG00000231246",
#   "ENSG00000231768", "ENSG00000235385", "ENSG00000241684",
#   "ENSG00000244675", "ENSG00000253177", "ENSG00000258884"
# )

genes_of_interest <- "ENSG00000214293"

# 🎯 List of selected metadata columns to plot
selected_meta <- c(
  "paper_age_at_initial_pathologic_diagnosis",
  "ajcc_pathologic_t",
  "paper_pathologic_stage",
  "definition",
  "gender",
  "paper_BRCA_Subtype_PAM50",
  "paper_vital_status",
  "primary_diagnosis",
  "race",
  "tumor_descriptor"
)

# Loop over each gene
for (gene in genes_of_interest) {
  
  row_match <- grep(paste0("^", gene, "\\."), rownames(vsd_mat), value = TRUE)
  if (length(row_match) == 0) {
    warning(paste("Gene", gene, "not found"))
    next
  }

  gene_expr <- data.frame(
    sample_id = colnames(vsd_mat),
    expression = vsd_mat[row_match, ]
  )

  meta_expr <- merge(gene_expr, meta, by.x = "sample_id", by.y = "X")

  output_dir <- file.path("/lyceum/cle1g21/smORFs/TCGA/Outputs/Figures/Metadata", gene)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Loop over selected metadata columns
  for (meta_col in selected_meta) {
    
    # Skip if metadata column doesn't exist
    if (!(meta_col %in% colnames(meta_expr))) {
      warning(paste("Skipping", meta_col, "- not found"))
      next
    }

    output_file <- file.path(output_dir, paste0(gene, "_metadata_", meta_col, ".png"))
    is_numeric <- is.numeric(meta_expr[[meta_col]])

    png(output_file, width = 16, height = 12, units = "in", res = 300)

    if (is_numeric) {
      # Scatter for numeric
      p <- ggplot(meta_expr, aes_string(x = meta_col, y = "expression")) +
        geom_point(alpha = 0.3, size = 1.2) +
        geom_smooth(method = "lm", se = TRUE, color = "blue", linetype = "dashed", linewidth = 1) +
        theme_minimal(base_size = 14) +
        labs(title = paste("Expression of", gene, "by", meta_col),
             y = "VST-normalized Expression", x = meta_col)
    } else {
      # Boxplot for categorical
      meta_expr[[meta_col]] <- factor(meta_expr[[meta_col]])  # ensure factor
      p <- ggplot(meta_expr, aes_string(x = meta_col, y = "expression", fill = meta_col)) +
        geom_boxplot() +
        theme_minimal(base_size = 14) +
        labs(title = paste("Expression of", gene, "by", meta_col),
             y = "VST-normalized Expression", x = meta_col) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12)) +
        guides(fill = "none")
    }

    print(p)
    dev.off()
  }  # end metadata loop
}  # end gene loop

# catagoring stages into:
#Stage 0, Stage I, Stage IA, Stage IB	Localised
#Stage II, Stage IIA, Stage IIB, Stage III, Stage IIIA, Stage IIIB, Stage IIIC	Lymph Node Involvement
#Stage IV	Metastatic
#Stage X, NA	Unknown

#reclassify stages
meta_expr$stage_group <- dplyr::case_when(
  meta_expr$ajcc_pathologic_stage %in% c("Stage 0", "Stage I", "Stage IA", "Stage IB") ~ "Localised",
  meta_expr$ajcc_pathologic_stage %in% c("Stage II", "Stage IIA", "Stage IIB", 
                                         "Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC") ~ "Lymph Node Involvement",
  meta_expr$ajcc_pathologic_stage %in% c("Stage IV") ~ "Metastasis",
  TRUE ~ "Unknown"  # This covers NA or Stage X
)

#plot
png("/lyceum/cle1g21/smORFs/TCGA/Outputs/ENSG00000029559.7_metadata_stage_classifications.png", width = 16, height = 12, units = "in", res = 300)

# Plot expression across stage
ggplot(meta_expr, aes(x = stage_group, y = expression, fill = stage_group)) + 
  geom_boxplot() +
  theme_minimal(base_size = 14) +
  labs(title = "Expression of ENSG00000029559.7 by Reclassified Stages",
       y = "VST-normalized Expression", x = "Stage") +
 theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    #legend.position = "bottom",         # moves the legend to the bottom
    #legend.box = "horizontal"           # lays it out horizontally
  ) +
  guides(fill = "none")             # removes legend

dev.off()