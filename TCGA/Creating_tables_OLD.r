# makes log2 fold change table of genes across all cancers - write TCGA/combined_log2FC_by_cancer_DEGsbeforefilters.csv
  # filters unfiltered log2fold change table to only include those that are nuORFs in lncRNAs goes from 54562 obs to 535 obs
# made 90% and 100% up/downregulated genes -> wrote into .csv files and used for visualisation:
  # of this those are are nuORF in lncRNAs -> obs goes to 28 - write TCGA/nuORFlncRNAs_significant_percent100_upregulatedDEGs_across_cancers.csv
# visualisation of LINC01614 (upregulated significantly in all the cancers) using bar plot - write TCGA/Outputs/log2foldchange_ENSG00000230838.png
# *tried heatmap but not very informative too many genes -> try again with nuORFs in lncRNAs
# makes padj table of genes across all cancers - write TCGA/combined_padj_by_cancer_DEGsbeforefilters.csv
# filtering the padj table by:
  # 1. highest significant genes at the top
  # 2. removing any genes less than 0.05 across all cancers - write TCGA/combined_padj_filteredbysignificance.csv
  # 3. use the table ordered by significance and genes less than 0.05 to see which of these are nuORFs in lncRNAs - write TCGA/nuORFs_lncRNAs_combined_padj_filteredbysignificance.csv

# Load libraries
library(dplyr)
library(purrr)
library(tibble)
library(stringr)
library(pheatmap)

# Log2 fold change table-----------------------------------------------------------------------------------------------------------------------------
# Get paths to all matched files under /TCGA/**/Outputs/ named DESeq2results_[cancer]_filtered_variants.txt
base_path <- "/lyceum/cle1g21/smORFs/TCGA"
file_paths <- list.files(
  path = base_path,
  pattern = "^DESeq2results_.*_filtered_variants\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)

# Read and pull out log2FoldChange
read_log2fc <- function(filepath) {
  # Extract cancer type from filename (e.g., "DESeq2results_BRCA_filtered_variants.txt")
  cancer <- str_extract(basename(filepath), "(?<=DESeq2results_)[A-Z]+(?=_filtered_variants)")
  # Select and rename column
  df <- read.delim(filepath, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  # Select and rename column
  df <- df %>%
    select(log2FoldChange) %>%
    rownames_to_column(var = "Gene_ID") %>%
    rename(!!cancer := log2FoldChange)
  return(df)
}

# Read all files and merge into one wide table
foldchange_table <- map(file_paths, read_log2fc) %>%
  reduce(full_join, by = "Gene_ID") %>%
  column_to_rownames(var = "Gene_ID")

View(foldchange_table)

write.csv(foldchange_table, "/lyceum/cle1g21/smORFs/TCGA/combined_log2FC_by_cancer_DEGsbeforefilters.csv")

# filter this list to only include nuORFs in lncRNAs
nuORFs_lncRNAs <- read.delim("/lyceum/cle1g21/smORFs/TCGA/BRCA/Outputs/lncRNAs_nuORFs_matched_BRCA.txt", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
View(nuORFs_lncRNAs)

# of foldchange_table unfiltered which are nuORF lncRNAs
nuORFlncRNAs_foldchange_table <- foldchange_table[rownames(foldchange_table) %in% rownames(nuORFs_lncRNAs), ]
View(nuORFlncRNAs_foldchange_table)

write.csv(nuORFlncRNAs_padjvalue_table, "/lyceum/cle1g21/smORFs/TCGA/nuORFs_lncRNAs_combined_padj_filteredbysignificance.csv")

# Read in and view the logfoldchange table
foldchange_table <- read.csv("/lyceum/cle1g21/smORFs/TCGA/combined_log2FC_by_cancer_DEGsbeforefilters.csv", row.names = 1, stringsAsFactors = FALSE)
View(foldchange_table)

# Which genes change in the same direction across the cancers?----------------------------------------------------------------------------
# create table counting the number of upregulation (> 0) or downregulation (< 0) across the cancers for each gene
direction_summary <- data.frame(
  up_count = rowSums(foldchange_table > 0, na.rm = TRUE),
  down_count = rowSums(foldchange_table < 0, na.rm = TRUE),
  total_nonNA = rowSums(!is.na(foldchange_table))
)
View(direction_summary)

# Genes that are upregulated in 100 of the cancers
percent100_up <- direction_summary %>%
  filter(up_count / total_nonNA >= 1.0, # e.g. 100% 
  total_nonNA >= 6)  # seen in at least 6 cancers

# Add log2FC values by rownames
percent100_up_with_fc <- cbind(
  percent100_up,
  foldchange_table[rownames(percent100_up), ]
)

View(percent100_up_with_fc)
write.csv(percent100_up_with_fc, "/lyceum/cle1g21/smORFs/TCGA/percent100_upregulatedDEGs_across_cancers.csv")

# of these which are significant 0.05
significant_percent100_up <- percent100_up_with_fc[rownames(percent100_up_with_fc) %in% rownames(significant_ordered_DEGs), ]
str(significant_percent100_up)

write.csv(significant_percent100_up, "/lyceum/cle1g21/smORFs/TCGA/significant_percent100_upregulatedDEGs_across_cancers.csv")
significant_percent100_up["ENSG00000230838", ]

# of these significant which are nuORFs in lncRNAs
nuORFlncRNAs_sig100upreg <- significant_percent100_up[rownames(significant_percent100_up) %in% rownames(nuORFs_lncRNAs), ]
str(nuORFlncRNAs_sig100upreg)
write.csv(nuORFlncRNAs_sig100upreg, "/lyceum/cle1g21/smORFs/TCGA/nuORFlncRNAs_significant_percent100_upregulatedDEGs_across_cancers.csv")

# Genes that are downregulated in 100% of the cancers
percent100_down <- direction_summary %>%
  filter(down_count / total_nonNA >= 1.0, # e.g. 90% or more
  total_nonNA >= 6)  # seen in at least 6 cancers

# Add log2FC values by rownames
percent100_down_with_fc <- cbind(
  percent100_down,
  foldchange_table[rownames(percent100_down), ]
)
View(percent100_down_with_fc)
write.csv(percent100_down_with_fc, "/lyceum/cle1g21/smORFs/TCGA/percent90_downregulatedDEGs_across_cancers.csv")

# of these which are significant 0.05
significant_percent100_down <- percent100_down_with_fc[rownames(percent100_down_with_fc) %in% rownames(significant_ordered_DEGs), ]

write.csv(significant_percent100_down, "/lyceum/cle1g21/smORFs/TCGA/significant_percent100_downregulatedDEGs_across_cancers.csv")

# of these significant which are nuORFs in lncRNAs
nuORFlncRNAs_sig100downreg <- significant_percent100_down[rownames(significant_percent100_down) %in% rownames(nuORFs_lncRNAs), ]
str(nuORFlncRNAs_sig100downreg)
write.csv(nuORFlncRNAs_sig100downreg, "/lyceum/cle1g21/smORFs/TCGA/nuORFlncRNAs_significant_percent100_downregulatedDEGs_across_cancers.csv")

# Visualisations----------------------------------------------------------------------------------------------------------------------------------
# Visualisation of specific significant gene dysregulated in 100% of the cancers
# For LINC01614:
# Grab log2FC values (assumes gene IDs are rownames)
gene_fc <- significant_percent100_up["ENSG00000230838", ]

# Drop any metadata columns (like up_count, down_count, total_nonNA)
gene_fc <- gene_fc[, !(colnames(gene_fc) %in% c("up_count", "down_count", "total_nonNA"))]

# make dataframe
gene_fc_df <- data.frame(
  Cancer = colnames(gene_fc),
  log2FC = as.numeric(gene_fc)
)

png("/lyceum/cle1g21/smORFs/TCGA/Outputs/log2foldchange_ENSG00000230838.png", width = 8, height = 6, units = "in", res = 300)

ggplot(gene_fc_df, aes(x = Cancer, y = log2FC, fill = log2FC > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("red"), guide = "none") +
  theme_minimal() +
  labs(
    title = "log2 Fold Change of ENSG00000230838 Across Cancers",
    y = "log2 Fold Change",
    x = "Cancer Type"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

# heatmap of log2 fold changes across all cancers
# Combine genes of significant up and down genes consistantly directional across 100% of the cancers

# Remove count columns to keep only cancer fold change columns
to_combineup <- nuORFlncRNAs_sig100upreg[, c("BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", "HNSC", "KICH")]

# Remove count columns to keep only cancer fold change columns
to_combinedown <- nuORFlncRNAs_sig100downreg[, c("BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", "HNSC", "KICH")]

# Combine row-wise
heatmap_matrix <- rbind(to_combineup, to_combinedown) # rbind as cancer columns same but rownames different

png("/lyceum/cle1g21/smORFs/TCGA/Outputs/log2FC_heatmap_sig100percent_upanddown_nuORFslncRNAs.png", width = 10, height = 12, units = "in", res = 300)

# Plot
pheatmap(heatmap_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "none",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Log2 Fold Change of Significant nuORFs in lncRNAs Across Cancers",
         fontsize_row = 6)

dev.off()

# P value table----------------------------------------------------------------------------------------------------------------------------------------

# Read and pull out padj value
read_log2fc <- function(filepath) {
  # Extract cancer type from filename (e.g., "DESeq2results_BRCA_filtered_variants.txt")
  cancer <- str_extract(basename(filepath), "(?<=DESeq2results_)[A-Z]+(?=_filtered_variants)")
  # Select and rename column
  df <- read.delim(filepath, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  # Select and rename column
  df <- df %>%
    select(padj) %>%
    rownames_to_column(var = "Gene_ID") %>%
    rename(!!cancer := padj)
  return(df)
}

# Read all files and merge into one wide table
padjvalue_table <- map(file_paths, read_log2fc) %>%
  reduce(full_join, by = "Gene_ID") %>%
  column_to_rownames(var = "Gene_ID")

write.csv(padjvalue_table, "/lyceum/cle1g21/smORFs/TCGA/combined_padj_by_cancer_DEGsbeforefilters.csv")

# Read in and view the logfoldchange table
padjvalue_table <- read.csv("/lyceum/cle1g21/smORFs/TCGA/combined_padj_by_cancer_DEGsbeforefilters.csv", stringsAsFactors = FALSE)
View(padjvalue_table)

# makes genes in 'X' column as rownames
rownames(padjvalue_table) <- padjvalue_table$X
# remove 'X' column
padjvalue_table$X <- NULL
View(padjvalue_table)

# checks LINC01614 is in the table
"ENSG00000230838" %in% rownames(padjvalue_table)  # Should return TRUE
# prints row of padj across cancers for LINC01614
padjvalue_table["ENSG00000230838", ]  

# Filtering the padj table by most significant genes --------------------------------------------------------------------------------------------------
# adds new column for min padj for that gene across the cancers
padjvalue_table$min_padj <- apply(padjvalue_table, 1, min, na.rm = TRUE) 

# reorders the genes by smallest padj -> most significant DEGs at the top 
padjvalue_table_ordered <- padjvalue_table[order(padjvalue_table$min_padj), ] 
View(padjvalue_table_ordered)

# remove any genes where all of the cancers are non-significant
# exclude 'min_padj'
padj_columns <- padjvalue_table_ordered[, !colnames(padjvalue_table_ordered) %in% "min_padj"] # selects only padj columns and ignores min padj

# Count how many padj values are < 0.05 per gene
significant_genes <- rowSums(padj_columns < 0.05, na.rm = TRUE) > 0

# Keep only genes with at least one significant padj across the cancers
significant_ordered_DEGs <- padjvalue_table_ordered[significant_genes, ]
View(significant_ordered_DEGs)

significant_ordered_DEGs["ENSG00000230838", ]

# filter this list to only include nuORFs in lncRNAs
res <- read.delim("/lyceum/cle1g21/smORFs/TCGA/BRCA/Outputs/lncRNAs_nuORFs_matched_BRCA.txt", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

write.csv(significant_ordered_DEGs, "/lyceum/cle1g21/smORFs/TCGA/combined_padj_filteredbysignificance.csv")

# filter this list to only include nuORFs in lncRNAs
nuORFs_lncRNAs <- read.delim("/lyceum/cle1g21/smORFs/TCGA/BRCA/Outputs/lncRNAs_nuORFs_matched_BRCA.txt", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
View(nuORFs_lncRNAs)

padjvalue_table <- read.csv("/lyceum/cle1g21/smORFs/TCGA/combined_padj_filteredbysignificance.csv", stringsAsFactors = FALSE)
View(padjvalue_table)

# make ensmbl ids row names
# Set the row names to Gene name
rownames(padjvalue_table) <- padjvalue_table[, "X"]

# Sort gene name in gene_ids object
gene_ids <- padjvalue_table[, "X"]

# Remove gene name column, row names are gene names
padjvalue_table <- padjvalue_table[, -1]

# of padjvalue_table which are nuORF lncRNAs
nuORFlncRNAs_padjvalue_table <- padjvalue_table[rownames(padjvalue_table) %in% rownames(nuORFs_lncRNAs), ]
View(nuORFlncRNAs_padjvalue_table)

write.csv(nuORFlncRNAs_padjvalue_table, "/lyceum/cle1g21/smORFs/TCGA/nuORFs_lncRNAs_combined_padj_filteredbysignificance.csv")
