# padj table of DEGs across all cancers filtered by at least 1 gene < 0.05 padj
    # went from 54562 obs. for all DEGs to 41861 obs. that met the < 0.05 padj filter

# filter significant DEGs to include only those which are nuORFs in lncRNAs 
    # went from 41861 obs. to 535 obs. Shows there are 535 nuORFs in lncRNAs significantly dysregulated across the cancers
    # wrote table to sig_nuORFlncRNAs_padj_across_all_cancers.csv

---------------------------------------------------------------------------------------------------------------------------------------------

# Log2 fold change table across all cancers filtered by at least 1 gene < 0.05 padj
    # went from 54562 obs. to 41861 obs. Same as padj table confirms it has been done correctly

# filter significant DEGs to include only those which are nuORFs in lncRNAs 
    # went from 41861 obs. to 535 obs. Shows the log2 fold change of the 535 nuORFs in lncRNAs significantly dysregulated across the cancers
    # wrote table to sig_nuORFlncRNAs_log2foldchange_across_all_cancers.csv

---------------------------------------------------------------------------------------------------------------------------------------------

# Binary scoring table where:
    # 1 = a nuORFs in lncRNAs meets the threshold of padj <0.05 and 10 fold change (log2(10) = 3.32) in at least 6 cancers showing the same direction of fold cahnge in those cancers
    # 0 = does not meet the threshold

# Visualisation using heatmap function
  # first heatmap of colour coded statistic table
  # second heatmap shows gradient of up and down regulation

# Visualisation of the expression(log fold change) of each gene with 4fold change in 6 or more cancers (bar graphs)
# Plot counts for each of the 9 genes in the heatmap across each cancer type in a grid

# Load libraries
library(dplyr)
library(purrr)
library(tibble)
library(stringr)
library(pheatmap)
library(org.Hs.eg.db)
library(ggplot2)
library(DESeq2)
library(dplyr)
# Get paths to all matched files under /TCGA/**/Outputs/ named DESeq2results_[cancer]_filtered_variants.txt
base_path <- "/lyceum/cle1g21/smORFs/TCGA"
file_paths <- list.files(
  path = base_path,
  pattern = "^DESeq2results_.*_filtered_variants\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)

nuORFs_lncRNAs <- read.delim("/lyceum/cle1g21/smORFs/TCGA/BRCA/Outputs/lncRNAs_nuORFs_matched_BRCA.txt", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

# padj value ---------------------------------------------------------------------------------------------------------------------------------------------
# Read and pull out padj value
read_log2fc <- function(filepath) {
  # Extract cancer type from filename (e.g., "DESeq2results_BRCA_filtered_variants.txt")
  cancer <- str_extract(basename(filepath), "(?<=DESeq2results_)[A-Z]+(?=_filtered_variants)")
  # Select and rename column
  df <- read.delim(filepath, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  # Select and rename column
  df <- df %>%
    dplyr::select(padj) %>%
    rownames_to_column(var = "Gene_ID") %>%
    dplyr::rename(!!cancer := padj)
  return(df)
}

# Read all files and merge into one wide table
padjvalue_table <- map(file_paths, read_log2fc) %>%
  purrr::reduce(full_join, by = "Gene_ID") %>%
  column_to_rownames(var = "Gene_ID")

head(padjvalue_table)

# checks LINC01614 is in the table
# "ENSG00000230838" %in% rownames(padjvalue_table)  # Should return TRUE
# prints row of padj across cancers for LINC01614
# padjvalue_table["ENSG00000230838", ]  

# filter by at least 1 gene < 0.05 padj across cancers
# Count how many padj values are < 0.05 per gene
significant_genes <- rowSums(padjvalue_table < 0.05, na.rm = TRUE) > 0

# Keep only DEGs with at least one significant padj across the cancers
sigDEGs_padjvalue_table <- padjvalue_table[significant_genes, ]
View(sigDEGs_padjvalue_table)

str(padjvalue_table) # 54562 obs.
str(sigDEGs_padjvalue_table) # 41861 obs.

# Of these significant DEGs in the padj table how many are nuORFs in lncRNAs?
signuORFs_lncRNAs_padjvalue_table <- sigDEGs_padjvalue_table[rownames(sigDEGs_padjvalue_table) %in% rownames(nuORFs_lncRNAs), ]
n_sig_nuORFs <- nrow(signuORFs_lncRNAs_padjvalue_table)
str(n_sig_nuORFs)
str(signuORFs_lncRNAs_padjvalue_table) # 535 obs.

write.csv(signuORFs_lncRNAs_padjvalue_table, "/lyceum/cle1g21/smORFs/TCGA/Outputs/Tables/sig_nuORFlncRNAs_padj_across_all_cancers.csv")

# log2 fold change ---------------------------------------------------------------------------------------------------------------------------------------------
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
  purrr::reduce(full_join, by = "Gene_ID") %>%
  column_to_rownames(var = "Gene_ID")

# Keep only DEGs with at least one significant padj (<0.05) across the cancers
sigDEGs_foldchange_table <- foldchange_table[rownames(foldchange_table) %in% rownames(sigDEGs_padjvalue_table), ]
View(sigfoldchange_table)

str(foldchange_table) # 54562 obs.
str(sigfoldchange_table) # 41861 obs.

# Checks filtering has been done correctly. Should return TRUE
# all(rownames(sigfoldchange_table) == rownames(sigDEGs_padjvalue_table))

# Of these significant DEGs in the log2 fold change table how many are nuORFs in lncRNAs?
signuORFs_lncRNAs_log2foldchange_table <- sigDEGs_foldchange_table[rownames(sigDEGs_foldchange_table) %in% rownames(nuORFs_lncRNAs), ]

str(signuORFs_lncRNAs_log2foldchange_table) # 535 obs.

# Checks filtering for nuORFs has been done correctly. Should return TRUE
# all(rownames(signuORFs_lncRNAs_padjvalue_table) == rownames(signuORFs_lncRNAs_log2foldchange_table))

write.csv(signuORFs_lncRNAs_log2foldchange_table, "/lyceum/cle1g21/smORFs/TCGA/Outputs/Tables/sig_nuORFlncRNAs_log2foldchange_across_all_cancers.csv")

# binary scoring table ---------------------------------------------------------------------------------------------------------------------------------------------
signuORFs_lncRNAs_log2foldchange_table <- read.csv("/lyceum/cle1g21/smORFs/TCGA/Outputs/Tables/sig_nuORFlncRNAs_log2foldchange_across_all_cancers.csv", row.names = 1, stringsAsFactors = FALSE)
signuORFs_lncRNAs_padjvalue_table <- read.csv("/lyceum/cle1g21/smORFs/TCGA/Outputs/Tables/sig_nuORFlncRNAs_padj_across_all_cancers.csv", row.names = 1, stringsAsFactors = FALSE)

# Set threshold
fc_threshold <- log2(4) # = 1 magnitude of change

# TRUE = passes the fold change threshold
threshold_foldchange <- abs(signuORFs_lncRNAs_log2foldchange_table) >= fc_threshold

# Combine to get "1/0" matrix where 1 if both conditions are TRUE
score_matrix <- (signuORFs_lncRNAs_padjvalue_table < 0.05) & threshold_foldchange

View(score_matrix)

# Is each nuORF lncRNA expressed in 6 or more cancers? eg. does it have a padj and log2 fold change values in at least 6 cancers
# for each gene it counts how many cancer types there is a number and keeps those are are expressed in 6 or more
expressed_in_6_or_more <- rowSums(!is.na(score_matrix)) >= 6

# keep only significant (padj <0.05) nuORFs in lncNRAs in at least 6 cancers showing the same direction of fold cahnge in those cancers
# Get the direction of change: +1 for upregulated, -1 for downregulated
direction_matrix <- sign(signuORFs_lncRNAs_log2foldchange_table)
 
View(direction_matrix)

# Only keep direction where both padj and FC thresholds are met
direction_matrix[!score_matrix] <- 0  # Set everything else to 0

# Count number of cancers with consistent upregulation or downregulation
up_count <- rowSums(direction_matrix == 1)
down_count <- rowSums(direction_matrix == -1)

# Is each nuORF lncRNA consistently up/down in 6 or more cancers?
#consistent_filter <- (up_count >= 6) | (down_count >= 6)

# Filter the table to get only consistently changing nuORFs AND expressed in at least 6 cancers 
consistent_filter <- ((up_count >= 6) | (down_count >= 6)) & expressed_in_6_or_more

# makes sure consistent_filter has gene names to match rows of direction_matrix
names(consistent_filter) <- rownames(direction_matrix)

# keep only consistantly changing nuORFs in consistent_nuORFs
consistent_nuORFs <- direction_matrix[which(consistent_filter), , drop = FALSE]

# sum(is.na(rownames(signuORFs_lncRNAs_padjvalue_table)))  # should be 0
# sum(is.na(rownames(signuORFs_lncRNAs_log2foldchange_table)))  # should be 0

print(consistent_nuORFs)

write.csv(consistent_nuORFs, "/lyceum/cle1g21/smORFs/TCGA/Outputs/Tables/sig_nuORFlncRNAs_log_6pluscancers.csv")

# visualise table ---------------------------------------------------------------------------------------------------------------------------------------------
# gradient to visualise strength of fold change 
# makes the fold change table only include genes from consitent_nuORFs (those who meet the threshold)
heatmap_matrix <- signuORFs_lncRNAs_log2foldchange_table[rownames(consistent_nuORFs), ]

# Grey out log2FC values that do not pass fold change threshold
heatmap_matrix[!score_matrix[rownames(consistent_nuORFs), ]] <- NA

# Map Ensembl IDs to Gene Symbols
gene_symbols <- mapIds(
    org.Hs.eg.db,           
    keys = row.names(heatmap_matrix),  
    column = "SYMBOL",        # Retrieve gene symbols
    keytype = "ENSEMBL",      # Use Ensembl IDs as input
    multiVals = "first"       # Use the first match if multiple exist
)

display_labels <- ifelse(is.na(gene_symbols), rownames(heatmap_matrix), gene_symbols)

color_palette <- colorRampPalette(c("blue", "white", "red"))(100) # shows gradient of up vs down regulated

png("/lyceum/cle1g21/smORFs/TCGA/Outputs/Figures/heatmap/heatmap_nuORFlncRNAs_4fold_6pluscancers_gradient_3.png", width = 12, height = 14, units = "in", res = 300)

# Check if we have enough genes for clustering
if(nrow(heatmap_matrix) >= 2) {
  cluster_rows_setting <- TRUE
} else {
  cluster_rows_setting <- FALSE
}

heatmap_matrix <- as.matrix(heatmap_matrix)

p2 <- pheatmap(heatmap_matrix,
         cluster_rows = cluster_rows_setting,
         cluster_cols = FALSE,
         color = color_palette,
         breaks = seq(-5, 5, length.out = 101),
         #main = "Consistently Regulated nuORFs in lncRNAs 4 Fold Change Across >= 6 Cancers",
         labels_row = display_labels,
         show_rownames = TRUE,
         legend_labels = c("Down", "Neutral", "Up"),
         fontsize_col = 16,        # FONT SIZE: Column labels (cancer types) = 16
         fontsize_row = 16,        # FONT SIZE: Row labels (gene names) = 14
         fontsize = 16,            # FONT SIZE: General text = 14
         angle_col = 45,
         na_col = "grey90")        # Grey out NAs of genes that do not meet the threshold

dev.off()

# signuORFs_lncRNAs_log2foldchange_table["ENSG00000258884", ] # to see specific fold change counts of a gene across the cancer types

# extract the row order from raw fold change clustering to keep gene order consistent
row_order <- rownames(heatmap_matrix)[p2$tree_row$order]

# block colours 1 (red) for upregulated, -1 (blue) downregulated, 0 (white) if did not meet threshold
ordered_consistent_nuORFs <- consistent_nuORFs[row_order, ]

color_palette <- colorRampPalette(c("blue", "white", "red"))(3)

ordered_display_labels <- display_labels[row_order]

png("/lyceum/cle1g21/smORFs/TCGA/Outputs/Figures/heatmap/heatmap_nuORFlncRNAs_4fold_6pluscancers_binary_2.png", width = 10, height = 12, units = "in", res = 300)

p3 <- pheatmap(ordered_consistent_nuORFs,
         cluster_rows = FALSE, # keep order from previous heatmap that uses true fold change for clustering
         cluster_cols = FALSE,
         color = color_palette,
         main = "Consistently Regulated nuORFs in lncRNAs 4 Fold Change Across >= 6 Cancers",
         labels_row = ordered_display_labels,
         show_rownames = TRUE,
         legend_breaks = c(-1, 0, 1),
         legend_labels = c("Down", "Neutral", "Up"),
         fontsize_col = 14,
         na_col = "grey90")

dev.off()

-------------------------------------------------------------------------------------------------------------------------------------------------------------

# Visualisation of log fold change of each gene with 4fold change in 6 or more cancers

logfc_table <- 
logfc_table <- read.csv("/lyceum/cle1g21/smORFs/TCGA/Outputs/Tables/sig_nuORFlncRNAs_log2foldchange_across_all_cancers.csv", row.names = 1, stringsAsFactors = FALSE)
consistent_nuORFs <- read.csv("/lyceum/cle1g21/smORFs/TCGA/Outputs/Tables/sig_nuORFlncRNAs_4fold_6pluscancers.csv", row.names = 1, stringsAsFactors = FALSE)
padj_table <- read.csv("/lyceum/cle1g21/smORFs/TCGA/Outputs/Tables/sig_nuORFlncRNAs_padj_across_all_cancers.csv", row.names = 1, stringsAsFactors = FALSE)
logfc_table <- read.csv("/lyceum/cle1g21/smORFs/TCGA/Outputs/Tables/sig_nuORFlncRNAs_log2foldchange_across_all_cancers.csv", row.names = 1, stringsAsFactors = FALSE)

# Convert all columns to numeric
padj_table[] <- lapply(padj_table, as.numeric)

# Get gene symbols
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = rownames(logfc_table),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Loop through each gene
for (gene_id in rownames(consistent_nuORFs)) {

  # Get data
  gene_fc <- logfc_table[gene_id, , drop = FALSE]
  gene_padj <- padj_table[gene_id, , drop = FALSE]
  if (nrow(gene_fc) == 0) next

  gene_fc_df <- data.frame(
    Cancer = colnames(gene_fc),
    log2FC = as.numeric(gene_fc[1, ]),
    padj = as.numeric(gene_padj[1, ])
  )

  # Add significance stars
  gene_fc_df$significance <- cut(
    gene_fc_df$padj,
    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
    labels = c("***", "**", "*", "")
  )

  gene_fc_df$direction <- ifelse(gene_fc_df$log2FC > 0, "Up", "Down")
  

  # Title text
  symbol <- gene_symbols[gene_id]
  title_text <- ifelse(is.na(symbol), gene_id, paste0(symbol, " (", gene_id, ")"))

  # Output path
  out_path <- paste0("/lyceum/cle1g21/smORFs/TCGA/Outputs/Figures/Bar_charts/log2foldchange_", gene_id, ".png")

  # Save plot
  png(out_path, width = 8, height = 6, units = "in", res = 300)

  print(
    ggplot(gene_fc_df, aes(x = Cancer, y = log2FC, fill = direction)) +
  geom_col() +
  geom_text(
  aes(label = significance),
  vjust = ifelse(gene_fc_df$log2FC >= 0, -0.5, 1.5)
) +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue"), guide = "none") +
      theme_minimal() +
      labs(
        title = paste("log2 Fold Change of", title_text, "Across Cancers"),
        y = "log2 Fold Change",
        x = "Cancer Type"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )

  dev.off()
}
ENSG00000176728

# Plot counts for each of the 9 genes across each cancer type
# Load 4 fold significant genes dysregulated in 6 or more cancers 
genes <- rownames(read.csv(
  "/lyceum/cle1g21/smORFs/TCGA/Outputs/Tables/sig_nuORFlncRNAs_4fold_6pluscancers.csv",
  row.names = 1
))

# Target gene
gene <- "ENSG00000174365"

# Get gene symbols once for all genes
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = genes,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# cancer types and paths to dds objects
cancer_types <- c("BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", "HNSC", "KICH")

# where are the DESeq2 dds objects are stored?
dds_path_base <- "/lyceum/cle1g21/smORFs/TCGA/{cancer}/Outputs/dds.rds"

# where to put the output?
output_dir <- "/lyceum/cle1g21/smORFs/TCGA/Outputs/Figures/Plot_Counts"

# Loop over cancers and genes to make plot counts
for (cancer in cancer_types) {
  # Load DESeq2 object
  dds_file <- gsub("{cancer}", cancer, dds_path_base, fixed = TRUE)
  dds <- readRDS(dds_file)

  # Clean gene IDs in dds to match your CSV
  rownames(dds) <- sub("\\..*", "", rownames(dds))

  for (gene in genes) {
    # Get symbol and format title
    symbol <- gene_symbols[gene]
    title_text <- ifelse(is.na(symbol), paste(gene, "-", cancer), paste0(symbol, " (", gene, ") - ", cancer))

    # Output file path
    outfile <- file.path(output_dir, paste0("Plotcount_", gene, "_", cancer, ".png"))

    # Save plot
    png(outfile, width = 800, height = 600)
    try({
      plotCounts(dds, gene = gene, intgroup = "Condition", main = title_text)
    }, silent = TRUE)
    dev.off()
  }
}


# For specific gene from MASCOT
# Specify gene of interest
gene_id <- "ENSG00000174365"

# Get gene data
gene_fc <- logfc_table[gene_id, , drop = FALSE]
gene_padj <- padj_table[gene_id, , drop = FALSE]

gene_fc_df <- data.frame(
  Cancer = colnames(gene_fc),
  log2FC = as.numeric(gene_fc[1, ]),
  padj = as.numeric(gene_padj[1, ])
)

# Add significance stars
gene_fc_df$significance <- cut(
  gene_fc_df$padj,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "")
)

gene_fc_df$direction <- ifelse(gene_fc_df$log2FC > 0, "Up", "Down")

# Get gene symbol
gene_symbol <- mapIds(
  org.Hs.eg.db,
  keys = gene_id,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

title_text <- ifelse(is.na(gene_symbol), gene_id, paste0(gene_symbol, " (", gene_id, ")"))

outfile <- file.path("/lyceum/cle1g21/smORFs/TCGA/ESCA/Outputs/Figures", paste0("log2_fold_change_MASCOT_", gene_id, ".png"))
  
png(outfile, width = 8, height = 6, units = "in", res = 300)
  
  par(mfrow = c(3, 3), mar = c(5, 4, 4, 2) + 0.1)  # Adjust margins

# Plot
ggplot(gene_fc_df, aes(x = Cancer, y = log2FC, fill = direction)) +
  geom_col() +
  geom_text(aes(label = significance), vjust = ifelse(gene_fc_df$log2FC >= 0, -0.5, 1.5)) +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue"), guide = "none") +
  theme_minimal() +
  labs(
    title = paste("log2 Fold Change of", title_text, "Across Cancers"),
    y = "log2 Fold Change",
    x = "Cancer Type"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

# Bar graph suggested by Owen. number of differentially expressed nuORFs (y-axis) that are consistently dysregulated (padj < 0.05) across 1 or more, 2 or more, 3 or more, 5 or more, 6 or more, 7 or more, 8 or more cancers (x-axis)
padj_table <- read.csv("/lyceum/cle1g21/smORFs/TCGA/Outputs/Tables/sig_nuORFlncRNAs_padj_across_all_cancers.csv", row.names = 1, stringsAsFactors = FALSE)
logfc_table <- read.csv("/lyceum/cle1g21/smORFs/TCGA/Outputs/Tables/sig_nuORFlncRNAs_log2foldchange_across_all_cancers.csv", row.names = 1, stringsAsFactors = FALSE)

# Define threshold
padj_threshold <- 0.05

# Calculate number of cancers with padj < 0.05 for each gene
sig_counts <- rowSums(padj_table < padj_threshold, na.rm = TRUE)
View(sig_counts)

# Define the thresholds you want to plot
thresholds <- c(1, 2, 3, 5, 6, 7, 8, 9)

# For each threshold, count how many genes meet or exceed it
bar_data <- data.frame(
  Cancers = paste0("\u2265", thresholds),  # \u2265 is the Unicode for "≥",
  Num_nuORFs = sapply(thresholds, function(t) sum(sig_counts >= t))
)

# Set output file path
outfile <- "/lyceum/cle1g21/smORFs/TCGA/Outputs/Figures/Bar_charts/nuORFs_barplot_2.png"

# Open PNG device
png(outfile, width = 8, height = 6, units = "in", res = 300)

# Plot using ggplot2
bar <- ggplot(bar_data, aes(x = Cancers, y = Num_nuORFs)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = Num_nuORFs), vjust = -0.5, size = 5) +  # <-- This adds the labels to the top fo each bar
  labs(
    title = " ",
    x = "Number of Cancers (padj < 0.05)",
    y = "Number of nuORFs"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16))

# Print the plot to the device
print(bar)

# Close device
dev.off()

# Find genes significant in all 9 cancers
nuORFs_in_all_9 <- names(sig_counts)[sig_counts == 9]
# Print the gene IDs
print(nuORFs_in_all_9)

nuORFs_in_all_9 <- names(sig_counts)[sig_counts == 9]
# Print the gene IDs
print(nuORFs_in_all_9)


# ggplot 3x3 grid plot counts for each gene of interest for appendices
gene <- "ENSG00000174365" # change for each gene
cancer_types <- c("BLCA","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KICH")
dds_path_base <- "/lyceum/cle1g21/smORFs/TCGA/{cancer}/Outputs/dds.rds"
output_dir <- "/lyceum/cle1g21/smORFs/TCGA/Outputs/Figures/Plot_Counts"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

gene_symbol <- mapIds(org.Hs.eg.db, keys = gene, column = "SYMBOL",
                      keytype = "ENSEMBL", multiVals = "first")

# extract counts per cancer
get_counts <- function(ca) {
  dds_file <- gsub("{cancer}", ca, dds_path_base, fixed = TRUE)
  dds <- readRDS(dds_file)
  rownames(dds) <- sub("\\..*", "", rownames(dds))  # drop Ensembl version suffixes
  if (!(gene %in% rownames(dds))) return(NULL)


  data <- plotCounts(dds, 
                   gene = gene, 
                   intgroup = "Condition", 
                   returnData = TRUE)

  # map Condition labels to N/T if possible (keeps original if unknown)
  levs <- levels(data$Condition)
  labs <- ifelse(grepl("normal|solid tissue normal|ctrl", tolower(levs)), "N",
          ifelse(grepl("tumou?r|primary", tolower(levs)), "T", levs))
  data$Condition <- factor(data$Condition, levels = levs, labels = labs)
  data$Condition <- factor(data$Condition, levels = c("N","T"))  # keep N left if both present

  data$count_adj <- pmax(data$count, 0.5)  # avoid log10 issues at 0
  data$Cancer <- ca
  data
}

# build one data frame for all cancers
df_all <- do.call(rbind, lapply(cancer_types, get_counts))

outfile <- file.path(output_dir, paste0("Plotcount_", gene, "_ALL.png"))

png(outfile, width = 10, height = 14, units = "in", res = 300) 

ggplot(df_all, aes(x=Condition, y=count_adj)) + 
  geom_jitter(width = 0.15, shape = 21, color = "black", fill = "black", alpha = 0.3) + 
  scale_y_log10() +
  facet_wrap(~ Cancer, ncol = 3) +
  labs(
    x = NULL,
    y = "Normalized counts",
    title = paste0(ifelse(is.na(gene_symbol), gene, paste0(gene_symbol, " (", gene, ")")))
  ) +
  theme_bw() + 
  theme(axis.text.x  = element_text(size = 14),              # x tick labels
    axis.text.y  = element_text(size = 14),       
    strip.text = element_text(size = 16),
    axis.title.y = element_text(size = 16, margin = margin(r = 8)))

dev.off()

# plot count for individual cancer types
# pick one cancer
one_cancer <- "ESCA"
gene <- "ENSG00000174365" # change for each gene

# subset (and drop unused factor levels)
df1 <- droplevels(subset(df_all, Cancer == one_cancer))

png("/lyceum/cle1g21/smORFs/TCGA/Outputs/Figures/Plot_Counts/test.png", width = 1800, height = 1800, res = 300)

ggplot(df1, aes(x = Condition, y = count_adj)) +
  geom_jitter(width = 0.15, shape = 21, color = "black", fill = "black", alpha = 0.3) +
  scale_y_log10() +
  labs(
    x = NULL,
    y = "Normalized counts",
    title = paste0(ifelse(is.na(gene_symbol), gene, paste0(gene_symbol, " (", gene, ")")),
                   " — ", one_cancer)
  ) +
  theme_bw() +
  theme(
    axis.text.x  = element_text(size = 14),              # x tick labels
    axis.text.y  = element_text(size = 14),       
    strip.text = element_text(size = 16),
    axis.title.y = element_text(size = 16, margin = margin(r = 8)))

dev.off()
