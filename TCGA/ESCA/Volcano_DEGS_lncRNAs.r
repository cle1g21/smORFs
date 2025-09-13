library(org.Hs.eg.db)
library(biomaRt)

# read in DESeq2 results for BLCA dataset
res <- read.delim("/lyceum/cle1g21/smORFs/TCGA/ESCA/Outputs/DESeq2results_ESCA_filtered_variants.txt", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
View(res)

# Filter res to include lncRNAs only from bioMart
# Connect to Ensembl using biomaRt
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Retrieve Ensembl IDs and gene symbols for lncRNAs
lncRNA_info <- getBM(
    attributes = 
    c("ensembl_gene_id",  # Ensembl IDs
    "external_gene_name", # Gene symbols
    "gene_biotype"),      # lncRNAs
    filters = "ensembl_gene_id",
    values = rownames(res),
    mart = ensembl
)

# Filter for lncRNAs only
lncRNA_info <- subset(lncRNA_info, gene_biotype == "lncRNA")

# Merge with DESeq2 results
res_lncRNA <- res[rownames(res) %in% lncRNA_info$ensembl_gene_id, ]

# Add gene symbols to lncRNA filtered DESeq2 results
res_lncRNA$gene_symbol <- lncRNA_info$external_gene_name[match(rownames(res_lncRNA), lncRNA_info$ensembl_gene_id)]

res_lncRNA$label <- ifelse(is.na(res_lncRNA$gene_symbol) | res_lncRNA$gene_symbol == "", rownames(res_lncRNA), res_lncRNA$gene_symbol)
View(res_lncRNA)

# remove gene_symbol column as content is now in label column
res_lncRNA$gene_symbol <- NULL

# Save Gene Symbol Mapping for count plots on full_TCGA_BRCA.r
gene_labels <- setNames(lncRNA_info$external_gene_name, lncRNA_info$ensembl_gene_id)
# Save as an RDS file (R's native format for saving objects)
saveRDS(gene_labels, file = "/lyceum/cle1g21/smORFs/TCGA/ESCA/Outputs/gene_labels.rds")

write.table(res_lncRNA, file = "/lyceum/cle1g21/smORFs/TCGA/ESCA/Outputs/DEGS_lncRNAs_ESCA.csv", sep = "\t", row.names = TRUE)

res <- read.delim("/lyceum/cle1g21/smORFs/TCGA/ESCA/Outputs/DEGS_lncRNAs_ESCA.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
View(res)

# Selects top 15 DEGs
top15_upreg <- res |> 
  filter(padj < 0.05, log2FoldChange < 0) |>         # Significant down-regulated genes in BRCA
  arrange(log2FoldChange) |>                        # Arrange by largest logFC
  slice_head(n = 15) |>                    # Top 15 genes
  mutate(Expression = "downregulated") |>                
  rownames_to_column(var = "Gene_ID")           

top15_downreg <- res |> 
  filter(padj < 0.05, log2FoldChange > 0) |>         # Significant up-regulated genes in BRCA
  arrange(desc(log2FoldChange)) |>                  # Arrange by largest logFC
  slice_head(n = 15) |>                    # Top 15 genes
  mutate(Expression = "upregualated") |>               
  rownames_to_column(var = "Gene_ID")

# Combining the top 15 DEGs for labelling
top30 <- bind_rows(top15_upreg, top15_downreg)
View(top30)

## Volcano plot ----
# Visualize volcano plot
volcano <- ggplot() +
  geom_point(
    data = filter(res, padj < 0.05, log2FoldChange < -1),  # Significant down-regulated genes in red
    aes(x = log2FoldChange, y = -log10(padj)),
    color = "blue", alpha = 0.5, size = 1.5 
  ) +
  geom_point(
    data = filter(res, padj < 0.05, log2FoldChange > 1),   # Significant up-regulated genes in green
    aes(x = log2FoldChange, y = -log10(padj)),
    color = "red", alpha = 0.5, size = 1.5 
  ) +
  geom_point(         
    data = filter(res, padj >= 0.05 | (padj < 0.05 & abs(log2FoldChange) <= 1)),  # Non-significant genes in gray
    aes(x = log2FoldChange, y = -log10(padj)),
    color = "grey", alpha= 0.5, size = 1.5                   
  ) +
  geom_text_repel(
    data = top30,                                           # Adds labels for top 30 significant genes
    aes(x = log2FoldChange, y = -log10(padj), label= label), 
        color = "black",
        size = 3, max.overlaps = 15
  ) +
  geom_vline(xintercept = c(-1,1), linetype = "dashed",     # Threshold
             color = "black", linewidth = 0.8) +
  geom_hline(yintercept = -log10(max(res$padj[res$log2FoldChange < 0.05])), 
           linetype = "dashed", color = "black", linewidth = 0.8) +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10(Padj)"
  ) +
  theme_minimal()

print(volcano)

ggsave("/lyceum/cle1g21/smORFs/TCGA/ESCA/Outputs/Figures/volcano_DEGS_only_lncRNA_1.png", plot = volcano, width = 8, height = 6, dpi = 300, bg = "white")
