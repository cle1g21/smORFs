# Same code as volcano for DEGs but instead used matched nuORFs differentially expressed in BRCA

# Load libraries
library(tidyverse)
library(ggrepel)
library(org.Hs.eg.db)

# read in DESeq2 results for BRCA dataset
res <- read.delim("/lyceum/cle1g21/smORFs/TCGA/BRCA/nuORFs_matched_BRCA.txt", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
View(res)

# Map Ensembl IDs to Gene Symbols
gene_symbols <- mapIds(
    org.Hs.eg.db,           
    keys = row.names(res),  
    column = "SYMBOL",        # Retrieve gene symbols
    keytype = "ENSEMBL",      # Use Ensembl IDs as input
    multiVals = "first"       # Use the first match if multiple exist
)

# Add gene symbols to the DESeq2 results (res)
res$gene_symbol <- gene_symbols

# Replace missing gene symbols with Ensembl IDs
res$label <- ifelse(is.na(res$gene_symbol), rownames(res), res$gene_symbol)

write.table(res, file = "/lyceum/cle1g21/smORFs/TCGA/BRCA/nuORFs_gene_symbols_BRCA.csv", sep = "\t", row.names = TRUE)


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

ggsave("/lyceum/cle1g21/smORFs/TCGA/BRCA/volcano_nuORFs_1.png", plot = volcano, width = 8, height = 6, dpi = 300, bg = "white")

# MYOM1 downregulated and highly significant seen to be commonly dysregulated in various cancers:
# LI, M., SUN, Q. & WANG, X. 2017. Transcriptional landscape of human cancers. Oncotarget, 8, 34534-34551.

