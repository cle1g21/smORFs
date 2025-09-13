library(tibble)
library(ggplot2)
library(ggrepel)

res <- read.delim("/lyceum/cle1g21/smORFs/TCGA/BLCA/Outputs/DEGS_lncRNAs_BLCA.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
View(res)

# Match the DESeq2 results to the nuORFs dataset
# Filters res by padj, keeping the rows from the res data frame where the padj value is less than 0.05
resOrdered <- res[which(res$padj < 0.05), ]
head(resOrdered)

# read in nuORFs dataset from https://smorfs.ddnetbio.com/
ncORFs <- read.csv("/lyceum/cle1g21/smORFs/nuORF_atlas.csv")
head(ncORFs)

# filter ncORFs to remove rows where NAs are found in rownames (gene name) of resOrdered and saved to a dataset called noNA_resOrdered. need to do this because the next line gives an error if NAs are found
noNA_resOrdered <- resOrdered[!is.na(rownames(resOrdered)), ]

# subset the resOrdered data frame to include only the rows where the Ensembl column is in ncORFs$Gene_id. so the resOrdered only include genes in the nuORFs list
resOrdered_nuORFs <- noNA_resOrdered[rownames(noNA_resOrdered) %in% ncORFs$Gene_id, ]
View(resOrdered_nuORFs)
str(resOrdered_nuORFs)

# save file of matching nuORFs
write.table(resOrdered_nuORFs, file = "/lyceum/cle1g21/smORFs/TCGA/BLCA/Outputs/lncRNAs_nuORFs_matched_BLCA.txt", sep = "\t", row.names = TRUE)

res <- read.delim("/lyceum/cle1g21/smORFs/TCGA/BLCA/Outputs/lncRNAs_nuORFs_matched_BLCA.txt", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
View(res)
nrow(res)
# Visualisation of lncRNA nuORFs in volcano plot
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
  mutate(Expression = "upregulated") |>               
  rownames_to_column(var = "Gene_ID")

write.csv(top15_downreg, "/lyceum/cle1g21/smORFs/TCGA/BLCA/Outputs/top15upreg_lncRNA_matched_nuORFs.csv")

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

ggsave("/lyceum/cle1g21/smORFs/TCGA/BLCA/Outputs/Figures/volcano_lncRNAs_nuORFs_2.png", plot = volcano, width = 8, height = 6, dpi = 300, bg = "white")
