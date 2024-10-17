library(GEOquery)
library(DESeq2)

count_table <- read.csv("GSE273464_raw_counts_GRCh38.p13_NCBI.tsv", header = TRUE, sep = "\t")
head(count_table)
as.data.frame(count_table)
dim(count_table)

design_table <- read.csv("design_tags.csv", header = TRUE)
colnames(design_table)

dds <- DESeqDataSetFromMatrix(countData = count_table,
                              colData = design_table,
                              design = ~rhPOSTN)

colnames(design_table) <- c("GeneID", "condition")
