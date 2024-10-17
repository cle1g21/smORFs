library(GEOquery)
library(DESeq2)

#reads raw counts for GSE273464 
count_table <- read.csv("GSE273464_raw_counts_GRCh38.p13_NCBI.tsv", header = TRUE, sep = "\t")
head(count_table)
#changes the raw counts to  data.frame
as.data.frame(count_table)
dim(count_table)

#reads design description for GSE273464
design_table <- read.csv("design_tags.csv", header = TRUE)
#changes row names in design_table to match column names in counts_table
colnames(design_table) <- c("GeneID", "condition")
colnames(design_table)

#removed gene_id column to stop error in DESeq matrix. put gene_ids in a new file
gene_ids <- count_table$GeneID
count_table_clean <- count_table[, -1]

#check column for counts matches row for data_info
ncol(count_table)
nrow(design_table)

#construct DeSeq dataset
dds <- DESeqDataSetFromMatrix(countData = count_table_clean,
                              colData = design_table,
                              design = ~condition)

dds

#setting factor level
dds$condition <- relevel(dds$condition, ref = "PBS")

dds$condition

dds <- DESeq(dds)
res <- results(dds)
