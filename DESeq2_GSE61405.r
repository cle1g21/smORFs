#** dataset GSE61405: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61405 **
#** dataset did not work. no significant genes found. need to check publication methodology (significance value, year of the paper, analysis process etc) and perform PCA before DESeq2 analysis **

library(GEOquery)
library(DESeq2)
library(dplyr)

#read GEO dataset GSE61405 raw counts and sample-info
rna <- read.delim("GSE61405_raw_counts.tsv", sep ="\t", header = TRUE)
rna
head(rna)

coldata <- read.table("sample_info.csv", sep = ",", header = TRUE)
coldata <-as.data.frame(apply(coldata, 2, as.factor))
head(coldata)
rownames(coldata)

#check that sample names match
colnames(coldata)[1] <- "GeneID"
colnames(coldata)
head(coldata)
colnames(rna)

#remove GeneID column
gene_ids <- rna$GeneID
rna_clean <- rna[, -1]
colnames(rna_clean)

#make DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = rna_clean,
                              colData = coldata,
                              design = ~Condition)
dds

#make sure comparison is against control
dds$Condition <-relevel(dds$Condition, ref = "control")

#run DESeq
dds <- DESeq(dds)
res <-results(dds, alpha = 0.05)
res

#check for significant genes
significant_genes <- rownames(res)[which(res$padj < 0.05)]
head(significant_genes)
length(significant_genes)
summary(res)
