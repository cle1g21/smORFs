library(GEOquery)
library(DESeq2)
library(dplyr)

setwd("/lyceum/cle1g21/smORFs/")

#read Owens raw counts and sample-info
rna <- read.delim("rna_counts.txt", sep = "\t", header = TRUE)

#read sample-info
coldata <- read.delim("sample_info.txt", sep = "\t")
coldata <-as.data.frame(apply(coldata, 2, as.factor))
head(coldata)

#check that sample names match
all(colnames(rna) == rownames(coldata))
colnames(rna)
rownames(coldata)

#make DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = rna,
                              colData = coldata,
                              design = ~Condition)
dds

#run DESeq
dds <- DESeq(dds)
res <- results(dds)
print(res)


################################################################# try using DeSeq2 with GEO dataset GSE61405
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

################################################################# try using DeSeq2 with GEO dataset GSE79666
#read GEO dataset GSE79666 raw counts and sample-info
rna <- read.delim("GSE79666_raw_counts_GRCh38.p13_NCBI.tsv", sep ="\t", header = TRUE)

coldata <- read.table("GSE79666_info.csv", sep = ",", header = TRUE)
coldata <-as.data.frame(apply(coldata, 2, as.factor))
colnames(coldata)[1] <- "GeneID"

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
dds$Condition <-relevel(dds$Condition, ref = "CTRL")
dds <- DESeq(dds)
res <-results(dds, alpha = 0.05)
res

#check for significant genes
significant_genes <- rownames(res)[which(res$padj < 0.05)]
head(significant_genes)
length(significant_genes)
summary(res)
