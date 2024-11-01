#** Datasets from https://github.com/SGDDNB/translational_regulation/tree/master/sample_data **
#** Replicate Owens findings using DESeq2 and match to nuORFs (filter by nc rna) from https://smorfs.ddnetbio.com/ **

library(GEOquery)
library(DESeq2)
library(dplyr)

setwd("/lyceum/cle1g21/smORFs/deltaTE_Owen")

#read Owens raw counts and sample-info
rna <- read.delim("rna_counts.txt", sep = "\t", header = TRUE)

#read sample-info
coldata <- read.delim("coldata.txt", sep = "\t")
coldata <-as.data.frame(apply(coldata, 2, as.factor))
head(coldata)

#changes to factor
coldata[,2] <-factor(coldata[,2])

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

#write results from DESeq2  into file
write.table(res, file = "significant_DESeq2results.txt", sep = "\t", row.names = FALSE)

#filters res by padj, keeping only those with padj < 0.05. then extracts rownames of those < 0.05 padj
resOrdered <- rownames(res[which(res$padj < 0.05),])
head(resOrdered)

#print LOC107985770 gene from res (checks resOrdered worked)
res[which(rownames(res) == "LOC107985770"),]

#read in nuORFs dataset from https://smorfs.ddnetbio.com/
ncORFs <- read.csv("/lyceum/cle1g21/smORFs/data-2024-10-25.csv")
head(ncORFs)

ncORFs$Gene_id

#filter resOrdered to only include genes in nuORFs
resOrdered_nuORFs <- resOrdered[which(resOrdered %in% ncORFs$Gene_id)]
resOrdered_nuORFs
str(resOrdered_nuORFs)

write.table(resOrdered_nuORFs, file = "matched_nuORFs.txt", sep = "\t", row.names = FALSE)

dds <- DESeqDataSetFromMatrix(countData = rna,
                              colData = coldata,
                              design = ~Condition)
dds

read.csv("/lyceum/cle1g21/smORFs/data-2024-10-25.csv")