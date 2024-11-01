#** Datasets from https://github.com/SGDDNB/translational_regulation/tree/master/sample_data **
#** Replicate Owens findings using DESeq2 and match to nuORFs (filter by nc rna) from https://smorfs.ddnetbio.com/ **

library(GEOquery)
library(DESeq2)
library(dplyr)

#set working directory to appropriate dataset folder within smORFs
setwd("/lyceum/cle1g21/smORFs/deltaTE_Owen")

#read Owens counts. read.delim returns a data.frame
cts <- read.delim("rna_counts.txt", sep = "\t", header = TRUE)

#read coldata
coldata <- read.delim("coldata.txt", sep = "\t")

#changes the second column in coldata to factors using factor() function
coldata[,2] <-factor(coldata[,2])
class(coldata$Condition)

#sets the row names of coldata (geneIDs) to the column names of cts (geneIDs)
rownames(coldata) <- colnames(cts)

#check that sample names match
all(colnames(rna) == rownames(coldata))
colnames(rna)
rownames(coldata)

#make DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~Condition)
dds

#PCA plot to check replicates
rld <- rlog(dds)
#png(filename = "pca.png", width = 600, height = 600) #save PCA plot
plotPCA(rld, intgroup = "Condition")
#dev.off()

#run DESeq
dds <- DESeq(dds)
res <- results(dds)
print(res) # prints the first few rows of the DESeq2 results
summary(res)  # provides a summary of the DESeq2 results 

#write results from DESeq2  into file
write.table(res, file = "DESeq2results.txt", sep = "\t", row.names = FALSE)

#filters res by padj, keeping only those with padj < 0.05. then extracts rownames of those < 0.05 padj
resOrdered <- res[which(res$padj < 0.05), ]
head(resOrdered)

#read in nuORFs dataset from https://smorfs.ddnetbio.com/
ncORFs <- read.csv("/lyceum/cle1g21/smORFs/data-2024-10-25.csv")
head(ncORFs)
colnames(resOrdered)

#filter resOrdered to only include genes in nuORFs. use rownames() if the gene ID is the row names and gene ID does not show as a column name in resOrdered
resOrdered_nuORFs <- resOrdered[rownames(resOrdered) %in% ncORFs$Gene_id, ]

# sets the row names of the resulting data frame to the original gene IDs. otherwsie the above line will replace the gene ID row names with standard number row names
rownames(resOrdered_nuORFs) <- rownames(resOrdered)[rownames(resOrdered) %in% ncORFs$Gene_id]

resOrdered_nuORFs
str(resOrdered_nuORFs)

# save file of matching nuORFs. row.names = TRUE writes the row names into the table (preserves the gene IDs as row names)
write.table(resOrdered_nuORFs, file = "nuORFs_matched.txt", sep = "\t", row.names = TRUE)

#inspect the DESeq2 results and the resOrdered_nuORFs data frame

DESeq2results <- read.table("DESeq2results.txt", header = TRUE)
nrow(DESeq2results)

nuORFs_matched <- read.table("nuORFs_matched.txt", header = TRUE)
nrow(nuORFs_matched)
head(nuORFs_matched)