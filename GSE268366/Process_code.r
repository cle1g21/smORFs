# ** dataset GSE268366: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268366 **
# ** current dataset I am working on to perform DESeq2 analysis. checked PCA and it is acceptable, in the process of transforming the data to run DESeq2 **
# ** problem matching gene ID to res from DESeq2 to get names of significant genes to match to nuORFs **
library(GEOquery)
library(DESeq2)
library(dplyr)

#set working directory to appropriate dataset folder within smORFs
setwd("/lyceum/cle1g21/smORFs/GSE268366")

#read GEO dataset GSE268366 counts
cts <- as.matrix(read.csv("GSE268366_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t", header = TRUE))
rownames(cts) <- cts[, "GeneID"]
head(cts)

#read coldata
coldata <- read.table("GSE268366_coldata.csv", header = TRUE, sep = ",")

#changes the second column in coldata to factors using factor() function. DESeq2 requires the design variables (columns in coldata) to be factors
coldata[,2] <-factor(coldata[,2])
head(coldata)
class(coldata$GeneID)
class(coldata$Condition)

colnames(coldata) <- c("SampleID", "Condition")
#remove GeneID column
gene_ids <- cts[, "GeneID"]
cts_clean <- cts[, -1]

#sets the row names of coldata (geneIDs) to the column names of cts_clean (geneIDs)
rownames(coldata) <- colnames(cts_clean)


#check that sample names match
all(colnames(cts_clean) == rownames(coldata))
rownames(coldata) <- coldata$SampleID
rownames(coldata) 
colnames(cts_clean)

dds <- DESeqDataSetFromMatrix(countData = cts_clean,
                              colData = coldata,
                              design = ~Condition)
dds

#PCA plot to check replicates
rld <- rlog(dds)
#png(filename = "pca.png", width = 600, height = 600) #save PCA plot
plotPCA(rld, intgroup = "Condition")
#dev.off()

dds <- DESeq(dds)
res <- results(dds)
print(res)

summary(res)  # Provides a summary of the results object

#filters res by padj, keeping only those with padj < 0.05. then extracts rownames of those < 0.05 padj
resOrdered <- res[which(res$padj < 0.05),]
resOrdered$GeneID <- gene_ids #ERROR length of gene_ids and DESeq2 reuslts do not match
head(resOrdered)

#filter resOrdered to only include genes in nuORFs
resOrdered_nuORFs <- resOrdered[which(resOrdered %in% ncORFs$Gene_id)]
resOrdered_nuORFs
str(resOrdered_nuORFs)

write.table(resOrdered_nuORFs, file = "matched_nuORFs.txt", sep = "\t", row.names = FALSE)


####need to go back over
#changing entrez to ensembl

BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")
# For human genes

library(AnnotationDbi)
library(org.Hs.eg.db)
# Convert Entrez IDs to Ensembl IDs
ensembl_ids <- mapIds(
    org.Hs.eg.db,
    keys = rownames(res),
    column ="ENSEMBL",
    keytype ="ENTREZID",
    multiVals ="first"
    # Use "first" to handle multiple mappings
    )
    # Add Ensembl IDs as a new column in the DESeq2 results
res$Ensembl <- ensembl_ids
 
# Print the updated result
head(res)

res

write.table(res, file = "significant_DESeq2results.txt", sep = "\t", row.names = FALSE)

significant_DESeq2results <- read.table("significant_DESeq2results.txt", header = TRUE)
nrow(significant_DESeq2results)

significant_DESeq2results_matched_nuORFs <- read.table("matched_nuORFs.txt", header = TRUE)
nrow(significant_DESeq2results_matched_nuORFs)

significant_DESeq2results_matched_nuORFs$Ensembl
colnames(significant_DESeq2results_matched_nuORFs)

