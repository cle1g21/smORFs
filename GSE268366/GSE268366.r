# ** dataset GSE268366: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268366 **
# ** current dataset I am working on to perform DESeq2 analysis. checked PCA and it is acceptable, in the process of transforming the data to run DESeq2 **
# ** problem matching gene ID to res from DESeq2 to get names of significant genes to match to nuORFs **
# file structure:
#          cleaning datasets check colnames(cts_clean) == rownames(coldata) ->   
#          DESeqDataSetFromMatrix ->
#          PCA ->
#          DESeq ->
#          save file of results ->
#          map significant genes to nuORFs ->
#          save file of matching nuORFs

library(GEOquery)
library(DESeq2)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)

#set working directory to appropriate dataset folder within smORFs
setwd("/lyceum/cle1g21/smORFs/GSE268366")

#read GEO dataset GSE268366 counts
cts <- as.matrix(read.csv("counts.tsv", sep = "\t", header = TRUE))

#sets the row names of the cts data frame to be the values in the "GeneID" column. This removes the row numbers and sets them to the geneIDs. common practice in bioinformatics 
rownames(cts) <- cts[, "GeneID"]
head(cts)

#remove GeneID column as these are now the row names in cts (also stored in gene_ids)
gene_ids <- cts[, "GeneID"]
cts_clean <- cts[, -1]

#read coldata
coldata <- read.table("coldata.csv", header = TRUE, sep = ",")

#changes the second column in coldata to factors using factor() function. DESeq2 requires the design variables (columns in coldata) to be factors
coldata[,2] <-factor(coldata[,2])
class(coldata$Condition)

#renames the columns in the coldata data frame to "SampleID" and "Condition"
colnames(coldata) <- c("SampleID", "Condition")
head(coldata)

#sets the row names of coldata (geneIDs) to the column names of cts_clean (geneIDs)
rownames(coldata) <- colnames(cts_clean)

#check that sample names match
all(colnames(cts_clean) == rownames(coldata))
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
print(res) # prints the first few rows of the DESeq2 results
summary(res)  # provides a summary of the DESeq2 results 

# this dataset uses Entrez IDs instead of Ensembl IDs. this is necessary to convert so that the nuORFs can be matched to the DESeq2 results
# maps Entrez IDs to Ensembl IDs using the org.Hs.eg.db database
ensembl_ids <- mapIds(        #mapIds() function returns a vector of Ensembl IDs that correspond to the Entrez IDs in the res data frame
    org.Hs.eg.db,             # org.Hs.eg.db database used to map Entrez IDs to Ensembl IDs for human genes
    keys = rownames(res),     # specifies that the entrez IDs are the row names of the DESeq2 results (res) data frame    
    column ="ENSEMBL",        # specifies that the Ensembl IDs should be retrieved from the "ENSEMBL" column of org.Hs.eg.db database
    keytype ="ENTREZID",      # specifies that the row names of res are Entrez IDs
    multiVals ="first"        # specifies that if there are multiple Ensembl IDs associated with a single Entrez ID, only the first one should be returned
    )  

# adds the mapped Ensembl IDs as a new column in the DESeq2 results (res)
res$Ensembl <- ensembl_ids
head(res)

# creates a file of the DESeq2 results for GSE268366 dataset
write.table(res, file = "DESeq2results.txt", sep = "\t", row.names = TRUE)

#filters res by padj, keeping the rows from the res data frame where the padj value is less than 0.05
resOrdered <- res[which(res$padj < 0.05), ]
head(resOrdered)

#read in nuORFs dataset from https://smorfs.ddnetbio.com/
ncORFs <- read.csv("/lyceum/cle1g21/smORFs/data-2024-10-25.csv")
head(ncORFs)

# filter ncORFs to remove rows where NAs are found in Ensembl column of resOrdered and saved to a dataset called noNA_resOrdered. need to do this because the next line gives an error if NAs are found
noNA_resOrdered <- resOrdered[!is.na(resOrdered$Ensembl), ]

# subset the resOrdered data frame to include only the rows where the Ensembl column is in ncORFs$Gene_id. so the resOrdered only include genes in the nuORFs list
resOrdered_nuORFs <- noNA_resOrdered[noNA_resOrdered$Ensembl %in% ncORFs$Gene_id, ]
resOrdered_nuORFs
str(resOrdered_nuORFs)

# save file of matching nuORFs
write.table(resOrdered_nuORFs, file = "nuORFs_matched.txt", sep = "\t", row.names = FALSE)

#inspect the DESeq2 results and the resOrdered_nuORFs data frame

DESeq2results <- read.table("DESeq2results.txt", header = TRUE)
nrow(DESeq2results)

nuORFs_matched <- read.table("nuORFs_matched.txt", header = TRUE)
nrow(nuORFs_matched)
nuORFs_matched$Ensembl