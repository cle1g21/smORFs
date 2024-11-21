# FAIL no nuORFs found in DESeq2 results of SRP012015 (M1 and M2)
library(recount3)
library(SummarizedExperiment)
library(DESeq2)
library(tidyverse)

# RECOUNT3 BIOCONDUCTOR TUTORIAL -----------------------------------------------
human_projects <- available_projects()

## Extract the RangedSummarizedExperiment ----
# find the project of interest
proj_info <- subset(
  human_projects,
  project == "SRP012015" & project_type == "data_sources"
)

# create RangedSummarizedExperiment object for the gene level expression data
rse_gene_SRP012015 <- create_rse(proj_info)

### Start note -----------------------------------------------------------------
# the above generates same RangedSummarizedExperiment as copying code from
# recount3 exploration tool after looking up SRP id (found on GEO)
# SRP012015.recount3 <- recount3::create_rse_manual(
#   project = "SRP012015",
#   project_home = "data_sources/sra",
#   organism = "human",
#   annotation = "gencode_v26",
#   type = "gene"
# )
### End note -------------------------------------------------------------------

# use transform_counts() on RangedSummarizedExperiment to transform the raw 
# coverage counts

assayNames(rse_gene_SRP012015) #check name is "counts"
#assay(rse_gene_SRP012015, "counts") <-
#  assay(rse_gene_SRP012015, "raw_counts") #renames "raw_counts" to "counts"

#transforms counts
transformed_counts <- transform_counts(rse_gene_SRP012015)
#replaces "counts" assay in RangedSummarizedExperiment with the transformed data
assay(rse_gene_SRP012015, "raw_counts") <- transformed_counts

# number of genes by number of samples
dim(rse_gene_SRP012015)
# information about the genes
rowRanges(rse_gene_SRP012015)

# DDS ----
## Prepare colData and create condition ----
# find condition in colData of RangedSummarizedExperiment
rse_gene_SRP012015@colData
colnames(colData(rse_gene_SRP012015)) #list all column names in the metadata

# here the conditions (M1 or M2) are included in sra.sample_attributes [24]
rse_gene_SRP012015@colData[24]

# create new column named "condition" matching each sample ID to M1 or M2
colnames(rse_gene_SRP012015) #lists sample ID to match order for each condition
colData(rse_gene_SRP012015)$condition <- c("M1", "M1", "M1", "M2", "M2", "M2")

# change "condition" to factor for DESeq2
colData(rse_gene_SRP012015)$condition <-
  factor(colData(rse_gene_SRP012015)$condition)

# check everything looks correct before generating DESeq2
# colnames(colData(rse_gene_SRP012015)) can be used to extract index number
rse_gene_SRP012015@colData[176] #index for "condition"

## Create colData and countData for DESeqDataSetFromMatrix() ----
# create coldata as a data frame
coldata <- as.data.frame(colData(rse_gene_SRP012015))
coldata$condition <- factor(coldata$condition) #make sure "condition" is a factor

# extract the counts matrix for countData
counts <- assay(rse_gene_SRP012015)

## Generate dds object for DESeq2 using countData and colData ----
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~condition)

#PCA plot to check replicates
rld <- rlog(dds)
#png(filename = "pca.png", width = 600, height = 600) #save PCA plot
plotPCA(rld, intgroup = "condition")
#dev.off()

## Run DESeq2 analysis
dds <- DESeq(dds)
res <- results(dds)
print(res) #prints the first few rows of the DESeq2 results
summary(res)  #provides a summary of the DESeq2 results 

# creates a file of the DESeq2 results for GSE268366 dataset
write.table(res, file = "DESeq2results.txt", sep = "\t", row.names = TRUE)

# filters res by padj, keeping only those with padj < 0.05. then extracts 
# rownames of those < 0.05 padj
resOrdered <- res[which(res$padj < 0.05), ]
head(resOrdered)

# MATCH RES TO NUORFS ----------------------------------------------------------
#read in nuORFs dataset from https://smorfs.ddnetbio.com/
ncORFs <- read.csv("data-2024-11-20.csv")
head(ncORFs)
colnames(resOrdered)

#filter resOrdered to only include genes in nuORFs. use rownames() if the gene
# ID is the row names and gene ID does not show as a column name in resOrdered
resOrdered_nuORFs <- resOrdered[rownames(resOrdered) %in% ncORFs$Gene_id, ]

# sets the row names of the resulting data frame to the original gene IDs.
# otherwise the above line will replace the gene ID row names with standard
# number row names
rownames(resOrdered_nuORFs) <-
  rownames(resOrdered)[rownames(resOrdered) %in% ncORFs$Gene_id]

resOrdered_nuORFs
str(resOrdered_nuORFs)

# save file of matching nuORFs
write.table(resOrdered_nuORFs, file = "nuORFs_matched.txt", sep = "\t", row.names = FALSE)
