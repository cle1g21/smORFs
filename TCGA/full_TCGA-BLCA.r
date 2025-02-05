# 15 mins to run and DESeq2:
#-- replacing outliers and refitting for 9480 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 
#-- original counts are preserved in counts(dds)

library(GEOquery)
library(DESeq2)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)

cts <- read.csv("/mainfs/ddnb/TCGA_data/full_TCGA-BLCA_matrix.csv")
#View(cts)

# Setup cts ----
# Set the row names to Gene name
rownames(cts) <- cts[, "X"]

# Sort gene name in gene_ids object
gene_ids <- cts[, "X"]
# Remove gene name column, row names are gene names
cts_clean <- cts[, -1]

# Check this has worked
#View(cts_clean)

# Create coldata (metadata)
# Take column names (sample) into a object
sample_ids <- colnames(cts[, -1])
# Create data frame with sample
coldata <- data.frame(SampleID = sample_ids)

# Setup coldata ----
# TCGA (The Cancer Genome Atlas) data stores condition in 5th and 6th characters of the barcode's sample
# Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29 http://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/

# Extract the sample type code (5th and 6th characters from SampleID) 
coldata$SampleType <- substr(coldata$SampleID, 14, 15)

# Map sample types to condition type. %02d makes sure always 2 digits (adds a 0 if only 1 digit)
coldata$Condition <- ifelse(coldata$SampleType %in% sprintf("%02d", 1:9), "Tumor",       # 01-09
                     ifelse(coldata$SampleType %in% sprintf("%02d", 10:19), "Normal",    # 10-19
                     ifelse(coldata$SampleType %in% sprintf("%02d", 20:29), "Control",  # 20-29
                     "Other")))  # For any other sample types


# Remove the SampleType column to keep only SampleID and Condition
coldata <- coldata[, c("SampleID", "Condition")]

# Check the updated coldata
print(coldata)
#View(coldata)

#changes the second column in coldata to factors. DESeq2 requires the design variables (columns in coldata) to be factors
coldata[,2] <-factor(coldata[,2])
class(coldata$Condition)

#sets the row names of coldata (geneIDs) to the column names of cts_clean (geneIDs)
rownames(coldata) <- colnames(cts_clean)

# DESeq2 ----
#check that sample names match
all(colnames(cts_clean) == rownames(coldata))
rownames(coldata) 
colnames(cts_clean)

# Create matrix
dds <- DESeqDataSetFromMatrix(countData = cts_clean,
                              colData = coldata,
                              design = ~Condition)
dds

#PCA plot to check replicates
#rld <- rlog(dds)
#png(filename = "pca.png", width = 600, height = 600) #save PCA plot
#plotPCA(rld, intgroup = "Condition")
#dev.off()

# Run DESeq2
dds <- DESeq(dds) # 15 minutes to run
res <- results(dds)
print(res) # prints the first few rows of the DESeq2 results
summary(res)  # provides a summary of the DESeq2 results 

# creates a file of the DESeq2 results
write.table(res, file = "/lyceum/cle1g21/smORFs/TCGA/DESeq2results_BLCA.txt", sep = "\t", row.names = TRUE)

#filters res by padj, keeping the rows from the res data frame where the padj value is less than 0.05
resOrdered <- res[which(res$padj < 0.05), ]
head(resOrdered)

#read in nuORFs dataset from https://smorfs.ddnetbio.com/
ncORFs <- read.csv("/lyceum/cle1g21/smORFs/nuORF_atlas.csv")
head(ncORFs)

# filter ncORFs to remove rows where NAs are found in rownames (gene name) of resOrdered and saved to a dataset called noNA_resOrdered. need to do this because the next line gives an error if NAs are found
noNA_resOrdered <- resOrdered[!is.na(rownames(resOrdered)), ]

# subset the resOrdered data frame to include only the rows where the Ensembl column is in ncORFs$Gene_id. so the resOrdered only include genes in the nuORFs list
resOrdered_nuORFs <- noNA_resOrdered[rownames(noNA_resOrdered) %in% ncORFs$Gene_id, ]
print(resOrdered_nuORFs)
str(resOrdered_nuORFs)

# save file of matching nuORFs
write.table(resOrdered_nuORFs, file = "/lyceum/cle1g21/smORFs/TCGA/nuORFs_matched_BLCA.txt", sep = "\t", row.names = TRUE)
