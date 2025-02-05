# > 30 minutes terminated DEGSeq2 before completed

library(GEOquery)
library(DESeq2)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)

cts <- read.csv("/mainfs/ddnb/TCGA_data/full_TCGA-BRCA_matrix.csv")
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
dds <- DESeq(dds) #  > 30 minutes terminated before completed