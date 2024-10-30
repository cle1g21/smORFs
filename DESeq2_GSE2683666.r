# ** dataset GSE2683666: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2683666 **
# ** current dataset I am working on to perform DESeq2 analysis. checked PCA and it is acceptable, in the process of transforming the data to run DESeq2 **

library(GEOquery)
library(DESeq2)
library(dplyr)

#read GEO dataset GSE2683666 counts
cts <- as.matrix(read.csv("GSE268366_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t", header = TRUE))
head(cts)

#read coldata
coldata <- read.csv("GSE268366_coldata.csv", sep = "\t")
