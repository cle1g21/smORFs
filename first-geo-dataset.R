library(GEOquery)
library(DESeq2)
library(dplyr)

#read Owens raw counts and sample-info
counts <- read.csv("rna_counts.txt", sep = "\t", header = TRUE)
df <- as.data.frame(counts)

info <- read.csv("sample_info.txt", sep = "\t")

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = info,
                              design = ~Condition)





#reading the .tab file and transforming so all the values are in separate columns
tab <- read.delim("GSE273142.tab", sep = "\t", header = TRUE, strip.white = TRUE)
#changing tab file downloaded from GEO to a .csv
write.csv(tab, file = "GSE273142.csv", row.names = FALSE)
head(tab)

#changed column names to match design_GSE273142
colnames(tab) <- c("GeneID", "GSM8422516", "GSM8422517", "GSM8422518", "GSM8422519", "GSM8422520", "GSM8422521")

#create data.frame
as.data.frame(tab)

#read design .csv 
design <- read.csv("design_GSE273142.csv", header = TRUE)

#removes GeneID column ready for DeSeq2
gene_ids <- tab$GeneID
count_table_clean <- tab[, -1]


#construct DeSeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_table_clean,
                              colData = design_clean,
                              design = ~condition)

rownames(design_clean)
rownames(count_table_clean)
