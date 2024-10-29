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

#coldata[9:16, ] selects for RNA-seq samples only

#changes to factor. 
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

#filters res by padj, keeping only those with padj < 0.05. then extracts rownames of those < 0.05 padj
resOrdered <- rownames(res[which(res$padj < 0.05),])
head(resOrdered)

#print LOC107985770 gene from res (checks resOrdered worked)
res[which(rownames(res) == "LOC107985770"),]

#read in nuORFs dataset from https://smorfs.ddnetbio.com/
ncORFs <- read.csv("data-2024-10-25.csv")
head(ncORFs)

ncORFs$Gene_id

#filter resOrdered to only include genes in nuORFs
resOrdered_nuORFs <- resOrdered[which(resOrdered %in% ncORFs$Gene_id)]
resOrdered_nuORFs

################################################################# try using DeSeq2 with GEO datasett using PCA to check sample groups first
#retrieves GEO dataset GSE45516 using GEOquery package
gset <- getGEO("GSE45516", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1

gset <- gset[[idx]]
head(gset)

#make data.frame from gset for coldata
coldata <- pData(gset)
head(coldata)
type(coldata)

GSE45516_series_matrix.txt

file_path <- "/lyceum/cle1g21/smORFs/GSE45516_series_matrix.txt"

#reads GSE45516_series_matrix.txt from !series_matrix_table_begin line to !series_matrix_table_end line
lines <- readLines(file_path)
begin_idx <- grep("!series_matrix_table_begin", lines)[1]
end_idx <- grep("!series_matrix_table_end", lines)[1]
rna <- read.table(file_path, skip = begin_idx, nrows = end_idx - begin_idx - 1)

head(data)

coldata <- read.csv("/lyceum/cle1g21/smORFs/coldata", sep = "\t", header = TRUE)


####continue from here

#read countdata
countdata <- read.delim("GSE268366_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t", header = TRUE)
head(countdata)
#read coldata
coldata <- read.csv("GSE268366_coldata.csv", sep = "\t")
#renames coldata column X...Gen to GeneID
coldata <- rename(coldata, "GeneID" = "X...Gen")

#changes the second column in coldata to factors using factor() function. DESeq2 requires the design variables (columns in coldata) to be factors
coldata[,2] <-factor(coldata[,2])
head(coldata)
class(coldata$GeneID)
class(coldata$Condition)

#check that sample names match
all(coldata$GeneID == colnames(countdata))

colnames(countdata)
rownames(coldata)

head(coldata)

length(coldata$GeneID)
length(colnames(countdata))

rownames(coldata) <- coldata$GeneID







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
colnames(coldata)[1] <- "SampleID"

#remove GeneID column
gene_ids <- rna$GeneID

#matches geneIDS
rownames(rna) = rna$GeneID
rna_clean <- rna[, -1]
colnames(rna_clean)

#make DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = rna_clean,
                              colData = coldata,
                              design = ~Condition)
dds

#PCA
rld <- rlog(dds)
#png(filename = "pca.png", width = 600, height = 600) #save PCA plot
plotPCA(rld, intgroup = "Condition")
# dev.off()

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
