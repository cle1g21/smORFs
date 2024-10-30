#** dataset GSE79666: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79666 **
#** dataset did not work. no genes with padj < 0.05 included in nuORFs from https://smorfs.ddnetbio.com/ **

library(GEOquery)
library(DESeq2)
library(dplyr)

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
png(filename = "pca.png", width = 600, height = 600) #save PCA plot
plotPCA(rld, intgroup = "Condition")
dev.off()

#make sure comparison is against control
dds$Condition <-relevel(dds$Condition, ref = "CTRL")

#perform DESeq2 analysis
dds <- DESeq(dds)
res <-results(dds, alpha = 0.05) 
res

#filters res by padj, keeping only those with padj < 0.05. then extracts rownames of those < 0.05 padj
resOrdered <- rownames(res[which(res$padj < 0.05),])
head(resOrdered)

#read in nuORFs dataset from https://smorfs.ddnetbio.com/
ncORFs <- read.csv("data-2024-10-25.csv")
head(ncORFs)

ncORFs$Gene_id

#filter resOrdered to only include genes in nuORFs
resOrdered_nuORFs <- resOrdered[which(resOrdered %in% ncORFs$Gene_id)]
resOrdered_nuORFs