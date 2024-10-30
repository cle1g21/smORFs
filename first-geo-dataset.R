
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

#remove GeneID column
gene_ids <- countdata$GeneID
countdata_clean <- countdata[, -1]
colnames(countdata_clean)

#check that sample names match
colnames(coldata)[1] <- "GeneID"
colnames(coldata)
head(coldata)
colnames(countdata_clean)

identical(colnames(countdata_clean), rownames(coldata)) #this worked
rownames(coldata) <- colnames(countdata_clean) #this worked

coldata <- data.frame(Condition = coldata$Condition)
rownames(coldata) <- rownames(coldata)  # Ensuring row names are still sample IDs

coldata <- coldata[, "Condition", drop = FALSE]


dds <- DESeqDataSetFromMatrix(countData = countdata_clean,
                              colData = coldata,
                              design = ~Condition)
dds

#PCA
rld <- rlog(dds)
#png(filename = "pca.png", width = 600, height = 600) #save PCA plot
plotPCA(rld, intgroup = "Condition")
#dev.off()

dds <- DESeq(dds)
res <- results(dds)
print(res)

#filters res by padj, keeping only those with padj < 0.05. then extracts rownames of those < 0.05 padj
resOrdered <- rownames(res[which(res$padj < 0.05),])
res_sig <- res[which(res$padj < 0.05),]
res_sig <- as.data.frame(res_sig)          # Convert to data frame
res_sig$GeneID <- rownames(res_sig)        # Add GeneID column with row names
head(res_sig)  # View the top rows

head(resOrdered)

min(res$padj, na.rm = TRUE)


















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
png(filename = "pca.png", width = 600, height = 600) #save PCA plot
plotPCA(rld, intgroup = "Condition")
dev.off()

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
