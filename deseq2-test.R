####DESeq2 workflow tutorial

install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("airway")
install.packages("tidyverse")
library(DESeq2)
library(tidyverse)
library(airway)

#script to get data from airway package for tutorial on GitHub
data(airway)
airway

sample_info <- as.data.frame(colData(airway))
sample_info <- sample_info[,c(2,3)]
sample_info$dex <- gsub('trt', 'treated', sample_info$dex)
sample_info$dex <- gsub('untrt', 'untreated', sample_info$dex)
names(sample_info) <- c('cellLine', 'dexamethasone')
write.table(sample_info, file = "sample_info.csv", sep = ',', col.names = T, row.names = T, quote = F)

countsData <- assay(airway)
write.table(countsData, file = "counts_data.csv", sep = ',', col.names = T, row.names = T, quote = F)

#read in count data
counts_data <- read.csv('counts_data.csv')
head(counts_data)

colData <- read.csv('sample_info.csv')

#make sure the row names in colData matches column names in counts_data, and make sure they are in the same order otherwise DeSeq will error
all(colnames(counts_data) %in% rownames(colData)) #same name as output TRUE
all(colnames(counts_data) == rownames(colData)) #same order as output TRUE

#construct DeSeq dataset
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~dexamethasone) #colData has only 1 design factor (dexamethasone), giving treated vs untreated samples

dds #shows it is a DESeqDataSet, with 63677 rows and 8 columns, rownames as genes(ENG0...)

#perform pre-filtering DeSeq dataset object
#remove rows with low gene counts across the samples (helps reduce size of DeSeq object and increase speed of computation)
keep <- rowSums(counts(dds)) >= 10 #FALSE on rows that do not have more than 10 genes across the samples
dds <- dds[keep,] #remove rows that have less than 10 genes across the samples
dds

#set factor level. comapre treated vs untreated (these are factors). want to tell DeSeq to use untreated as the reference level
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")

dds$dexamethasone #shows now 2 levels (untreated and treated)

###NOTE: if the dataset has technical replicates you will need to collapse them before the DeSeq. **Never** collapse the biological replicates

#run DeSeq function
dds <- DESeq(dds)
res <-results(dds)
res #log2fold change calculated from the treated vs untreated. means whatever value for logfold change is for the treated compared to the untreated positive log2fold change = upregulated genes in the treated condition
#basemean average of normalized counts over all samples
#lfcSE standard estimates for the log2foldchange
#stat are the world test value for each gene
#padj is the corrected p value for multiple testing. need to correct because p value 0.05 so 5% of DEGs due to random change (false positive -> drug has no real effect in treated)

#explore results
summary(res) #shows how many genes are upregulated(LFCup) and downregulated(LFCdown), currently using padj 0.1

#can change the padj. this changes the values of LFC upreg and downreg
res0.01 <- results(dds, alpha = 0.01)
res

resultsNames(dds) #shows summary of comparison of treated vs untreated with dexamethasone
#if you have multiple levels to the dataset eg. treated_4hrs, treated_8hrs, untreated. want to compare each with untreated reference level:
#results(dds, contrast = c("dexamethasone", "treated_4hrs", "untreated")) #where first level will be the design factor (dexamethasome), and the second the level you want to compare
#allows comparison between levels. we do not need to use this here, as we are only comapring 1 level to the reference

#visualization of data eg. volcano plot or MA plot to visualize up and downregulation of genes                                                                    
#MA plot: scatterplot of logfold change compared to the mean of normalized counts. shows differentially expressed genes, with those highlighted differentially expressed with padj 0.05. smaller triangles higher fold changes and direction of arrow shows direction of fold change.want to see points in upper or lower right showing high logfold change




