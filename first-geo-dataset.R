library(GEOquery)
library(DESeq2)

library(dplyr)

#read Owens raw counts and sample-info
rna <- read.delim("rna_counts.txt", sep = "\t", header = TRUE)

coldata <- read.delim("sample_info.txt", sep = "\t")
coldata <-as.data.frame(apply(coldata, 2, as.factor))
head(coldata)

all(colnames(rna) == rownames(coldata))
colnames(rna)
rownames(coldata)

dds <- DESeqDataSetFromMatrix(countData = rna,
                              colData = coldata,
                              design = ~Condition)


#changed column names to match design_GSE273142
colnames(rna) <- c("GeneID", "GSM8422516", "GSM8422517", "GSM8422518", "GSM8422519", "GSM8422520", "GSM8422521")

colnames(coldata)|>
    select(c("Pao24RIBO24hpA5R1001trimmedL29NoAbundant.fastq.gzAligned.sortedByCoord.out.bam", "Pao27RIBO24hpA10R1001trimmedL29NoAbundant.fastq.gzAligned.sortedByCoord.out.bam	2", "Pao34RIBO24hpA15R1001trimmedL29NoAbundant.fastq.gzAligned.sortedByCoord.out.bam", "Pao41RIBO24hpA20R1001trimmedL29NoAbundant.fastq.gzAligned.sortedByCoord.out.bam", "Pao24RIBOBSLpA1R1001trimmedL29NoAbundant.fastq.gzAligned.sortedByCoord.out.bam", "Pao27RIBOBSLpA6R1001trimmedL29NoAbundant.fastq.gzAligned.sortedByCoord.out.bam", "Pao34RIBOBSLpA11R1001trimmedL29NoAbundant.fastq.gzAligned.sortedByCoord.out.bam", "Pao41RIBOBSLpA16R1001trimmedL29NoAbundant.fastq.gzAligned.sortedByCoord.out.bam"))

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
