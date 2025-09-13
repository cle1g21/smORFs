# .csv files for matched peptides for:
# all nuorf lncrnas differentially expressed in BRCA
# top 15 upregulated nuorf lncrnas differentially expressed in BRCA
# all upregulated nuorf lncrnas differentially expressed in BRCA
# Created FASTA files from the .csv files

# Load required packages
library(dplyr)
library(biomaRt)

# peptides for all nuorf lncrnas differentially expressed in BRCA--------------------------------------------------------------------------------------------------------------------
# Load the files
nuorf_atlas <- read.csv("/lyceum/cle1g21/smORFs/nuORF_atlas.csv")
nuorf_lncrna_matched <- read.delim("/lyceum/cle1g21/smORFs/TCGA/BRCA/Outputs/lncRNAs_nuORFs_matched_BRCA.txt", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

View(nuorf_atlas)
View(nuorf_lncrna_matched)

# Convert rownames into a new column
nuorf_lncrna_matched$Gene_id <- rownames(nuorf_lncrna_matched)

# Creates a merged dataframe of lncRNAs that are nuORFs differentially expressed in BRCA with their predicted peptide sequences
merged_data <- nuorf_lncrna_matched %>%
  left_join(nuorf_atlas[, c("Gene_id", "Peptide.seq")], by = "Gene_id")

# Checks merged_data has worked correctly. Should return TRUE
"ENSG00000230838" %in% merged_data$Gene_id
View(merged_data)

merged_data$Gene_id <- make.unique(merged_data$Gene_id)

write.csv(merged_data, "/lyceum/cle1g21/smORFs/TCGA/BRCA/Outputs/FASTA_files/nuORF_lncRNAs_with_peptides.csv", row.names = FALSE)

# peptides for top 15 upregulated nuorf lncrnas differentially expressed in BRCA-----------------------------------------------------------------------------------------------------
nuORF_lncRNAs_with_peptides <- read.csv("/lyceum/cle1g21/smORFs/TCGA/BRCA/Outputs/FASTA_files/nuORF_lncRNAs_with_peptides.csv")
top15_nuorf_lncrna_degs <- read.csv("/lyceum/cle1g21/smORFs/TCGA/BRCA/Outputs/top15upreg_lncRNA_matched_nuORFs.csv")

View(nuORF_lncRNAs_with_peptides)
View(top15_nuorf_lncrna_degs)

# Filter to keep only genes in the top15 list
filtered_top15_nuorf_lncRNAs <- nuORF_lncRNAs_with_peptides %>%
  filter(Gene_id %in% top15_nuorf_lncrna_degs$Gene_ID)

View(filtered_top15_nuorf_lncRNAs)
# Step 3: Save the filtered data to a new CSV
write.csv(filtered_top15_nuorf_lncRNAs, "/lyceum/cle1g21/smORFs/TCGA/BRCA/Outputs/FASTA_files/filtered_top15_nuorf_lncRNAs_with_peptides.csv", row.names = FALSE)

# peptides for upregulated nuorf lncrnas differentially expressed in BRCA------------------------------------------------------------------------------------------------------------
upregulated_nuorf_lncrnas <- nuORF_lncRNAs_with_peptides %>%
  filter(padj < 0.05, log2FoldChange > 0)

View(upregulated_nuorf_lncrnas)
# Should return TRUE
"ENSG00000230838" %in% nuORF_lncRNAs_with_peptides$Gene_id
"ENSG00000230838" %in% upregulated_nuorf_lncrnas$Gene_id

View(nuORF_lncRNAs_with_peptides)

# Save the file
write.csv(upregulated_nuorf_lncrnas, "/lyceum/cle1g21/smORFs/TCGA/BRCA/Outputs/FASTA_files/upregulated_nuorf_lncRNAs_with_peptides.csv", row.names = FALSE)

# Create FASTA files from .csv files------------------------------------------------------------------------------------------------------------
# For all nuorf lncrnas differentially expressed in BRCA

# Filter for non-empty peptide sequences
all_nuORFs <- merged_data %>% filter(Peptide.seq != "" & !is.na(Peptide.seq))

# Open a file connection for writing
fasta_file <- file("/lyceum/cle1g21/smORFs/TCGA/BRCA/Outputs/FASTA_files/nuORF_lncRNAs_with_peptides.fasta", open = "w")

# Loop over each row and write to the file in FASTA format
for (i in 1:nrow(all_nuORFs)) {
  cat(paste0(">", all_nuORFs$Gene_id[i], "\n", all_nuORFs$Peptide.seq[i], "\n"), file = fasta_file)
}

# Close the file
close(fasta_file)

# For top 15 upregulated nuorf lncrnas differentially expressed in BRCA

# Filter for non-empty peptide sequences
filtered_top15_nuorf_lncRNAs <- filtered_top15_nuorf_lncRNAs %>% filter(Peptide.seq != "" & !is.na(Peptide.seq))

# Open a file connection for writing
fasta_file <- file("/lyceum/cle1g21/smORFs/TCGA/BRCA/Outputs/FASTA_files/filtered_top15_nuorf_lncRNAs_with_peptides.fasta", open = "w")

# Loop over each row and write to the file in FASTA format
for (i in 1:nrow(filtered_top15_nuorf_lncRNAs)) {
  cat(paste0(">", filtered_top15_nuorf_lncRNAs$Gene_id[i], "\n", filtered_top15_nuorf_lncRNAs$Peptide.seq[i], "\n"), file = fasta_file)
}

# Close the file
close(fasta_file)

# For all upregulated nuorf lncrnas differentially expressed in BRCA

# Filter for non-empty peptide sequences
upregulated_nuorf_lncrnas <- upregulated_nuorf_lncrnas %>% filter(Peptide.seq != "" & !is.na(Peptide.seq))

# Open a file connection for writing
fasta_file <- file("/lyceum/cle1g21/smORFs/TCGA/BRCA/Outputs/FASTA_files/upregulated_nuorf_lncRNAs_with_peptides.fasta", open = "w")

# Loop over each row and write to the file in FASTA format
for (i in 1:nrow(upregulated_nuorf_lncrnas)) {
  cat(paste0(">", upregulated_nuorf_lncrnas$Gene_id[i], "\n", upregulated_nuorf_lncrnas$Peptide.seq[i], "\n"), file = fasta_file)
}

# Close the file
close(fasta_file)

# make human_ref fasta file same as nuORFs -> remove pipe and change to ensembl ids
human_ref_nuORFs <- readLines("/lyceum/cle1g21/smORFs/TCGA/BRCA/Outputs/FASTA_files/nuORF_humanRef_BRCA.fasta")

for (i in seq_along(human_ref_nuORFs)) {
  if (startsWith(human_ref_nuORFs[i], ">")) {
    # Check if it's a UniProt header (starts with tr| or sp|, both are UniProt)
    if (grepl("^>tr\\|[^|]+\\|", human_ref_nuORFs[i]) || grepl("^>sp\\|[^|]+\\|", human_ref_nuORFs[i])) {
      # Extract only the UniProt ID
      uniprot_id <- sub("^>[^|]+\\|([^|]+)\\|.*", "\\1", human_ref_nuORFs[i])
      human_ref_nuORFs[i] <- paste0(">", uniprot_id)
    }
    # Otherwise leave header as-is (e.g., Ensembl)
  }
}

View(human_ref_nuORFs)

# Write to a new FASTA file
writeLines(human_ref_nuORFs, "/lyceum/cle1g21/smORFs/TCGA/BRCA/Outputs/FASTA_files/nuORF_humanRef_BRCA_modified.fasta")
