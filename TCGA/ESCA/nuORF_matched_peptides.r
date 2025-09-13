# Load required packages
library(dplyr)
library(biomaRt)

# peptides for all nuorf lncrnas differentially expressed in ESCA--------------------------------------------------------------------------------------------------------------------
# Load the files
nuorf_atlas <- read.csv("/lyceum/cle1g21/smORFs/nuORF_atlas.csv")
nuorf_lncrna_matched <- read.delim("/lyceum/cle1g21/smORFs/TCGA/ESCA/Outputs/lncRNAs_nuORFs_matched_ESCA.txt", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

View(nuorf_atlas)
View(nuorf_lncrna_matched)
rownames(nuorf_lncrna_matched)

# Convert rownames into a new column
nuorf_lncrna_matched$Gene_id <- rownames(nuorf_lncrna_matched)

# Creates a merged dataframe of lncRNAs that are nuORFs differentially expressed in BRCA with their predicted peptide sequences
merged_data <- nuorf_lncrna_matched %>%
  left_join(nuorf_atlas[, c("Gene_id", "Peptide.seq")], by = "Gene_id")

# Checks merged_data has worked correctly. Should return TRUE
"ENSG00000230838" %in% merged_data$Gene_id
View(merged_data)

merged_data$Gene_id <- make.unique(merged_data$Gene_id)

write.csv(merged_data, "/lyceum/cle1g21/smORFs/TCGA/ESCA/Outputs/FASTA_files/nuORF_lncRNAs_with_peptides_ESCA.csv", row.names = FALSE)

# Create FASTA files from .csv files------------------------------------------------------------------------------------------------------------
# For all nuorf lncrnas differentially expressed in BRCA

# Filter for non-empty peptide sequences
all_nuORFs <- merged_data %>% filter(Peptide.seq != "" & !is.na(Peptide.seq))

# Open a file connection for writing
fasta_file <- file("/lyceum/cle1g21/smORFs/TCGA/ESCA/Outputs/FASTA_files/nuORF_lncRNAs_with_peptides_E_peaks2.fasta", open = "w")

# # Loop over each row and write to the file in FASTA format
# for (i in 1:nrow(all_nuORFs)) {
#   cat(paste0(">", all_nuORFs$Gene_id[i], "\n", all_nuORFs$Peptide.seq[i], "\n"), file = fasta_file)
# }

for (i in 1:nrow(all_nuORFs)) {
  gene_id <- all_nuORFs$Gene_id[i]
  peptide <- all_nuORFs$Peptide.seq[i]
  
  # Compose UniProt-like custom header for PEAKS2
  fasta_header <- paste0(
    ">tr|XX", gene_id, "|XX", gene_id, "_HUMAN nuORF ",
    "XX=Homo sapiens OX=9606 GN=", gene_id, " PE=1 SV=1"
  )
  # Write to file
  cat(fasta_header, "\n", peptide, "\n", sep = "", file = fasta_file)
}

# Close the file
close(fasta_file)

# make human_ref fasta file same as nuORFs -> remove pipe and change to ensembl ids
human_ref_nuORFs <- readLines("/lyceum/cle1g21/smORFs/TCGA/ESCA/Outputs/FASTA_files/nuORF_humanRef_ESCA.fasta")

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
writeLines(human_ref_nuORFs, "/lyceum/cle1g21/smORFs/TCGA/ESCA/Outputs/FASTA_files/nuORF_humanRef_ESCA_modified.fasta")
