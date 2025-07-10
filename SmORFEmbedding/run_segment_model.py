"""
run_segment_model.py
-------------------
Extracts flanking regions for smORF-containing transcripts, runs a token classification model (e.g., Nucleotide Transformer), and saves conservation predictions for each region.
"""
import pandas as pd
import gtfparse
import pyfaidx
from tqdm import tqdm
import torch
from transformers import AutoTokenizer, AutoModelForTokenClassification

def get_flanking_regions(gtf_path, smorf_transcripts, genome_fa, upstream=1000, downstream=1000):
    """
    Extract upstream and downstream flanking regions for a list of transcripts from a GTF and genome FASTA.

    Args:
        gtf_path (str): Path to GTF file.
        smorf_transcripts (list): List of transcript IDs.
        genome_fa (str): Path to genome FASTA file.
        upstream (int): Number of bases upstream to extract.
        downstream (int): Number of bases downstream to extract.

    Returns:
        dict: {"<transcript_id>_upstream": sequence, "<transcript_id>_downstream": sequence}
    """
    print("Parsing GTF file...")
    gtf = gtfparse.read_gtf(gtf_path)
    transcripts_gtf = gtf[gtf['feature'] == 'transcript']
    flanking_seqs = {}
    print("Loading genome...")
    genome = pyfaidx.Fasta(genome_fa)
    print(f"Processing {len(smorf_transcripts)} transcripts...")
    for transcript_id in tqdm(smorf_transcripts):
        transcript_info = transcripts_gtf[transcripts_gtf['transcript_id'] == transcript_id]
        if transcript_info.empty:
            print(f"Warning: Transcript {transcript_id} not found in GTF file.")
            continue
        transcript_info = transcript_info.iloc[0]
        chrom = transcript_info['seqname']
        start = int(transcript_info['start'])
        end = int(transcript_info['end'])
        strand = transcript_info['strand']
        # Calculate flanking region coordinates based on strand
        if strand == '+':
            up_start = start - 1 - upstream
            up_end = start - 1
            down_start = end
            down_end = end + downstream
        else: # strand == '-'
            up_start = end
            up_end = end + upstream
            down_start = start - 1 - downstream
            down_end = start - 1
        # pyfaidx uses 1-based indexing for start/end
        # Handle chromosome names like 'chr1' vs '1'
        chrom_name_in_fasta = f"chr{chrom}" if "chr" not in chrom else chrom
        if chrom_name_in_fasta not in genome:
            chrom_name_in_fasta = chrom.replace("chr", "")
        if chrom_name_in_fasta not in genome:
             print(f"Warning: Chromosome {chrom} not found in FASTA file.")
             continue
        try:
            if strand == '+':
                up_seq = genome[chrom_name_in_fasta][up_start:up_end].seq
                down_seq = genome[chrom_name_in_fasta][down_start:down_end].seq
            else: # On minus strand, the "upstream" region is downstream in coordinate terms
                  # and we need the reverse complement
                up_seq = genome[chrom_name_in_fasta][up_start:up_end].reverse.complement.seq
                down_seq = genome[chrom_name_in_fasta][down_start:down_end].reverse.complement.seq
            flanking_seqs[f"{transcript_id}_upstream"] = up_seq
            flanking_seqs[f"{transcript_id}_downstream"] = down_seq
        except pyfaidx.FastaIndexingError as e:
            print(f"Error extracting sequence for {transcript_id}: {e}")
    return flanking_seqs

def predict_conservation(sequences, model, tokenizer):
    """
    Predict conservation scores for a batch of sequences using a token classification model.

    Args:
        sequences (list): List of DNA sequences (str).
        model: Hugging Face model for token classification.
        tokenizer: Hugging Face tokenizer.

    Returns:
        torch.Tensor: Mean logits per sequence (conservation score).
    """
    inputs = tokenizer(sequences, return_tensors="pt", padding=True, truncation=True, max_length=1000)
    with torch.no_grad():
        outputs = model(**inputs)
    # The output format will depend on the model. 
    # For token classification, it's usually logits.
    # We can average them to get a per-sequence score.
    return outputs.logits.mean(dim=1)

def main():
    """
    Main workflow: loads smORF transcript list, extracts flanking regions, runs model, and saves results.
    """
    # --- Configuration ---
    SMORF_CSV_PATH = 'data/smorfs_v2.csv'
    GTF_PATH = 'data/Homo_sapiens.GRCh38.114.gtf.gz'
    # IMPORTANT: You need to provide the path to the human genome FASTA file.
    # You can download it from Ensembl or GENCODE.
    # e.g. https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/
    GENOME_FASTA_PATH = 'path/to/your/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
    MODEL_ID = "InstaDeepAI/nucleotide-transformer-2.5b-multi-species"
    # --- Script ---
    print("Loading smORF data...")
    smorfs_df = pd.read_csv(SMORF_CSV_PATH)
    smorf_transcripts = smorfs_df['transcript'].unique()
    print(f"Found {len(smorf_transcripts)} unique transcripts with smORFs.")
    # For testing, let's just use a small subset
    smorf_transcripts_subset = smorf_transcripts[:10]
    flanking_regions = get_flanking_regions(GTF_PATH, smorf_transcripts_subset, GENOME_FASTA_PATH)
    if not flanking_regions:
        print("No flanking regions could be extracted. Exiting.")
        return
    print("Loading model and tokenizer...")
    tokenizer = AutoTokenizer.from_pretrained(MODEL_ID)
    model = AutoModelForTokenClassification.from_pretrained(MODEL_ID)
    print("Running predictions...")
    results = {}
    # Process sequences in batches for efficiency
    sequences_to_process = list(flanking_regions.values())
    sequence_ids = list(flanking_regions.keys())
    batch_size = 4 
    for i in tqdm(range(0, len(sequences_to_process), batch_size)):
        batch_ids = sequence_ids[i:i+batch_size]
        batch_seqs = sequences_to_process[i:i+batch_size]
        predictions = predict_conservation(batch_seqs, model, tokenizer)
        for seq_id, pred in zip(batch_ids, predictions):
            results[seq_id] = pred.numpy()
    print("Saving results...")
    results_df = pd.DataFrame.from_dict(results, orient='index')
    results_df.to_csv('flanking_region_conservation.csv')
    print("Done. Results saved to flanking_region_conservation.csv")

if __name__ == "__main__":
    main() 