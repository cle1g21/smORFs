"""
setup_offline.py
----------------
Pre-downloads all required Hugging Face models/tokenizers and Ensembl reference files for offline use.
Run this script before moving your pipeline to an offline/HPC environment.
"""
import os
from transformers.utils.hub import snapshot_download
import requests

# --- Hugging Face models/tokenizers to pre-download ---
MODELS = [
    # DNABERT
    "zhihan1996/DNABERT-2-117M",
    "zhihan1996/DNABERT-6",
    # Nucleotide Transformer
    "InstaDeepAI/nucleotide-transformer-500m-human-ref",
    "InstaDeepAI/nucleotide-transformer-2.5b-multi-species",
    # ESM2
    "facebook/esm2_t6_8M_UR50D",
    # HyenaDNA
    "LongSafari/hyenadna-small-32k-seqlen-hf",
]

# --- Ensembl files to pre-download (as referenced in code) ---
ENSEMBL_FILES = [
    # GTFs
    ("https://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/Homo_sapiens.GRCh38.114.gtf.gz", "data/Homo_sapiens.GRCh38.114.gtf.gz"),
    # FASTA (genome)
    ("https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz", "data/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"),
    # cDNA and ncRNA (for transcript FASTA generation)
    ("https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz", "data/Homo_sapiens.GRCh38.cdna.all.fa.gz"),
    ("https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz", "data/Homo_sapiens.GRCh38.ncrna.fa.gz"),
]

def download_ensembl_file(url, dest):
    """
    Download a file from Ensembl FTP if it does not already exist locally.
    Args:
        url (str): URL to download from.
        dest (str): Local destination path.
    """
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    if os.path.exists(dest):
        print(f"[Ensembl] File {dest} already exists. Skipping download.")
        return
    print(f"[Ensembl] Downloading {url} ...")
    r = requests.get(url, stream=True)
    r.raise_for_status()
    with open(dest, 'wb') as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)
    print(f"[Ensembl] Downloaded to {dest}")

def download_hf_model(model_id):
    """
    Download and cache a Hugging Face model/tokenizer for offline use.
    Args:
        model_id (str): Hugging Face model repo ID.
    """
    print(f"[HF] Downloading and caching model/tokenizer: {model_id}")
    snapshot_download(repo_id=model_id, local_files_only=False, ignore_patterns=["*.safetensors"])
    print(f"[HF] Cached: {model_id}")

def main():
    """
    Download all required Hugging Face models and Ensembl files for offline use.
    """
    print("--- Downloading Hugging Face models/tokenizers ---")
    for model in MODELS:
        try:
            download_hf_model(model)
        except Exception as e:
            print(f"[HF] Error downloading {model}: {e}")

    print("\n--- Downloading Ensembl GTF/FASTA files ---")
    for url, dest in ENSEMBL_FILES:
        try:
            download_ensembl_file(url, dest)
        except Exception as e:
            print(f"[Ensembl] Error downloading {url}: {e}")

    print("\nAll downloads complete. You should now be able to run the pipeline offline.")

if __name__ == "__main__":
    main() 