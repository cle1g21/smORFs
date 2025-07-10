#!/usr/bin/env python3
"""
embed_protein_fasta_transformer.py
----------------------------------
Embeds all known protein sequences from Ensembl and smORF protein sequences from smorfs_v2.csv using ProtTrans models (ProtBERT, ProtT5).
Saves mean and max pooled embeddings and metadata to HDF5.

Features:
- Downloads Ensembl protein FASTA if not present
- Extracts unique protein sequences from Ensembl and smORFs
- Embeds using selected ProtTrans model
- Saves embeddings and metadata
- Visualises UMAP (Ensembl vs smORF)
- Finds nearest Ensembl proteins for each smORF in embedding space
- Test mode: only 100 smORFs and 100 Ensembl proteins
"""
import os
import argparse
import gzip
import h5py
import numpy as np
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
import requests
import torch
from transformers import AutoTokenizer, AutoModel
import umap
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.neighbors import NearestNeighbors

ENSEMBL_PEP_URL = "https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz"
ENSEMBL_PEP_PATH = "data/Homo_sapiens.GRCh38.pep.all.fa.gz"
SMORF_CSV_PATH = "data/smorfs_v2.csv"

MODEL_CONFIGS = {
    'protbert': {
        'default_model': 'Rostlab/prot_bert',
        'max_len': 1024,
        'tokenize': lambda seq: ' '.join(list(seq)),
        'preprocess': lambda seq: seq.upper().replace('U', 'X'),
    },
    'prott5': {
        'default_model': 'Rostlab/prot_t5_xl_uniref50',
        'max_len': 1024,
        'tokenize': lambda seq: ' '.join(list(seq)),
        'preprocess': lambda seq: seq.upper().replace('U', 'X'),
    },
}

def download_ensembl_pep():
    """Download Ensembl protein FASTA if not present."""
    os.makedirs(os.path.dirname(ENSEMBL_PEP_PATH), exist_ok=True)
    if os.path.exists(ENSEMBL_PEP_PATH):
        print(f"[Ensembl] Protein FASTA already exists: {ENSEMBL_PEP_PATH}")
        return
    print(f"[Ensembl] Downloading protein FASTA from {ENSEMBL_PEP_URL} ...")
    r = requests.get(ENSEMBL_PEP_URL, stream=True)
    r.raise_for_status()
    with open(ENSEMBL_PEP_PATH, 'wb') as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)
    print(f"[Ensembl] Downloaded to {ENSEMBL_PEP_PATH}")

def parse_ensembl_pep(fasta_path):
    """
    Parse Ensembl protein FASTA and return dict: {protein_id: sequence}
    Args:
        fasta_path (str): Path to gzipped Ensembl protein FASTA
    Returns:
        dict: protein_id -> sequence
    """
    seqs = {}
    with gzip.open(fasta_path, 'rt') as f:
        for record in SeqIO.parse(f, 'fasta'):
            seqs[record.id] = str(record.seq)
    return seqs

def get_smorfs_protein_seqs_with_ids(csv_path):
    """
    Extract (releasev45_id, sequence_aa_MS) pairs from smorfs_v2.csv.
    Args:
        csv_path (str): Path to smorfs_v2.csv
    Returns:
        list of [releasev45_id, sequence_aa_MS]
    """
    df = pd.read_csv(csv_path, usecols=lambda c: c.strip().startswith('sequence_aa_MS') or c.strip() == 'releasev45_id')
    col = [c for c in df.columns if c.strip().startswith('sequence_aa_MS')][0]
    seqs = df[["releasev45_id", col]].dropna().astype(str)
    seqs[col] = seqs[col].str.strip()
    return seqs.drop_duplicates().values.tolist()

def embed_sequences(seqs, model_type, model_name=None, max_len=None, batch_size=8, device=None):
    """
    Embed a list of protein sequences using the selected ProtTrans model.
    Args:
        seqs (list): List of protein sequences (str)
        model_type (str): 'protbert' or 'prott5'
        model_name (str): Optional Hugging Face model ID
        max_len (int): Max tokens per sequence
        batch_size (int): Batch size for embedding
        device (str): 'cuda' or 'cpu'
    Returns:
        mean_embeds (np.ndarray): Mean pooled embeddings
        max_embeds (np.ndarray): Max pooled embeddings
        meta (list): List of sequences (str)
    """
    config = MODEL_CONFIGS[model_type]
    model_name = model_name or config['default_model']
    max_len = max_len or config['max_len']
    preprocess = config['preprocess']
    tokenize = config['tokenize']
    print(f"Loading model: {model_name}")
    tokenizer = AutoTokenizer.from_pretrained(model_name, do_lower_case=False)
    model = AutoModel.from_pretrained(model_name)
    device = device or ('cuda' if torch.cuda.is_available() else 'cpu')
    model = model.to(device)
    mean_embeds = []
    max_embeds = []
    meta = []
    for i in tqdm(range(0, len(seqs), batch_size), desc="Embedding proteins"):
        batch = seqs[i:i+batch_size]
        batch_proc = [preprocess(s) for s in batch]
        batch_tok = [tokenize(s) for s in batch_proc]
        inputs = tokenizer(batch_tok, return_tensors='pt', padding=True, truncation=True, max_length=max_len)
        inputs = {k: v.to(device) for k, v in inputs.items()}
        with torch.no_grad():
            outputs = model(**inputs)
            if hasattr(outputs, 'last_hidden_state'):
                emb = outputs.last_hidden_state.cpu().numpy()
            elif hasattr(outputs, 'logits'):
                emb = outputs.logits.cpu().numpy()
            else:
                raise ValueError("Model output does not have last_hidden_state or logits.")
        # Pool over sequence length (axis=1)
        mean_embeds.append(emb.mean(axis=1))
        max_embeds.append(emb.max(axis=1))
        meta.extend(batch)
    mean_embeds = np.concatenate(mean_embeds, axis=0)
    max_embeds = np.concatenate(max_embeds, axis=0)
    return mean_embeds, max_embeds, meta

def plot_umap(embeddings, labels, out_path):
    """
    Visualise embeddings in 2D UMAP, coloured by label (Ensembl or smORF).
    Args:
        embeddings (np.ndarray): Embedding matrix
        labels (list): List of 'Ensembl' or 'smORF'
        out_path (str): Output PNG path
    """
    reducer = umap.UMAP(random_state=42)
    umap_emb = reducer.fit_transform(embeddings)
    plt.figure(figsize=(10, 8))
    sns.scatterplot(x=umap_emb[:,0], y=umap_emb[:,1], hue=labels, s=10, alpha=0.7, legend='full')
    plt.title('UMAP of Protein Embeddings (Ensembl vs smORF)')
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()
    print(f"UMAP plot saved to {out_path}")

def create_nearest_neighbor_network(mean_embeds, meta, smorf_set, ensembl_set, out_path, k=5, smorf_id_map=None, ensembl_id_map=None):
    """
    For each smORF, find the k nearest Ensembl proteins in embedding space.
    Save as CSV: smORF_sequence, smORF_releasev45_id, ensembl_sequence, ensembl_protein_id, distance
    Args:
        mean_embeds (np.ndarray): Mean embeddings for all sequences
        meta (list): List of sequences (str)
        smorf_set (set): Set of smORF sequences
        ensembl_set (set): Set of Ensembl protein sequences
        out_path (str): Output CSV path
        k (int): Number of nearest neighbours
        smorf_id_map (dict): smORF sequence -> releasev45_id
        ensembl_id_map (dict): Ensembl sequence -> protein_id
    """
    smorf_indices = [i for i, seq in enumerate(meta) if seq in smorf_set]
    ensembl_indices = [i for i, seq in enumerate(meta) if seq in ensembl_set]
    if not smorf_indices or not ensembl_indices:
        print("No smORF or Ensembl sequences found for network analysis.")
        return
    smorf_embeds = mean_embeds[smorf_indices]
    ensembl_embeds = mean_embeds[ensembl_indices]
    ensembl_seqs = [meta[i] for i in ensembl_indices]
    smorf_seqs = [meta[i] for i in smorf_indices]
    nn = NearestNeighbors(n_neighbors=min(k, len(ensembl_embeds)), metric='euclidean')
    nn.fit(ensembl_embeds)
    distances, indices = nn.kneighbors(smorf_embeds)
    rows = []
    for i, smorf_seq in enumerate(smorf_seqs):
        smorf_id = smorf_id_map.get(smorf_seq, '') if smorf_id_map else ''
        for j in range(indices.shape[1]):
            ensembl_seq = ensembl_seqs[indices[i, j]]
            ensembl_id = ensembl_id_map.get(ensembl_seq, '') if ensembl_id_map else ''
            dist = distances[i, j]
            rows.append({'smORF_sequence': smorf_seq, 'smORF_releasev45_id': smorf_id, 'ensembl_sequence': ensembl_seq, 'ensembl_protein_id': ensembl_id, 'distance': dist})
    df = pd.DataFrame(rows)
    df.to_csv(out_path, index=False)
    print(f"Nearest-neighbor network saved to {out_path}")

def main():
    """
    Main workflow:
    - Download Ensembl protein FASTA if needed
    - Parse Ensembl and smORF protein sequences
    - (Test mode: only 100 smORFs and 100 Ensembl proteins)
    - Embed all unique sequences
    - Save mean/max embeddings and metadata
    - Plot UMAP (Ensembl vs smORF)
    - Find k nearest Ensembl proteins for each smORF
    """
    parser = argparse.ArgumentParser(description="Embed all known protein sequences from Ensembl and smORFs using ProtTrans models.")
    parser.add_argument('--type', type=str, choices=['protbert', 'prott5'], default='protbert', help='Model type')
    parser.add_argument('--model', type=str, help='Model checkpoint (overrides default)')
    parser.add_argument('--outdir', type=str, default='embeddings', help='Output directory')
    parser.add_argument('--test', action='store_true', help='Embed only 100 smORFs and 100 Ensembl proteins for testing')
    args = parser.parse_args()

    # Download Ensembl protein FASTA if needed
    download_ensembl_pep()
    print("Parsing Ensembl protein FASTA ...")
    ensembl_proteins = parse_ensembl_pep(ENSEMBL_PEP_PATH)
    print(f"Loaded {len(ensembl_proteins)} Ensembl protein sequences.")
    print("Extracting smORF protein sequences ...")
    smorf_id_seqs = get_smorfs_protein_seqs_with_ids(SMORF_CSV_PATH)
    smorf_proteins = [seq for _, seq in smorf_id_seqs]
    smorf_id_map = {seq: rid for rid, seq in smorf_id_seqs}
    # Ensembl protein ID map
    ensembl_id_map = {v: k for k, v in ensembl_proteins.items()}
    # TEST MODE: select only 100 smORFs and 100 Ensembl proteins
    if args.test:
        print("[TEST MODE] Only embedding 100 smORFs and 100 Ensembl proteins.")
        smorf_id_seqs = smorf_id_seqs[:100]
        smorf_proteins = [seq for _, seq in smorf_id_seqs]
        smorf_id_map = {seq: rid for rid, seq in smorf_id_seqs}
        ensembl_items = list(ensembl_proteins.items())[:100]
        ensembl_proteins = dict(ensembl_items)
        ensembl_id_map = {v: k for k, v in ensembl_proteins.items()}
    # Combine all unique protein sequences
    all_proteins = list(set(ensembl_proteins.values()) | set(smorf_proteins))
    print(f"Total unique protein sequences: {len(all_proteins)}")
    os.makedirs(args.outdir, exist_ok=True)
    out_path = os.path.join(args.outdir, f"protein_{args.type}_embeddings.h5" if not args.test else f"test_protein_{args.type}_embeddings.h5")
    mean_embeds, max_embeds, meta = embed_sequences(all_proteins, args.type, args.model)
    # Save embeddings and metadata
    with h5py.File(out_path, 'w') as h5:
        h5.create_dataset('mean', data=mean_embeds)
        h5.create_dataset('max', data=max_embeds)
        h5.create_dataset('sequence', data=np.array(meta, dtype='S'))
    print(f"Saved embeddings to {out_path}")
    # UMAP visualisation
    smorf_set = set(smorf_proteins)
    ensembl_set = set(ensembl_proteins.values())
    labels = [('smORF' if seq in smorf_set else 'Ensembl') for seq in meta]
    umap_path = os.path.join(args.outdir, f"protein_{args.type}_umap.png" if not args.test else f"test_protein_{args.type}_umap.png")
    plot_umap(mean_embeds, labels, umap_path)
    # Nearest-neighbour network
    network_path = os.path.join(args.outdir, f"protein_{args.type}_smorf_ensembl_network.csv" if not args.test else f"test_protein_{args.type}_smorf_ensembl_network.csv")
    create_nearest_neighbor_network(mean_embeds, meta, smorf_set, ensembl_set, network_path, k=5, smorf_id_map=smorf_id_map, ensembl_id_map=ensembl_id_map)

if __name__ == "__main__":
    main() 