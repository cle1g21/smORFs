#!/usr/bin/env python3
import os
import argparse
import h5py
import numpy as np
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
import torch
from transformers import AutoTokenizer, AutoModel
import umap
import matplotlib.pyplot as plt
import seaborn as sns
import re

MODEL_CONFIGS = {
    'dnabert': {
        'default_model': 'zhihan1996/DNABERT-2-117M',
        'k': 6,
        'max_len': 512,
        'tokenize': lambda seq, k: [seq[i:i+k] for i in range(len(seq)-k+1)],
        'preprocess': lambda seq: seq.upper().replace('U', 'T'),
    },
    'nucleotide-transformer': {
        'default_model': 'InstaDeepAI/nucleotide-transformer-500m-human-ref',
        'k': 6,
        'max_len': 512,
        'tokenize': lambda seq, k: [seq[i:i+k] for i in range(len(seq)-k+1)],
        'preprocess': lambda seq: seq.upper().replace('U', 'T'),
    },
    'esm2': {
        'default_model': 'facebook/esm2_t6_8M_UR50D',
        'k': 1,
        'max_len': 1022,
        'tokenize': lambda seq, k: list(seq),
        'preprocess': lambda seq: seq.upper().replace('U', 'T'),
    },
    'hyenadna': {
        'default_model': 'LongSafari/hyenadna-small-32k-seqlen-hf',
        'k': 1,
        'max_len': 32768,
        'tokenize': lambda seq, k: list(seq),
        'preprocess': lambda seq: seq.upper().replace('U', 'T'),
    },
}

def fasta_to_df(fasta_path):
    """
    Parse a FASTA file and return a DataFrame with transcript_id, biotype, and sequence.

    Args:
        fasta_path (str): Path to the input FASTA file.

    Returns:
        pd.DataFrame: DataFrame with columns ['transcript_id', 'biotype', 'seq']
    """
    records = []
    for record in SeqIO.parse(fasta_path, 'fasta'):
        header = record.description
        m = re.search(r'biotype=([^\\s]+)', header)
        biotype = m.group(1) if m else 'NA'
        m_tid = re.search(r'transcript_id=([\w\.\-]+)', header)
        transcript_id = m_tid.group(1) if m_tid else record.id
        records.append({
            'transcript_id': transcript_id,
            'biotype': biotype,
            'seq': str(record.seq)
        })
    return pd.DataFrame(records)

def scan_and_embed(seq, tokenizer, model, device, model_type, k, max_len):
    """
    Embed a sequence by splitting it into windows, tokenizing, and running through the model.
    Returns mean and max pooled embeddings for the sequence.

    Args:
        seq (str): Input sequence (DNA/RNA string)
        tokenizer: Hugging Face tokenizer
        model: Hugging Face model
        device: torch.device
        model_type (str): Model type key (e.g. 'dnabert', 'nucleotide-transformer', 'esm2', 'hyenadna')
        k (int): k-mer size (or 1 for models that don't use k-mers)
        max_len (int): Maximum tokens per window

    Returns:
        mean_pool (np.ndarray): Mean pooled embedding for the sequence
        max_pool (np.ndarray): Max pooled embedding for the sequence
    """
    config = MODEL_CONFIGS[model_type]
    seq = config['preprocess'](seq)
    tokens = config['tokenize'](seq, k)
    all_embeds = []
    for start in range(0, len(tokens), max_len):
        window = tokens[start:start+max_len]
        if len(window) < k:
            continue
        if model_type in ['dnabert', 'nucleotide-transformer']:
            input_str = ' '.join(window)
        elif model_type in ['esm2', 'hyenadna']:
            input_str = ''.join(window)
        else:
            raise ValueError(f"Unknown model_type: {model_type}")
        inputs = tokenizer(input_str, return_tensors='pt')
        inputs = {k: v.to(device) for k, v in inputs.items()}
        with torch.no_grad():
            outputs = model(**inputs)
            if hasattr(outputs, 'last_hidden_state'):
                emb = outputs.last_hidden_state.squeeze(0).cpu().numpy()
            elif hasattr(outputs, 'logits'):
                emb = outputs.logits.squeeze(0).cpu().numpy()
            else:
                raise ValueError("Model output does not have last_hidden_state or logits.")
            all_embeds.append(emb)
    if not all_embeds:
        return None, None
    all_embeds = np.concatenate(all_embeds, axis=0)
    mean_pool = np.mean(all_embeds, axis=0)
    max_pool = np.max(all_embeds, axis=0)
    return mean_pool, max_pool

def main():
    """
    Main workflow for embedding FASTA sequences using a selected transformer model.
    Handles argument parsing, model loading, embedding, saving results, and UMAP plotting.
    """
    parser = argparse.ArgumentParser(description='Embed FASTA sequences using DNABERT, Nucleotide Transformer, or ESM-2.')
    parser.add_argument('-i', '--input', type=str, default='data/all_transcripts.fasta', help='Input FASTA file')
    parser.add_argument('-o', '--outdir', type=str, default='embeddings', help='Output directory for embeddings')
    parser.add_argument('--type', type=str, choices=['dnabert', 'nucleotide-transformer', 'esm2', 'hyenadna'], default='dnabert', help='Model type')
    parser.add_argument('--model', type=str, help='Model checkpoint (overrides default for type)')
    parser.add_argument('--k', type=int, help='k-mer size (overrides default for type)')
    parser.add_argument('--max_len', type=int, help='Max tokens per window (overrides default for type)')
    parser.add_argument('--test', action='store_true', help='If set, only embed the first 100 sequences for testing')
    args = parser.parse_args()

    # Ensure output directory exists
    os.makedirs(args.outdir, exist_ok=True)
    basename = os.path.splitext(os.path.basename(args.input))[0]
    test_prefix = "test_" if args.test else ""
    h5_path = os.path.join(args.outdir, f'{test_prefix}{basename}_{args.type}_embeddings.h5')
    config = MODEL_CONFIGS[args.type]
    model_name = args.model or config['default_model']
    k = args.k or config['k']
    max_len = args.max_len or config['max_len']

    # If embeddings already exist, load them
    if os.path.exists(h5_path):
        print(f"Embeddings already exist at {h5_path}. Loading...")
        with h5py.File(h5_path, 'r') as h5:
            mean_embeds = h5['mean'][:]
            max_embeds = h5['max'][:]
            meta = pd.DataFrame({k: h5[k][:].astype(str) for k in ['transcript_id', 'biotype']})
        # Try to load all window embeddings if the separate file exists
        all_window_embeds = None
        all_window_biotypes = None
        all_window_ids = None
        all_h5_path = os.path.join(args.outdir, f'{test_prefix}{basename}_{args.type}_all_windows.h5')
        if os.path.exists(all_h5_path):
            with h5py.File(all_h5_path, 'r') as h5w:
                all_window_embeds = h5w['embeddings'][:]
                all_window_biotypes = h5w['biotypes'][:]
                all_window_ids = h5w['window_ids'][:]
                # window_starts and window_ends are also available if needed
    else:
        print(f"Loading sequences from {args.input} ...")
        df = fasta_to_df(args.input)
        if args.test:
            print("[TEST MODE] Only embedding the first 100 sequences.")
            df = df.iloc[:100]
        print(f"Loaded {len(df)} sequences.")
        print(f"Loading model: {model_name} ...")
        tokenizer = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
        model = AutoModel.from_pretrained(model_name, trust_remote_code=True)
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        model = model.to(device)
        mean_embeds = []
        max_embeds = []
        meta = []
        all_window_embeds = []
        all_window_biotypes = []
        all_window_ids = []
        all_window_starts = []
        all_window_ends = []
        for idx, row in tqdm(df.iterrows(), total=len(df), desc='Embedding sequences'):
            # For each sequence, split into windows, embed, and collect window-level and pooled embeddings
            config = MODEL_CONFIGS[args.type]
            seq = config['preprocess'](row['seq'])
            tokens = config['tokenize'](seq, k)
            window_embeds = []
            for start in range(0, len(tokens), max_len):
                window = tokens[start:start+max_len]
                if len(window) < k:
                    continue
                if args.type in ['dnabert', 'nucleotide-transformer']:
                    input_str = ' '.join(window)
                elif args.type in ['esm2', 'hyenadna']:
                    input_str = ''.join(window)
                else:
                    raise ValueError(f"Unknown model_type: {args.type}")
                inputs = tokenizer(input_str, return_tensors='pt')
                inputs = {k: v.to(device) for k, v in inputs.items()}
                with torch.no_grad():
                    outputs = model(**inputs)
                    if hasattr(outputs, 'last_hidden_state'):
                        emb = outputs.last_hidden_state.squeeze(0).cpu().numpy()
                    elif hasattr(outputs, 'logits'):
                        emb = outputs.logits.squeeze(0).cpu().numpy()
                    else:
                        raise ValueError("Model output does not have last_hidden_state or logits.")
                    window_embeds.append(emb)
                    all_window_embeds.append(emb.mean(axis=0))  # mean pool for UMAP
                    all_window_biotypes.append(row['biotype'])
                    all_window_ids.append(f"{row['transcript_id']}_window{start//max_len}")
                    all_window_starts.append(start)
                    all_window_ends.append(start + max_len)
            if window_embeds:
                window_embeds_concat = np.concatenate(window_embeds, axis=0)
                mean_pool = np.mean(window_embeds_concat, axis=0)
                max_pool = np.max(window_embeds_concat, axis=0)
                mean_embeds.append(mean_pool)
                max_embeds.append(max_pool)
                meta.append((row['transcript_id'], row['biotype']))
        mean_embeds = np.stack(mean_embeds)
        max_embeds = np.stack(max_embeds)
        meta = pd.DataFrame(meta, columns=['transcript_id', 'biotype'])
        # Save mean and max pooled embeddings for each sequence
        with h5py.File(h5_path, 'w') as h5:
            h5.create_dataset('mean', data=mean_embeds)
            h5.create_dataset('max', data=max_embeds)
            h5.create_dataset('transcript_id', data=meta['transcript_id'].astype('S'))
            h5.create_dataset('biotype', data=meta['biotype'].astype('S'))
        print(f"Saved embeddings to {h5_path}")
        # Save all window embeddings and metadata to a separate HDF5 file
        all_h5_path = os.path.join(args.outdir, f'{test_prefix}{basename}_{args.type}_all_windows.h5')
        all_window_embeds = np.stack(all_window_embeds)
        all_window_biotypes = np.array(all_window_biotypes, dtype='S')
        all_window_ids = np.array(all_window_ids, dtype='S')
        all_window_starts = np.array(all_window_starts)
        all_window_ends = np.array(all_window_ends)
        with h5py.File(all_h5_path, 'w') as h5w:
            h5w.create_dataset('embeddings', data=all_window_embeds)
            h5w.create_dataset('biotypes', data=all_window_biotypes)
            h5w.create_dataset('window_ids', data=all_window_ids)
            h5w.create_dataset('window_starts', data=all_window_starts)
            h5w.create_dataset('window_ends', data=all_window_ends)
        print(f"Saved all window embeddings to {all_h5_path}")

    # UMAP visualization (mean embedding)
    print("Running UMAP and plotting (mean embedding)...")
    reducer = umap.UMAP(random_state=42)
    umap_emb = reducer.fit_transform(mean_embeds)
    meta['biotype'] = meta['biotype'].astype(str)
    plt.figure(figsize=(10, 8))
    sns.scatterplot(x=umap_emb[:,0], y=umap_emb[:,1], hue=meta['biotype'], s=10, alpha=0.7, legend='full')
    plt.title(f'UMAP of {args.type} Mean Embeddings (colored by biotype)')
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, f'{test_prefix}{basename}_{args.type}_umap.png'))
    plt.close()
    print(f"UMAP plot saved to {args.outdir}/{test_prefix}{basename}_{args.type}_umap.png")

    # UMAP visualization (all window embeddings)
    if all_window_embeds is not None and len(all_window_embeds) > 0:
        print("Running UMAP and plotting (all window embeddings)...")
        reducer = umap.UMAP(random_state=42)
        umap_emb = reducer.fit_transform(all_window_embeds)
        all_window_biotypes_str = [b.decode() if isinstance(b, bytes) else b for b in all_window_biotypes]
        plt.figure(figsize=(10, 8))
        sns.scatterplot(x=umap_emb[:,0], y=umap_emb[:,1], hue=all_window_biotypes_str, s=10, alpha=0.7, legend='full')
        plt.title(f'UMAP of {args.type} Window Embeddings (colored by biotype)')
        plt.tight_layout()
        plt.savefig(os.path.join(args.outdir, f'{test_prefix}{basename}_{args.type}_all_windows_umap.png'))
        plt.close()
        print(f"UMAP plot (all windows) saved to {args.outdir}/{test_prefix}{basename}_{args.type}_all_windows_umap.png")

if __name__ == '__main__':
    main() 