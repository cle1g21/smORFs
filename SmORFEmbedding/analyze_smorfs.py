"""
analyze_smorfs.py
-----------------
Analyze smORFs and compare characteristics by primary set, for both ORF and full transcript.
Extracts features, performs statistical tests, and generates plots.
"""
#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from scipy.stats import mannwhitneyu
import re
from Bio import SeqIO

def gc_content(seq):
    """
    Calculate GC content of a sequence.
    Args:
        seq (str): Nucleotide sequence.
    Returns:
        float: GC content (0-1).
    """
    seq = seq.upper()
    gc = seq.count('G') + seq.count('C')
    return gc / len(seq) if len(seq) > 0 else 0

def sequence_entropy(seq):
    """
    Calculate Shannon entropy of a sequence.
    Args:
        seq (str): Nucleotide sequence.
    Returns:
        float: Entropy value.
    """
    seq = seq.upper()
    counts = Counter(seq)
    total = len(seq)
    entropy = -sum((c/total) * np.log2(c/total) for c in counts.values() if c > 0)
    return entropy

def has_kozak(seq):
    """
    Check for presence of a Kozak consensus sequence around ATG start codons.
    Args:
        seq (str): Nucleotide sequence.
    Returns:
        int: 1 if Kozak found, 0 otherwise.
    """
    seq = seq.upper()
    for m in re.finditer('ATG', seq):
        start = m.start()
        if start >= 6 and start+3+3 <= len(seq):
            context = seq[start-6:start+3+3]
            if re.match(r'GCC[AG]CCATGG', context):
                return 1
    return 0

def kmer_freq(seq, k=3):
    """
    Calculate k-mer frequencies for a sequence.
    Args:
        seq (str): Nucleotide sequence.
        k (int): k-mer size.
    Returns:
        dict: k-mer -> frequency
    """
    seq = seq.upper()
    total = len(seq) - k + 1
    counts = Counter([seq[i:i+k] for i in range(total)])
    for kmer in counts:
        counts[kmer] /= total
    return counts

def load_transcript_seqs(fasta_path):
    """
    Load transcript sequences from a FASTA file.
    Args:
        fasta_path (str): Path to FASTA file.
    Returns:
        dict: transcript_id -> sequence
    """
    seqs = {}
    for record in SeqIO.parse(fasta_path, 'fasta'):
        # Extract transcript_id=ENST... from the header
        m = re.search(r'transcript_id=([\w\.\-]+)', record.description)
        if m:
            tid = m.group(1).split('.')[0]  # match on base ID
            seqs[tid] = str(record.seq)
    return seqs

def extract_features(df, k=3, transcript_seqs=None):
    """
    Extract features for each smORF and (optionally) its transcript.
    Args:
        df (pd.DataFrame): smORF dataframe.
        k (int): k-mer size.
        transcript_seqs (dict): transcript_id -> sequence.
    Returns:
        pd.DataFrame: Feature matrix.
    """
    features = []
    all_kmers = set()
    for seq in df['sequence_nt']:
        all_kmers.update(kmer_freq(seq, k).keys())
    if transcript_seqs is not None:
        for tid in df['transcript']:
            tseq = transcript_seqs.get(tid.split('.')[0])
            if tseq:
                all_kmers.update(kmer_freq(tseq, k).keys())
    all_kmers = sorted(all_kmers)
    for _, row in df.iterrows():
        feats = {}
        seq = row['sequence_nt']
        feats['orf_length'] = len(seq)
        feats['orf_gc_content'] = gc_content(seq)
        feats['orf_entropy'] = sequence_entropy(seq)
        feats['orf_kozak'] = has_kozak(seq)
        kmer_counts = kmer_freq(seq, k)
        for kmer in all_kmers:
            feats[f'orf_kmer_{kmer}'] = kmer_counts.get(kmer, 0)
        # Transcript features
        tseq = transcript_seqs.get(row['transcript'].split('.')[0]) if transcript_seqs else None
        if tseq:
            feats['tx_length'] = len(tseq)
            feats['tx_gc_content'] = gc_content(tseq)
            feats['tx_entropy'] = sequence_entropy(tseq)
            feats['tx_kozak'] = has_kozak(tseq)
            t_kmer_counts = kmer_freq(tseq, k)
            for kmer in all_kmers:
                feats[f'tx_kmer_{kmer}'] = t_kmer_counts.get(kmer, 0)
        else:
            feats['tx_length'] = np.nan
            feats['tx_gc_content'] = np.nan
            feats['tx_entropy'] = np.nan
            feats['tx_kozak'] = np.nan
            for kmer in all_kmers:
                feats[f'tx_kmer_{kmer}'] = np.nan
        feats['primary_set'] = row['primary_set']
        features.append(feats)
    return pd.DataFrame(features)

def plot_and_test(features, outdir, k=3):
    """
    Generate violin plots and perform Mann-Whitney U tests for features by primary set.
    Args:
        features (pd.DataFrame): Feature matrix.
        outdir (str): Output directory for plots and stats.
        k (int): k-mer size.
    """
    os.makedirs(outdir, exist_ok=True)
    results = []
    for prefix in ['orf', 'tx']:
        for col in ['length', 'gc_content', 'entropy', 'kozak']:
            colname = f'{prefix}_{col}'
            if colname not in features:
                continue
            plt.figure(figsize=(8, 5))
            sns.violinplot(x=features['primary_set'], y=features[colname])
            plt.title(f'{col.replace("_", " ").title()} by Primary Set ({prefix})')
            plt.tight_layout()
            plt.savefig(os.path.join(outdir, f'{colname}_by_primary_set.png'))
            plt.close()
            # Mann-Whitney U test and effect size
            group1 = features[features['primary_set'] == 'yes'][colname].dropna()
            group2 = features[features['primary_set'] == 'no'][colname].dropna()
            if len(group1) > 0 and len(group2) > 0:
                stat, p = mannwhitneyu(group1, group2, alternative='two-sided')
                n1, n2 = len(group1), len(group2)
                # Rank-biserial correlation (effect size)
                effect_size = 1 - (2 * stat) / (n1 * n2)
                results.append(f'{colname}: U={stat:.2f}, p={p:.3e}, effect_size={effect_size:.3f}')
            else:
                results.append(f'{colname}: Not enough data for test')
        # Top k-mers
        kmer_cols = [col for col in features.columns if col.startswith(f'{prefix}_kmer_')]
        if kmer_cols:
            mean_kmers_yes = features[features['primary_set']=='yes'][kmer_cols].mean().sort_values(ascending=False)[:10]
            mean_kmers_no = features[features['primary_set']=='no'][kmer_cols].mean().sort_values(ascending=False)[:10]
            plt.figure(figsize=(10, 5))
            sns.barplot(x=mean_kmers_yes.values, y=mean_kmers_yes.index, color='blue', label='Primary Set')
            plt.title(f'Top 10 {k}-mer Frequencies (Primary Set, {prefix})')
            plt.xlabel('Mean Frequency')
            plt.tight_layout()
            plt.savefig(os.path.join(outdir, f'top10_{k}mer_primary_set_{prefix}.png'))
            plt.close()
            plt.figure(figsize=(10, 5))
            sns.barplot(x=mean_kmers_no.values, y=mean_kmers_no.index, color='red', label='Non-Primary Set')
            plt.title(f'Top 10 {k}-mer Frequencies (Non-Primary Set, {prefix})')
            plt.xlabel('Mean Frequency')
            plt.tight_layout()
            plt.savefig(os.path.join(outdir, f'top10_{k}mer_non_primary_set_{prefix}.png'))
            plt.close()
    # Save test results
    with open(os.path.join(outdir, 'statistical_tests.txt'), 'w') as f:
        for line in results:
            f.write(line + '\n')
    print('\n'.join(results))

def main():
    """
    Main workflow: loads smORF CSV, extracts features, generates plots, and performs statistical tests.
    """
    parser = argparse.ArgumentParser(description='Analyze smORFs and compare characteristics by primary set, for both ORF and full transcript.')
    parser.add_argument('-i', '--input', type=str, default='data/smorfs_v2.csv', help='Input smORF CSV file')
    parser.add_argument('-k', '--kmer', type=int, default=3, help='k-mer size')
    parser.add_argument('-o', '--outdir', type=str, default='plots/smorfs', help='Output directory for plots and stats')
    parser.add_argument('-t', '--transcript_fasta', type=str, default='data/all_transcripts.fasta', help='Transcriptome FASTA file')
    args = parser.parse_args()

    print(f'Loading {args.input} ...')
    df = pd.read_csv(args.input, low_memory=False)
    print(f'Loaded {len(df)} smORFs.')
    print(f'Loading transcriptome FASTA: {args.transcript_fasta} ...')
    transcript_seqs = load_transcript_seqs(args.transcript_fasta)
    print(f'Loaded {len(transcript_seqs)} transcript sequences.')

    # Log missing transcript IDs
    missing = []
    for tid in df['transcript']:
        tid_short = tid.split('.')[0]
        if tid_short not in transcript_seqs:
            missing.append(tid)
    os.makedirs(args.outdir, exist_ok=True)
    missing_path = os.path.join(args.outdir, 'missing_transcripts.txt')
    with open(missing_path, 'w') as f:
        for tid in missing:
            f.write(tid + '\n')
    print(f"{len(missing)} transcript IDs from smORFs did not match any sequence in the FASTA. See {missing_path}")

    print('Extracting features ...')
    features = extract_features(df, k=args.kmer, transcript_seqs=transcript_seqs)
    print('Plotting and testing ...')
    plot_and_test(features, args.outdir, k=args.kmer)
    print(f'Plots and stats saved to {args.outdir}')

if __name__ == '__main__':
    main() 