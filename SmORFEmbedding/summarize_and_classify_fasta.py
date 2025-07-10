"""
summarize_and_classify_fasta.py
------------------------------
Summarizes sequence statistics and classifies FASTA sequences by biotype. Extracts features, visualizes distributions, and performs machine learning classification.
"""
#!/usr/bin/env python3
import argparse
import re
from collections import Counter, defaultdict
from Bio import SeqIO
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import confusion_matrix
import os
from datetime import datetime


def parse_fasta_with_biotype(fasta_path):
    """
    Parse a FASTA file and extract biotype from the header.
    Args:
        fasta_path (str): Path to FASTA file.
    Returns:
        pd.DataFrame: DataFrame with columns ['id', 'seq', 'biotype', 'header']
    """
    records = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        header = record.description
        # Extract biotype from header
        m = re.search(r"biotype=([^\s]+)", header)
        biotype = m.group(1) if m else "NA"
        records.append({
            'id': record.id,
            'seq': str(record.seq),
            'biotype': biotype,
            'header': header
        })
    return pd.DataFrame(records)

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
    # Normalize
    for kmer in counts:
        counts[kmer] /= total
    return counts

def has_kozak(seq):
    """
    Check for presence of a Kozak consensus sequence around ATG start codons.
    Args:
        seq (str): Nucleotide sequence.
    Returns:
        int: 1 if Kozak found, 0 otherwise.
    """
    # Kozak consensus: GCCRCCATGG (R = A/G), centered on ATG start codon
    # Search for GCC[AG]CCATGG, allowing for some flexibility around the first ATG
    seq = seq.upper()
    for m in re.finditer('ATG', seq):
        start = m.start()
        # Check for Kozak context: 6 bases before and 3 after ATG
        if start >= 6 and start+3+3 <= len(seq):
            context = seq[start-6:start+3+3]
            if re.match(r'GCC[AG]CCATGG', context):
                return 1
    return 0

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

def extract_features(df, k=3):
    """
    Extract features for each sequence in the DataFrame.
    Args:
        df (pd.DataFrame): DataFrame with 'seq' and 'biotype'.
        k (int): k-mer size.
    Returns:
        pd.DataFrame: Feature matrix.
    """
    features = []
    all_kmers = set()
    # First pass to get all kmers
    for seq in df['seq']:
        all_kmers.update(kmer_freq(seq, k).keys())
    all_kmers = sorted(all_kmers)
    for _, row in df.iterrows():
        feats = {}
        feats['length'] = len(row['seq'])
        feats['gc_content'] = gc_content(row['seq'])
        feats['entropy'] = sequence_entropy(row['seq'])
        feats['kozak'] = has_kozak(row['seq'])
        kmer_counts = kmer_freq(row['seq'], k)
        for kmer in all_kmers:
            feats[f'kmer_{kmer}'] = kmer_counts.get(kmer, 0)
        feats['biotype'] = row['biotype']
        features.append(feats)
    return pd.DataFrame(features)

def main():
    """
    Main workflow: summarizes and classifies FASTA sequences by biotype, extracts features, visualizes, and performs ML classification.
    """
    parser = argparse.ArgumentParser(description="Summarize and classify FASTA sequences by biotype.")
    parser.add_argument('-i', '--input', type=str, default='data/canonical_genes.fasta', help='Input FASTA file')
    parser.add_argument('-k', '--kmer', type=int, default=3, help='k-mer size')
    parser.add_argument('-p', '--plotdir', type=str, default='plots', help='Directory to save plots and summary')
    args = parser.parse_args()

    os.makedirs(args.plotdir, exist_ok=True)
    run_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    readme_path = os.path.join(args.plotdir, 'README.txt')
    with open(readme_path, 'w') as f:
        f.write(f"Plots generated from: {args.input}\n")
        f.write(f"Run date: {run_time}\n")

    print(f"Reading {args.input} ...")
    df = parse_fasta_with_biotype(args.input)
    print(f"Loaded {len(df)} sequences.")

    print("Extracting features ...")
    features = extract_features(df, k=args.kmer)
    print(features.describe(include='all'))

    # Visualize characteristics by biotype
    print("\nVisualizing sequence characteristics by biotype ...")
    plt.figure(figsize=(10, 6))
    sns.violinplot(x=df['biotype'], y=features['length'])
    plt.title('Sequence Length by Biotype')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(args.plotdir, 'length_by_biotype.png'))
    plt.close()

    plt.figure(figsize=(10, 6))
    sns.violinplot(x=df['biotype'], y=features['gc_content'])
    plt.title('GC Content by Biotype')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(args.plotdir, 'gc_content_by_biotype.png'))
    plt.close()

    plt.figure(figsize=(10, 6))
    sns.violinplot(x=df['biotype'], y=features['entropy'])
    plt.title('Sequence Entropy by Biotype')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(args.plotdir, 'entropy_by_biotype.png'))
    plt.close()

    plt.figure(figsize=(10, 6))
    sns.barplot(x=df['biotype'], y=features['kozak'], estimator=np.mean)
    plt.title('Fraction with Kozak Sequence by Biotype')
    plt.xticks(rotation=45, ha='right')
    plt.ylabel('Fraction with Kozak')
    plt.tight_layout()
    plt.savefig(os.path.join(args.plotdir, 'kozak_by_biotype.png'))
    plt.close()

    # Visualize top k-mers (mean frequency across all sequences)
    kmer_cols = [col for col in features.columns if col.startswith('kmer_')]
    if kmer_cols:
        mean_kmers = features[kmer_cols].mean().sort_values(ascending=False)[:10]
        plt.figure(figsize=(8, 5))
        sns.barplot(x=mean_kmers.values, y=mean_kmers.index)
        plt.title(f'Top 10 {args.kmer}-mer Frequencies (mean across all sequences)')
        plt.xlabel('Mean Frequency')
        plt.tight_layout()
        plt.savefig(os.path.join(args.plotdir, f'top10_{args.kmer}mer_freq.png'))
        plt.close()

    # ML classification
    print("\nClassifying by biotype ...")
    X = features.drop('biotype', axis=1)
    y = features['biotype']
    # Remove biotypes with <2 samples
    biotype_counts = y.value_counts()
    valid_biotypes = biotype_counts[biotype_counts >= 2].index
    mask = y.isin(valid_biotypes)
    X = X[mask]
    y = y[mask]
    # Only classify if there are at least 2 classes
    if len(set(y)) > 1:
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)
        clf = RandomForestClassifier(n_estimators=100, random_state=42)
        clf.fit(X_train, y_train)
        y_pred = clf.predict(X_test)
        print(classification_report(y_test, y_pred))

        # Confusion matrix
        cm = confusion_matrix(y_test, y_pred, labels=clf.classes_)
        plt.figure(figsize=(8, 6))
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=clf.classes_, yticklabels=clf.classes_)
        plt.xlabel('Predicted')
        plt.ylabel('True')
        plt.title('Confusion Matrix')
        plt.tight_layout()
        plt.savefig(os.path.join(args.plotdir, 'confusion_matrix.png'))
        plt.close()

        # Feature importances
        importances = clf.feature_importances_
        indices = np.argsort(importances)[-20:][::-1]  # Top 20 features
        feature_names = X.columns[indices]
        plt.figure(figsize=(10, 6))
        sns.barplot(x=importances[indices], y=feature_names)
        plt.xlabel('Importance')
        plt.title('Top 20 Feature Importances')
        plt.tight_layout()
        plt.savefig(os.path.join(args.plotdir, 'top20_feature_importances.png'))
        plt.close()
    else:
        print("Not enough biotype classes for classification.")

if __name__ == "__main__":
    main() 