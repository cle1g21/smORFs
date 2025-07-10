"""
postprocess_segment_results.py
-----------------------------
Post-processes the region_analysis/all_transcripts_summary.txt file and generates summary plots for protein coding probabilities by primary set and smORF biotype.
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Load the summary file
df = pd.read_csv('region_analysis/all_transcripts_summary.txt', sep='\t')

# Merge in smORF biotype from smorfs_v2.csv
smorfs = pd.read_csv('data/smorfs_v2.csv', low_memory=False)
# Drop duplicates to avoid multiple matches per transcript
smorfs_biotype = smorfs[['transcript', 'gene_biotype']].drop_duplicates()
df = df.merge(smorfs_biotype, left_on='transcript_id', right_on='transcript', how='left')
df = df.rename(columns={'gene_biotype': 'smorf_biotype'})
df.drop(columns=['transcript'], inplace=True)
print(f"Merged smORF biotype: {df['smorf_biotype'].nunique()} unique biotypes found.")

# Filter for lncRNA biotype and print
print("\n--- smORFs from lncRNAs ---")
lncrna_smorfs = df[df['smorf_biotype'] == 'lncRNA'].copy()
if lncrna_smorfs.empty:
    print("No smORFs from lncRNAs found in the summary file.")
else:
    with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.width', 1000):
        print(lncrna_smorfs[['transcript_id', 'primary_set', 'smorf_biotype', 'protein_coding_mean', 'smorf_protein_coding_mean']])
print("-----------------------------\n")

# Ensure 'primary_set' is string and handle possible NaNs
df['primary_set'] = df['primary_set'].astype(str).fillna('no')

# Only keep rows with valid protein_coding_mean values
df = df[pd.to_numeric(df['protein_coding_mean'], errors='coerce').notnull()]
df['protein_coding_mean'] = df['protein_coding_mean'].astype(float)

# Plot distribution
plt.figure(figsize=(8, 6))
sns.histplot(data=df, x='protein_coding_mean', hue='primary_set', bins=50, kde=True, stat='density', common_norm=False, palette='Set2')
plt.title('Distribution of Protein Coding Probabilities by Primary Set')
plt.xlabel('Protein Coding Probability (Transcript Region Mean)')
plt.ylabel('Density')
plt.legend(title='Primary Set', labels=['No', 'Yes'])
plt.tight_layout()

os.makedirs('plots', exist_ok=True)
plt.savefig('plots/protein_coding_probability_distribution_by_primary_set.png')
plt.close()

print('Plot saved to plots/protein_coding_probability_distribution_by_primary_set.png')

# Filter for plotting smorf data
df_smorf = df[pd.to_numeric(df['smorf_protein_coding_mean'], errors='coerce').notnull()].copy()
df_smorf['smorf_protein_coding_mean'] = df_smorf['smorf_protein_coding_mean'].astype(float)

# Plot distribution for smORF region
plt.figure(figsize=(8, 6))
sns.histplot(data=df_smorf, x='smorf_protein_coding_mean', hue='primary_set', bins=50, kde=True, stat='density', common_norm=False, palette='Set2')
plt.title('Distribution of Protein Coding Probabilities in smORF Region by Primary Set')
plt.xlabel('Protein Coding Probability (smORF Region Mean)')
plt.ylabel('Density')
plt.legend(title='Primary Set', labels=['No', 'Yes'])
plt.tight_layout()

plt.savefig('plots/smorf_protein_coding_probability_distribution_by_primary_set.png')
plt.close()

print('Plot saved to plots/smorf_protein_coding_probability_distribution_by_primary_set.png')

# Plot distribution for smORF region by biotype
if 'smorf_biotype' in df_smorf.columns:
    plt.figure(figsize=(12, 8))
    # Filter out NaNs for plotting
    plot_df = df_smorf.dropna(subset=['smorf_protein_coding_mean', 'smorf_biotype'])
    
    # Get top 10 biotypes by count to avoid overly cluttered plot
    top_biotypes = plot_df['smorf_biotype'].value_counts().nlargest(10).index
    plot_df_top = plot_df[plot_df['smorf_biotype'].isin(top_biotypes)]

    sns.histplot(data=plot_df_top, x='smorf_protein_coding_mean', hue='smorf_biotype', bins=50, kde=True, stat='density', common_norm=False, palette='tab20')
    plt.title('Distribution of Protein Coding Probabilities in smORF Region by Biotype (Top 10)')
    plt.xlabel('Protein Coding Probability (smORF Region Mean)')
    plt.ylabel('Density')
    plt.legend(title='smORF Biotype', loc='upper left')
    plt.tight_layout()
    plt.savefig('plots/smorf_protein_coding_probability_distribution_by_biotype.png')
    plt.close()
    print('Plot saved to plots/smorf_protein_coding_probability_distribution_by_biotype.png') 