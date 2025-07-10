import os
import pandas as pd
import pyfaidx
import gtfparse
import polars as pl
from tqdm import tqdm
import gc

import nucleotide_transformer


if "COLAB_TPU_ADDR" in os.environ:
    from jax.tools import colab_tpu

    colab_tpu.setup_tpu()

from Bio import SeqIO
import gzip
import haiku as hk
import jax
import jax.numpy as jnp
import numpy as np
import seaborn as sns
from typing import List
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from nucleotide_transformer.enformer.pretrained import get_pretrained_segment_enformer_model
from nucleotide_transformer.enformer.features import FEATURES

jax.config.update("jax_platform_name", "cpu")

backend = "cpu"
devices = jax.devices(backend)
num_devices = len(devices)
print(f"Devices found: {devices}")

# seaborn settings
sns.set_style("whitegrid")
sns.set_context(
    "notebook",
    font_scale=1,
    rc={
        "font.size": 14,
        "axes.titlesize": 18,
        "axes.labelsize": 18,
        "xtick.labelsize": 16,
        "ytick.labelsize": 16,
        "legend.fontsize": 16,
        }
)

plt.rcParams['xtick.bottom'] = True
plt.rcParams['ytick.left'] = True

# set colors
colors = sns.color_palette("Set2").as_hex()
colors2 = sns.color_palette("husl").as_hex()

# Rearrange order of the features to match Fig.3 from the paper.
features_rearranged = [
 'protein_coding_gene',
 'lncRNA',
 '5UTR',
 '3UTR',
 'exon',
 'intron',
 'splice_donor',
 'splice_acceptor',
 'promoter_Tissue_specific',
 'promoter_Tissue_invariant',
 'enhancer_Tissue_specific',
 'enhancer_Tissue_invariant',
 'CTCF-bound',
 'polyA_signal',
]

def plot_features(
    predicted_probabilities_all,
    seq_length: int,
    features: List[str],
    order_to_plot: List[str],
    fig_width=8,
    transcript_id: str = "",
    transcript_location=None,
    other_annotations=None,
    chrom: str = "",
    window_start: int = 0,
    window_end: int = 0,
    save_prefix: str = "segment_plot",
    smorf_location=None
):
    """
    Plot predicted probabilities for genomic features and annotate transcript, smORF, and gene context.

    Args:
        predicted_probabilities_all: Probabilities per genomic feature for each nucleotide in the DNA sequence.
        seq_length: DNA sequence length.
        features: List of genomic features to plot.
        order_to_plot: Order in which to plot the genomic features.
        fig_width: Width of the figure.
        transcript_id: ID of the main transcript for the title.
        transcript_location: Tuple (start, end) of the transcript location relative to the sequence.
        other_annotations: List of dicts for other genes to annotate.
        chrom: Chromosome name.
        window_start: Genomic start coordinate of the window.
        window_end: Genomic end coordinate of the window.
        save_prefix: Prefix for saved image files (including directory).
        smorf_location: Tuple (start, end) of the smORF location relative to the sequence.
    """
    # Ensure plots directory exists
    os.makedirs(os.path.dirname(save_prefix), exist_ok=True)

    sc = 1.8
    n_panels = 8  # 7 for probabilities, 1 for gene track

    fig, axes = plt.subplots(
        n_panels, 1,
        figsize=(fig_width * sc, (n_panels + 0.5) * sc),
        sharex=True,
        gridspec_kw={'height_ratios': [2, 2, 2, 2, 2, 2, 2, 1.2]}
    )
    fig.subplots_adjust(hspace=0.1)

    # Plot probabilities on the first 7 axes
    prob_axes = axes[:7]
    
    # Add vertical span for the main transcript's location on all probability plots
    if transcript_location:
        prob_axes[0].axvspan(
            transcript_location[0], transcript_location[1],
            color="red", alpha=0.15, zorder=0, label="Transcript Location"
        )
        for ax in prob_axes[1:]:
            ax.axvspan(
                transcript_location[0], transcript_location[1],
                color="red", alpha=0.15, zorder=0
            )
    # Add vertical span for the smORF location on all probability plots
    if smorf_location:
        prob_axes[0].axvspan(
            smorf_location[0], smorf_location[1],
            color="green", alpha=0.18, zorder=0, label="smORF Location"
        )
        for ax in prob_axes[1:]:
            ax.axvspan(
                smorf_location[0], smorf_location[1],
                color="green", alpha=0.18, zorder=0
            )

    for n, feat in enumerate(order_to_plot):
        feat_id = features.index(feat)
        prob_dist = predicted_probabilities_all[:, feat_id]
        ax = prob_axes[n // 2]
        try:
            id_color = colors[feat_id]
        except IndexError:
            id_color = colors2[feat_id - 8]
        ax.plot(
            prob_dist,
            color=id_color,
            label=feat,
            linestyle="-",
            linewidth=1.5,
        )
        ax.grid(False)
        ax.spines['bottom'].set_color('black')
        ax.spines['top'].set_color('black')
        ax.spines['right'].set_color('black')
        ax.spines['left'].set_color('black')
        ax.set_ylim(0, 1.05)
        ax.set_ylabel("Prob.")
        ax.legend(loc="upper left", bbox_to_anchor=(1, 1), borderaxespad=0)

    # --- Gene Track Annotation ---
    gene_ax = axes[7]
    gene_ax.set_yticks([])
    gene_ax.spines['left'].set_visible(False)
    gene_ax.spines['right'].set_visible(False)
    gene_ax.spines['top'].set_visible(False)

    y_positions = [0.6, 0.3, 0.0]  # For staggering gene labels

    # Draw main transcript
    if transcript_location:
        start, end = transcript_location
        rect = patches.Rectangle((start, 0.6), end - start, 0.25, facecolor='red', edgecolor='black', clip_on=False)
        gene_ax.add_patch(rect)
        gene_ax.text(start + (end - start) / 2, 0.9, transcript_id, ha='center', va='bottom', fontsize=9, fontweight='bold')
    # Draw smORF
    if smorf_location:
        s_start, s_end = smorf_location
        rect = patches.Rectangle((s_start, 0.6), s_end - s_start, 0.25, facecolor='green', edgecolor='black', clip_on=False)
        gene_ax.add_patch(rect)
        gene_ax.text(s_start + (s_end - s_start) / 2, 0.5, 'smORF', ha='center', va='top', fontsize=8, fontweight='bold', color='green')

    # Draw other genes
    if other_annotations:
        for i, gene in enumerate(other_annotations):
            start, end = gene['start'], gene['end']
            y_pos = y_positions[i % len(y_positions)]
            rect = patches.Rectangle((start, y_pos), end - start, 0.2, facecolor='skyblue', edgecolor='black', clip_on=False)
            gene_ax.add_patch(rect)
            arrow = " ->" if gene['strand'] == '+' else "<- "
            gene_ax.text(start + (end - start) / 2, y_pos - 0.05, f"{gene['name']}{arrow}", ha='center', va='top', fontsize=8)

    # --- Axis labels and title ---
    gene_ax.set_xlabel(f"Genomic Position on {chrom} ({window_start:,} - {window_end:,})")
    axes[0].set_title(f"Probabilities for features around {transcript_id}", fontweight="bold")

    fig.tight_layout()
    fig.savefig(f"{save_prefix}_{transcript_id}.png", dpi=200)
    plt.close(fig)  # Close main figure

    # --- Save zoomed-in plot on transcript of interest ---
    if transcript_location and chrom and window_start and window_end:
        zoom_margin = 2000  # 2kb margin on each side
        zoom_start = max(0, transcript_location[0] - zoom_margin)
        zoom_end = min(seq_length, transcript_location[1] + zoom_margin)

        fig_zoom, axes_zoom = plt.subplots(
            n_panels, 1,
            figsize=(fig_width * sc, (n_panels + 0.5) * sc),
            sharex=True,
            gridspec_kw={'height_ratios': [2, 2, 2, 2, 2, 2, 2, 1.2]}
        )
        fig_zoom.subplots_adjust(hspace=0.1)

        # Probability plots
        prob_axes_zoom = axes_zoom[:7]
        if transcript_location:
            prob_axes_zoom[0].axvspan(
                transcript_location[0], transcript_location[1],
                color="red", alpha=0.15, zorder=0, label="Transcript Location"
            )
            for ax in prob_axes_zoom[1:]:
                ax.axvspan(
                    transcript_location[0], transcript_location[1],
                    color="red", alpha=0.15, zorder=0
                )
        if smorf_location:
            prob_axes_zoom[0].axvspan(
                smorf_location[0], smorf_location[1],
                color="green", alpha=0.18, zorder=0, label="smORF Location"
            )
            for ax in prob_axes_zoom[1:]:
                ax.axvspan(
                    smorf_location[0], smorf_location[1],
                    color="green", alpha=0.18, zorder=0
                )
        for n, feat in enumerate(order_to_plot):
            feat_id = features.index(feat)
            prob_dist = predicted_probabilities_all[:, feat_id]
            ax = prob_axes_zoom[n // 2]
            try:
                id_color = colors[feat_id]
            except IndexError:
                id_color = colors2[feat_id - 8]
            ax.plot(
                np.arange(zoom_start, zoom_end),
                prob_dist[zoom_start:zoom_end],
                color=id_color,
                label=feat,
                linestyle="-",
                linewidth=1.5,
            )
            ax.grid(False)
            ax.spines['bottom'].set_color('black')
            ax.spines['top'].set_color('black')
            ax.spines['right'].set_color('black')
            ax.spines['left'].set_color('black')
            ax.set_ylim(0, 1.05)
            ax.set_ylabel("Prob.")
            ax.legend(loc="upper left", bbox_to_anchor=(1, 1), borderaxespad=0)

        # Gene track
        gene_ax_zoom = axes_zoom[7]
        gene_ax_zoom.set_yticks([])
        gene_ax_zoom.spines['left'].set_visible(False)
        gene_ax_zoom.spines['right'].set_visible(False)
        gene_ax_zoom.spines['top'].set_visible(False)
        y_positions = [0.6, 0.3, 0.0]
        # Main transcript
        if transcript_location:
            start, end = transcript_location
            rect = patches.Rectangle((start, 0.6), end - start, 0.25, facecolor='red', edgecolor='black', clip_on=False)
            gene_ax_zoom.add_patch(rect)
            gene_ax_zoom.text(start + (end - start) / 2, 0.9, transcript_id, ha='center', va='bottom', fontsize=9, fontweight='bold')
        # smORF
        if smorf_location:
            s_start, s_end = smorf_location
            rect = patches.Rectangle((s_start, 0.6), s_end - s_start, 0.25, facecolor='green', edgecolor='black', clip_on=False)
            gene_ax_zoom.add_patch(rect)
            gene_ax_zoom.text(s_start + (s_end - s_start) / 2, 0.5, 'smORF', ha='center', va='top', fontsize=8, fontweight='bold', color='green')
        # Other genes
        if other_annotations:
            for i, gene in enumerate(other_annotations):
                start, end = gene['start'], gene['end']
                y_pos = y_positions[i % len(y_positions)]
                rect = patches.Rectangle((start, y_pos), end - start, 0.2, facecolor='skyblue', edgecolor='black', clip_on=False)
                gene_ax_zoom.add_patch(rect)
                arrow = " ->" if gene['strand'] == '+' else "<- "
                gene_ax_zoom.text(start + (end - start) / 2, y_pos - 0.05, f"{gene['name']}{arrow}", ha='center', va='top', fontsize=8)
        # X label
        genomic_zoom_start = window_start + zoom_start
        genomic_zoom_end = window_start + zoom_end - 1
        gene_ax_zoom.set_xlabel(f"Genomic Position on {chrom} ({genomic_zoom_start:,} - {genomic_zoom_end:,})")
        axes_zoom[0].set_title(f"Zoomed: Probabilities for features around {transcript_id}", fontweight="bold")
        # Set xlim for all axes
        for ax in axes_zoom:
            ax.set_xlim(zoom_start, zoom_end)
        # Save
        fig_zoom.tight_layout()
        fig_zoom.savefig(f"{save_prefix}_{transcript_id}_zoom.png", dpi=200)
        plt.close(fig_zoom)  # Close zoomed-in figure

def predict_gene_model_from_probabilities(
    probabilities,
    region_start,
    region_end,
    window_start,
    chrom,
    transcript_id,
    strand='+',
    threshold=0.75,
    output_dir="predicted_gtf"
):
    """
    Predict gene model from model probabilities and save as GTF (IGV compatible).
    Only applies to the region [region_start, region_end) in the window.
    Merges adjacent/overlapping blocks of the same feature.

    Args:
        probabilities: Model output probabilities (array-like).
        region_start: Start index of region in window.
        region_end: End index of region in window.
        window_start: Genomic start coordinate of window.
        chrom: Chromosome name.
        transcript_id: Transcript ID for annotation.
        strand: Genomic strand ('+' or '-').
        threshold: Probability threshold for calling protein-coding.
        output_dir: Directory to save GTF file.
    """
    os.makedirs(output_dir, exist_ok=True)
    # Get feature indices
    idx_protein = FEATURES.index('protein_coding_gene')
    idx_exon = FEATURES.index('exon')
    idx_intron = FEATURES.index('intron')
    idx_5utr = FEATURES.index('5UTR')
    idx_3utr = FEATURES.index('3UTR')

    # Restrict to region
    region_probs = probabilities[region_start:region_end, :]
    region_offset = window_start + region_start

    # Find high-confidence protein-coding regions
    protein_track = region_probs[:, idx_protein]
    is_protein = protein_track > threshold

    # Find contiguous blocks of protein-coding
    from itertools import groupby
    from operator import itemgetter
    blocks = []
    for k, g in groupby(enumerate(is_protein), key=itemgetter(1)):
        if k:  # True = protein-coding
            group = list(g)
            start = group[0][0]
            end = group[-1][0]
            blocks.append((start, end))

    # For each block, assign features by argmax among exon/intron/5UTR/3UTR
    feature_indices = [idx_exon, idx_intron, idx_5utr, idx_3utr]
    feature_names = ['exon', 'intron', '5UTR', '3UTR']
    feature_blocks = []  # (start, end, feature_name)
    gene_id = f"pred_{transcript_id}"
    transcript_id_str = f"pred_{transcript_id}.1"
    for block_start, block_end in blocks:
        block_probs = region_probs[block_start:block_end+1, :]
        feature_track = np.argmax(block_probs[:, feature_indices], axis=1)
        # Collapse contiguous features
        last_feat = None
        feat_start = block_start
        for i, feat_idx in enumerate(feature_track):
            abs_pos = block_start + i
            if last_feat is None:
                last_feat = feat_idx
                feat_start = abs_pos
            elif feat_idx != last_feat:
                feat_name = feature_names[last_feat]
                feature_blocks.append((feat_start, abs_pos-1, feat_name))
                last_feat = feat_idx
                feat_start = abs_pos
        # Last feature in block
        feat_name = feature_names[last_feat]
        feature_blocks.append((feat_start, block_end, feat_name))

    # Merge adjacent/overlapping blocks of the same feature
    merged_blocks = []
    if feature_blocks:
        feature_blocks.sort(key=lambda x: (x[2], x[0]))  # sort by feature, then start
        from collections import defaultdict
        by_feature = defaultdict(list)
        for start, end, feat in feature_blocks:
            by_feature[feat].append((start, end))
        for feat, blocks in by_feature.items():
            blocks.sort()
            merged = []
            for s, e in blocks:
                if not merged:
                    merged.append([s, e])
                else:
                    last = merged[-1]
                    if s <= last[1] + 1:  # adjacent or overlapping
                        last[1] = max(last[1], e)
                    else:
                        merged.append([s, e])
            for s, e in merged:
                merged_blocks.append((s, e, feat))
    merged_blocks.sort(key=lambda x: x[0])  # sort by start

    # Write merged blocks to GTF
    gtf_lines = []
    for s, e, feat_name in merged_blocks:
        gtf_lines.append(
            f"{chrom}\tSegmentNT\t{feat_name}\t{region_offset+s+1}\t{region_offset+e+1}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id_str}\";"
        )
    # Add gene and transcript features for the union of all blocks
    if merged_blocks:
        gene_start = min(s for s, _, _ in merged_blocks)
        gene_end = max(e for _, e, _ in merged_blocks)
        gtf_lines.append(
            f"{chrom}\tSegmentNT\tgene\t{region_offset+gene_start+1}\t{region_offset+gene_end+1}\t.\t{strand}\t.\tgene_id \"{gene_id}\";"
        )
        gtf_lines.append(
            f"{chrom}\tSegmentNT\ttranscript\t{region_offset+gene_start+1}\t{region_offset+gene_end+1}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id_str}\";"
        )
    gtf_path = os.path.join(output_dir, f"{transcript_id}_predicted.gtf")
    with open(gtf_path, 'w') as gtf:
        for line in gtf_lines:
            gtf.write(line + '\n')
    print(f"Predicted gene model saved to {gtf_path}")

# --- Main workflow ---
"""
Main workflow:
- Loads smORF table and GTF/genome reference files
- For each transcript, extracts a fixed-size window, runs the model, and analyses the region
- Generates plots, region analysis files, summary, and predicted GTFs
"""
# --- Step 1: Get transcript coordinates from GTF ---
# Load smORF table
smorfs_df = pd.read_csv('data/smorfs_v2.csv')
all_transcript_ids = smorfs_df['transcript'].unique()

# Load GTF and genome FASTA ONCE
GTF_PATH = "data/Homo_sapiens.GRCh38.114.gtf.gz"
gtf = gtfparse.read_gtf(GTF_PATH)
MODEL_INPUT_LENGTH = 196_608
GENOME_FASTA = "data/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
genome = pyfaidx.Fasta(GENOME_FASTA)

# Prepare summary output
os.makedirs('region_analysis', exist_ok=True)
summary_header = 'transcript_id\tprimary_set\tregion_start_in_window\tregion_end_in_window\tprotein_coding_mean\tprotein_coding_min\tprotein_coding_max\tlncrna_mean\tlncrna_min\tlncrna_max\tsmorf_region_start_in_window\tsmorf_region_end_in_window\tsmorf_protein_coding_mean\tsmorf_protein_coding_min\tsmorf_protein_coding_max\tsmorf_lncrna_mean\tsmorf_lncrna_min\tsmorf_lncrna_max\n'

# Write header once at the start
with open('region_analysis/all_transcripts_summary.txt', 'w') as f:
    f.write(summary_header)

# Load model and tokenizer ONCE
parameters, state, forward_fn, tokenizer, config = get_pretrained_segment_enformer_model()
forward_fn = hk.transform_with_state(forward_fn)
apply_fn = jax.pmap(forward_fn.apply, devices=devices)
random_key = jax.random.PRNGKey(seed=0)
keys = jax.device_put_replicated(random_key, devices=devices)
parameters = jax.device_put_replicated(parameters, devices=devices)
state = jax.device_put_replicated(state, devices=devices)

for TRANSCRIPT_ID in tqdm(all_transcript_ids, desc='Transcripts'):
    print(f"Processing {TRANSCRIPT_ID}")
    # Get primary_set info
    primary_set = ''
    smorf_row = smorfs_df[smorfs_df['transcript'] == TRANSCRIPT_ID]
    if not smorf_row.empty and 'primary_set' in smorf_row.columns:
        primary_set = str(smorf_row.iloc[0]['primary_set'])

    # GTF extraction
    transcript_df = gtf.filter((gtf['feature'] == 'transcript') & (gtf['transcript_id'] == TRANSCRIPT_ID))
    if transcript_df.is_empty():
        print(f"Transcript {TRANSCRIPT_ID} not found in GTF. Skipping.")
        continue
    chrom = transcript_df['seqname'][0]
    start = int(transcript_df['start'][0])
    end = int(transcript_df['end'][0])
    strand = transcript_df['strand'][0]
    main_gene_name = transcript_df['gene_name'][0]

    # --- Step 2: Extract a fixed-size window centered on the transcript ---
    # Calculate transcript center
    transcript_center = start + (end - start) // 2
    # Define extraction window
    window_start = transcript_center - (MODEL_INPUT_LENGTH // 2)
    window_start = max(1, window_start)
    window_end = window_start + MODEL_INPUT_LENGTH

    # Calculate transcript's position relative to the window for plotting
    if strand == '+':
        relative_start = start - window_start
        relative_end = end - window_start
    else:  # strand == '-'
        relative_start = MODEL_INPUT_LENGTH - (end - window_start)
        relative_end = MODEL_INPUT_LENGTH - (start - window_start)
    relative_start = max(0, relative_start)
    relative_end = min(MODEL_INPUT_LENGTH, relative_end)
    transcript_location_in_window = (relative_start, relative_end)

    # pyfaidx uses 0-based, end-exclusive indexing for slices.
    sequence = genome[chrom][window_start - 1 : window_start - 1 + MODEL_INPUT_LENGTH].seq
    if strand == '-':
        sequence = pyfaidx.Sequence(sequence).reverse.complement.seq
    if len(sequence) < MODEL_INPUT_LENGTH:
        padding_needed = MODEL_INPUT_LENGTH - len(sequence)
        sequence = 'N' * padding_needed + sequence
    sequences = [sequence]

    tokens_ids = [b[1] for b in tokenizer.batch_tokenize(sequences)]
    tokens = jnp.stack([jnp.asarray(tokens_ids, dtype=jnp.int32)]*num_devices, axis=0)
    outs, state = apply_fn(parameters, state, keys, tokens)
    logits = outs["logits"]
    probabilities = np.asarray(jax.nn.softmax(logits, axis=-1))[...,-1]

    # --- Step 3: Find other genes in the window for annotation ---
    other_annotations = []
    genes_in_region = gtf.filter(
        (pl.col('feature') == 'gene') &
        (pl.col('seqname') == chrom) &
        (pl.col('start') < window_end) & 
        (pl.col('end') > window_start)
    )
    if not genes_in_region.is_empty():
        genes_list = genes_in_region.to_dicts()
        for gene in genes_list:
            if gene['gene_name'] == main_gene_name:
                continue
            g_start, g_end, g_strand = gene['start'], gene['end'], gene['strand']
            if strand == '+':
                rel_g_start = g_start - window_start
                rel_g_end = g_end - window_start
            else:
                rel_g_start = MODEL_INPUT_LENGTH - (g_end - window_start)
                rel_g_end = MODEL_INPUT_LENGTH - (g_start - window_start)
            rel_g_start = max(0, rel_g_start)
            rel_g_end = min(MODEL_INPUT_LENGTH, rel_g_end)
            if rel_g_end > rel_g_start:
                other_annotations.append({
                    'name': gene['gene_name'],
                    'start': rel_g_start,
                    'end': rel_g_end,
                    'strand': g_strand
                })

    # --- smORF region analysis from smorfs_v2.csv ---
    smorf_location_in_window = None
    smorf_summary = {
        'smorf_region_start_in_window': '',
        'smorf_region_end_in_window': '',
        'smorf_protein_coding_mean': '',
        'smorf_protein_coding_min': '',
        'smorf_protein_coding_max': '',
        'smorf_lncrna_mean': '',
        'smorf_lncrna_min': '',
        'smorf_lncrna_max': ''
    }
    if not smorf_row.empty:
        coord_str = smorf_row.iloc[0]['genomic_coordinates (1-\nbased)'] if 'genomic_coordinates (1-\nbased)' in smorf_row.columns else smorf_row.iloc[0]['genomic_coordinates (1-based)']
        import re
        match = re.match(r'chr?(\w+):(\d+)-(\d+)', str(coord_str).replace(",", ""))
        if match:
            smorf_chrom, smorf_start_1b, smorf_end_1b = match.groups()
            smorf_chrom = smorf_chrom.replace('chr', '')
            smorf_start_1b = int(smorf_start_1b)
            smorf_end_1b = int(smorf_end_1b)
            if smorf_chrom == chrom.replace('chr', ''):
                smorf_start_in_window = smorf_start_1b - window_start
                smorf_end_in_window = smorf_end_1b - window_start
                smorf_start_in_window = max(0, smorf_start_in_window)
                smorf_end_in_window = min(MODEL_INPUT_LENGTH, smorf_end_in_window)
                if smorf_end_in_window > smorf_start_in_window:
                    smorf_location_in_window = (smorf_start_in_window, smorf_end_in_window)
                    protein_coding_idx = FEATURES.index('protein_coding_gene')
                    lncrna_idx = FEATURES.index('lncRNA')
                    smorf_region_probs = probabilities[0, 0, smorf_start_in_window:smorf_end_in_window, :]
                    smorf_protein_coding_probs = smorf_region_probs[:, protein_coding_idx]
                    smorf_lncrna_probs = smorf_region_probs[:, lncrna_idx]
                    smorf_summary = {
                        'smorf_region_start_in_window': smorf_start_in_window,
                        'smorf_region_end_in_window': smorf_end_in_window,
                        'smorf_protein_coding_mean': f"{np.mean(smorf_protein_coding_probs):.6f}",
                        'smorf_protein_coding_min': f"{np.min(smorf_protein_coding_probs):.6f}",
                        'smorf_protein_coding_max': f"{np.max(smorf_protein_coding_probs):.6f}",
                        'smorf_lncrna_mean': f"{np.mean(smorf_lncrna_probs):.6f}",
                        'smorf_lncrna_min': f"{np.min(smorf_lncrna_probs):.6f}",
                        'smorf_lncrna_max': f"{np.max(smorf_lncrna_probs):.6f}"
                    }

    # --- Save region analysis for smORF region ---
    protein_coding_idx = FEATURES.index('protein_coding_gene')
    lncrna_idx = FEATURES.index('lncRNA')
    smorf_start, smorf_end = transcript_location_in_window
    region_probs = probabilities[0, 0, smorf_start:smorf_end, :]
    protein_coding_probs = region_probs[:, protein_coding_idx]
    lncrna_probs = region_probs[:, lncrna_idx]
    summary = {
        'transcript_id': TRANSCRIPT_ID,
        'primary_set': primary_set,
        'region_start_in_window': smorf_start,
        'region_end_in_window': smorf_end,
        'protein_coding_mean': float(np.mean(protein_coding_probs)),
        'protein_coding_min': float(np.min(protein_coding_probs)),
        'protein_coding_max': float(np.max(protein_coding_probs)),
        'lncrna_mean': float(np.mean(lncrna_probs)),
        'lncrna_min': float(np.min(lncrna_probs)),
        'lncrna_max': float(np.max(lncrna_probs)),
    }
    summary_line = f"{summary['transcript_id']}\t{summary['primary_set']}\t{summary['region_start_in_window']}\t{summary['region_end_in_window']}\t{summary['protein_coding_mean']:.6f}\t{summary['protein_coding_min']:.6f}\t{summary['protein_coding_max']:.6f}\t{summary['lncrna_mean']:.6f}\t{summary['lncrna_min']:.6f}\t{summary['lncrna_max']:.6f}\t{smorf_summary['smorf_region_start_in_window']}\t{smorf_summary['smorf_region_end_in_window']}\t{smorf_summary['smorf_protein_coding_mean']}\t{smorf_summary['smorf_protein_coding_min']}\t{smorf_summary['smorf_protein_coding_max']}\t{smorf_summary['smorf_lncrna_mean']}\t{smorf_summary['smorf_lncrna_min']}\t{smorf_summary['smorf_lncrna_max']}\n"

    # Write per-transcript region analysis file
    with open(f'region_analysis/{TRANSCRIPT_ID}_region_probabilities.txt', 'w') as f:
        f.write(summary_header)
        f.write(summary_line)
    # Write summary line incrementally
    with open('region_analysis/all_transcripts_summary.txt', 'a') as f:
        f.write(summary_line)

    # --- Plotting and GTF ---
    plot_features(
        probabilities[0, 0],
        probabilities.shape[-2],
        fig_width=20,
        features=FEATURES,
        order_to_plot=features_rearranged,
        transcript_id=TRANSCRIPT_ID,
        transcript_location=transcript_location_in_window,
        other_annotations=other_annotations,
        chrom=chrom,
        window_start=window_start,
        window_end=(window_start + MODEL_INPUT_LENGTH - 1),
        save_prefix=f"plots/segment_plot",
        smorf_location=smorf_location_in_window
    )
    region_start = max(0, relative_start - 2000)
    region_end = min(MODEL_INPUT_LENGTH, relative_end + 2000)
    predict_gene_model_from_probabilities(
        probabilities[0, 0],
        region_start,
        region_end,
        window_start,
        chrom,
        TRANSCRIPT_ID,
        strand=strand,
        threshold=0.75,
        output_dir="predicted_gtf"
    )
    # Explicitly free memory after all uses of probabilities
    del sequence, probabilities, outs, region_probs, protein_coding_probs, lncrna_probs
    gc.collect()
