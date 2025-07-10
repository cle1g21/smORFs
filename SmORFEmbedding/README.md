# SmORFEmbedding

A toolkit for extracting, analyzing, and embedding transcript and smORF sequences from Ensembl/GENCODE GTF and FASTA files, with downstream statistical and machine learning analysis.

## Environment Setup

Install dependencies with conda:

```bash
conda env create -f environment.yml
conda activate canonical-fasta
```

## Offline Setup (Pre-download All Models and Reference Files)

To run the pipeline fully offline (e.g. on an HPC without internet access), you can pre-download all required Hugging Face models/tokenizers and Ensembl reference files using the provided script:

```bash
python setup_offline.py
```

**What does setup_offline.py download and why?**

| File/Model                                      | Type         | Why is it needed?                                                                                 |
|-------------------------------------------------|--------------|---------------------------------------------------------------------------------------------------|
| `zhihan1996/DNABERT-2-117M`                     | Model        | DNABERT model for DNA sequence embedding (6-mer, 512 tokens). Used in `embed_fasta_transformer.py`|
| `zhihan1996/DNABERT-6`                          | Model        | DNABERT tokenizer and variant. Required for compatibility and tokenization.                       |
| `InstaDeepAI/nucleotide-transformer-500m-human-ref` | Model    | Nucleotide Transformer (500M params, 6-mer, up to 1000 tokens). Used for high-accuracy DNA embedding. |
| `InstaDeepAI/nucleotide-transformer-2.5b-multi-species` | Model | Larger, multi-species Nucleotide Transformer. Optional for advanced users or cross-species work.  |
| `facebook/esm2_t6_8M_UR50D`                     | Model        | ESM2 model for protein/nucleotide embedding (1-mer, 1022 tokens). Used for alternative embedding. |
| `LongSafari/hyenadna-small-32k-seqlen-hf`       | Model        | HyenaDNA model for very long DNA sequences (up to 32,768 nt). Used for full-length or large regions. |
| `data/Homo_sapiens.GRCh38.114.gtf.gz`           | Ensembl GTF  | Gene annotation file (Ensembl/GENCODE). Needed for transcript/gene structure, biotypes, etc.      |
| `data/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz` | Ensembl FASTA | Human reference genome sequence. Used for extracting transcript/genomic sequences.                |
| `data/Homo_sapiens.GRCh38.cdna.all.fa.gz`       | Ensembl FASTA | All cDNA transcript sequences. Used for generating transcriptome FASTA files.                     |
| `data/Homo_sapiens.GRCh38.ncrna.fa.gz`          | Ensembl FASTA | All ncRNA transcript sequences. Used for generating transcriptome FASTA files.                    |

- **Hugging Face models** are required for embedding DNA or protein sequences using transformer models. Each model supports different sequence lengths and biological use cases.
- **Ensembl files** provide the reference genome and gene/transcript annotations needed for extracting, summarising, and embedding transcript sequences.

**Note:** `setup_offline.py` uses `from transformers.utils.hub import snapshot_download` and requires `transformers>=4.29.0`.
If you get an `ImportError` for `snapshot_download`, check your transformers version with:
```bash
python -c "import transformers; print(transformers.__version__)"
```
If needed, upgrade with:
```bash
pip install --upgrade 'transformers>=4.29.0'
```

This will:
- Download and cache all referenced Hugging Face models and tokenizers (DNABERT, Nucleotide Transformer, ESM2, HyenaDNA)
- Download all referenced Ensembl GTF and FASTA files into your `data/` directory

You only need to run this once (on a machine with internet access). Afterward, you can copy the environment and data to an offline system and run all analysis scripts without internet access.

## Directory Structure

- `data/` — Downloaded GTF/FASTA files, generated FASTA, and summary tables
- `plots/` — All generated plots
- `embeddings/` — All generated sequence embeddings

## Main Scripts

### 1. `generate_canonical_fasta.py`
- **Purpose:** Download latest Ensembl GTF, cDNA, and ncRNA files. Create a FASTA with one transcript per gene: canonical for protein-coding, longest for non-coding.
- **Output:** `data/canonical_genes.fasta`, `data/biotype_summary.txt`
- **Usage:**
  ```bash
  python generate_canonical_fasta.py
  ```
  Optional: `--gtf`, `--cdna`, `--ncrna`, `-o` for custom files/output.

### 2. `generate_all_transcripts_fasta.py`
- **Purpose:** Download latest Ensembl GTF, cDNA, and ncRNA files. Create a FASTA with all transcript sequences.
- **Output:** `data/all_transcripts.fasta`
- **Usage:**
  ```bash
  python generate_all_transcripts_fasta.py
  ```
  Optional: `--gtf`, `--cdna`, `--ncrna`, `-o` for custom files/output.

### 3. `summarize_and_classify_fasta.py`
- **Purpose:** Compute sequence statistics (length, GC content, entropy, Kozak, k-mer frequencies) for each transcript in a FASTA. Visualize distributions, classify by biotype, and plot feature importances.
- **Output:** Plots and summary stats in `plots/`.
- **Usage:**
  ```bash
  python summarize_and_classify_fasta.py -i data/canonical_genes.fasta
  ```
  Optional: `-k` for k-mer size, `-p` for plot directory.

### 4. `analyze_smorfs.py`
- **Purpose:** Analyze smORFs from `data/smorfs_v2.csv`, compute statistics for both ORF and full transcript (using a transcriptome FASTA), compare primary vs. non-primary sets, and plot/test differences.
- **Output:** Plots and stats in `plots/smorfs/`, missing transcript log.
- **Usage:**
  ```bash
  python analyze_smorfs.py
  ```
  Optional: `-i` for input CSV, `-t` for transcriptome FASTA, `-o` for output dir, `-k` for k-mer size.

### 5. `embed_fasta_transformer.py`
- **Purpose:** Embed all sequences in a FASTA using DNABERT, Nucleotide Transformer, ESM-2, or HyenaDNA. Computes mean and max pooled embeddings for each sequence, saves to HDF5, and plots a UMAP colored by biotype.
- **Output:** HDF5 file and UMAP plot in `embeddings/`.
- **Usage:**
  ```bash
  # DNABERT
  python embed_fasta_transformer.py --type dnabert --input data/all_transcripts.fasta

  # Nucleotide Transformer
  python embed_fasta_transformer.py --type nucleotide-transformer --input data/all_transcripts.fasta

  # ESM-2 (for protein, use translated sequences)
  python embed_fasta_transformer.py --type esm2 --input data/all_transcripts.fasta

  # HyenaDNA
  python embed_fasta_transformer.py --type hyenadna --input data/all_transcripts.fasta
  ```
  Optional: `--test` for first 100 seqs, `--model` for model checkpoint, `--k` for k-mer size, `--max_len` for window size.

#### Detailed Model and Embedding Strategy

| Model Type              | Default Model Name                                 | k-mer Size | Max Window Length | Tokenization Strategy         |
|------------------------|----------------------------------------------------|------------|-------------------|------------------------------|
| `dnabert`              | zhihan1996/DNABERT-2-117M                          | 6          | 512               | Overlapping 6-mers, joined by space |
| `nucleotide-transformer`| InstaDeepAI/nucleotide-transformer-500m-human-ref  | 6          | 512 (default), **up to 1000** | Overlapping 6-mers, joined by space |
| `esm2`                 | facebook/esm2_t6_8M_UR50D                         | 1          | 1022              | Single nucleotides, joined as string |
| `hyenadna`             | LongSafari/hyenadna-small-32k-seqlen-hf            | 1          | 32768             | Single nucleotides, joined as string |

- **k-mer size**: For DNABERT and Nucleotide Transformer, sequences are split into overlapping k-mers (default k=6). For ESM2 and HyenaDNA, each nucleotide is treated as a token (k=1).
- **Windowing**: Sequences longer than the model's max window length are split into non-overlapping windows of up to `max_len` tokens. Each window is embedded separately.
- **Tokenization**:
  - DNABERT/Nucleotide Transformer: Each window is a space-separated string of k-mers (e.g., "ATGCGA TGCGA..."), as required by the model.
  - ESM2/HyenaDNA: Each window is a contiguous string of nucleotides (e.g., "ATGCGA...").

#### Embedding Combination Strategy

- For each window, the model produces a sequence of hidden states (one per token).
- For each window, the mean of the hidden states is computed to yield a single window embedding.
- For each full sequence (transcript):
  - All window embeddings are concatenated.
  - **Mean pooling**: The mean across all window embeddings is computed (per-dimension) to yield the sequence's mean embedding.
  - **Max pooling**: The max across all window embeddings is computed (per-dimension) to yield the sequence's max embedding.
- Both mean and max pooled embeddings are saved for each sequence.
- All window embeddings (with window start/end positions and biotype) are also saved to a separate HDF5 file for downstream analysis (e.g., UMAP of all windows).

#### Output Files

- `embeddings/<basename>_<model>_embeddings.h5`: Contains mean and max pooled embeddings for each sequence, plus transcript IDs and biotypes.
- `embeddings/<basename>_<model>_all_windows.h5`: Contains all window embeddings, window IDs, biotypes, and window start/end positions.
- `embeddings/<basename>_<model>_umap.png`: UMAP plot of mean sequence embeddings, coloured by biotype.
- `embeddings/<basename>_<model>_all_windows_umap.png`: UMAP plot of all window embeddings, coloured by biotype.

#### Model Selection and Parameter Tuning

- **DNABERT**: Good for short/medium DNA sequences, uses 6-mers. Max window: 512 tokens (i.e., 517 nt for k=6).
- **Nucleotide Transformer**: Similar to DNABERT but larger, also uses 6-mers. Max window: 512 tokens (default), **up to 1000** (Hugging Face model card).
- **ESM2**: Designed for protein, but can be used for nucleotide sequences (k=1). Max window: 1022 tokens.
- **HyenaDNA**: Handles very long sequences (up to 32,768 nt). Use for full-length transcripts or large regions.
- You can override k-mer size (`--k`) and window length (`--max_len`) for advanced use, but defaults are recommended for most cases.
- For large-scale embedding, use a GPU for best performance.
- For offline use, pre-download all models with `python setup_offline.py`.

#### Parameter Effects on Compute Time, Memory, and Embedding Quality

| Parameter         | Description                                                      | Effect on Compute Time/Memory         | Effect on Embedding Quality/Accuracy         |
|------------------|------------------------------------------------------------------|---------------------------------------|---------------------------------------------|
| `--type`         | Model type: dnabert, nucleotide-transformer, esm2, hyenadna      | Larger models (nucleotide-transformer, hyenadna) require more memory and compute; hyenadna can process much longer sequences in one go | Choice of model affects biological relevance and accuracy; larger models may capture more complex features |
| `--model`        | Hugging Face model checkpoint (overrides default for type)        | Some checkpoints may be larger/slower | May improve accuracy if using a more specialised or larger model checkpoint |
| `--k`            | k-mer size (for dnabert, nucleotide-transformer)                  | Larger k increases token count, may increase compute time slightly | Affects the granularity of sequence features captured; default (6) is recommended for most cases |
| `--max_len`      | Max tokens per window                                            | Smaller values increase number of windows (slower, more memory); larger values reduce window count (faster, but may hit model limits). For Nucleotide Transformer, can be set up to 1000. | Too small: may lose context; too large: may cause model errors or memory issues. Use model default unless needed. For Nucleotide Transformer, up to 1000 is supported. |
| `--input`        | Input FASTA file                                                 | Larger files take longer to process   | No effect on accuracy, but more data may improve downstream analyses |
| `--outdir`       | Output directory for embeddings                                  | No effect                            | No effect                                    |
| `--test`         | Only embed first 100 sequences                                   | Much faster, for debugging/testing    | Not representative of full dataset           |

##### Practical Recommendations

- **For large datasets:**
  - Use a GPU for best performance.
  - Stick to default `--max_len` and `--k` unless you have a specific reason to change them.
  - Use `--test` to quickly check pipeline functionality before running on the full dataset.
- **For long sequences:**
  - Use `hyenadna` if you need to embed very long transcripts or genomic regions (up to 32,768 nt per window).
  - For DNABERT/Nucleotide Transformer, sequences longer than 512 tokens will be split into multiple windows, increasing compute time.
- **For best biological accuracy:**
  - Use the default model for each type unless you have a reason to use a different checkpoint.
  - Larger models (e.g., Nucleotide Transformer, HyenaDNA) may capture more complex sequence features but require more resources.
- **For quick tests or debugging:**
  - Use `--test` to process only 100 sequences.
  - Use a small input FASTA file.
- **Memory/compute constraints:**
  - Reduce `--max_len` if you encounter out-of-memory errors, but be aware this will increase the number of windows and total compute time.
  - Avoid using the largest models on CPU for large datasets.

#### DNABERT Model Selection

- By default, when you use `--type dnabert`, the script uses the model `zhihan1996/DNABERT-2-117M` (117M parameters) for embedding.
- You can override this by specifying `--model`, e.g. `--model zhihan1996/DNABERT-6` to use a different DNABERT variant.
- Both `zhihan1996/DNABERT-2-117M` (the main model) and `zhihan1996/DNABERT-6` (tokenizer/variant) are downloaded by `setup_offline.py` for compatibility, but only the one specified in the script or via `--model` is actually used for embedding.
- If you want to use a different DNABERT model by default, you can change the `default_model` entry in the `MODEL_CONFIGS` dictionary in `embed_fasta_transformer.py`.

#### Nucleotide Transformer Model Selection

- By default, when you use `--type nucleotide-transformer`, the script uses the model `InstaDeepAI/nucleotide-transformer-500m-human-ref` (500M parameters, human reference genome) for embedding.
- You can override this by specifying `--model`, e.g. `--model InstaDeepAI/nucleotide-transformer-2.5b-multi-species` to use a larger or multi-species variant.
- Both `InstaDeepAI/nucleotide-transformer-500m-human-ref` (the main model) and `InstaDeepAI/nucleotide-transformer-2.5b-multi-species` (larger/multi-species variant) are downloaded by `setup_offline.py` for compatibility, but only the one specified in the script or via `--model` is actually used for embedding.
- If you want to use a different Nucleotide Transformer model by default, you can change the `default_model` entry in the `MODEL_CONFIGS` dictionary in `embed_fasta_transformer.py`.

### 6. `test_segment.py`
- **Purpose:**
  - For each transcript in `data/smorfs_v2.csv`, extract a 196,608 bp window centered on the transcript, run the Nucleotide Transformer SegmentNT model, and analyse the region and any smORF(s) present.
  - Generate:
    - Probability plots (full and zoomed, with transcript and smORF locations, and gene context) in `plots/`
    - Per-transcript region analysis files in `region_analysis/`
    - A concatenated summary file for all transcripts in `region_analysis/all_transcripts_summary.txt`
    - Predicted gene models (GTF format, IGV-compatible) in `predicted_gtf/`
- **Input requirements:**
  - `data/smorfs_v2.csv` (with columns: transcript, primary_set, genomic_coordinates (1-based), ...)
  - `data/Homo_sapiens.GRCh38.114.gtf.gz` (Ensembl/GENCODE GTF)
  - `data/Homo_sapiens.GRCh38.dna.primary_assembly.fa` (genome FASTA)
- **Outputs:**
  - `plots/segment_plot_<transcript_id>.png` and `plots/segment_plot_<transcript_id>_zoom.png`: Annotated probability plots for each transcript
  - `region_analysis/<transcript_id>_region_probabilities.txt`: Per-transcript summary (transcript and smORF region statistics)
  - `region_analysis/all_transcripts_summary.txt`: All summary lines concatenated (tab-delimited)
  - `predicted_gtf/<transcript_id>_predicted.gtf`: Predicted gene model for transcript ±2kb, based on model probabilities
- **Usage:**
  ```bash
  python test_segment.py
  ```
- **Features:**
  - Batch processing of all transcripts in `smorfs_v2.csv` (with progress bar)
  - Model and reference files loaded only once for efficiency
  - Plots show transcript, smORF, and other gene locations
  - Region analysis includes both transcript and smORF region statistics, plus `primary_set` info
  - GTF output merges adjacent/overlapping features for IGV compatibility
  - All outputs are written to their respective folders for easy downstream analysis

### 7. `generate_qsub_scripts.py`
- **Purpose:** Automatically generates qsub job scripts for each main analysis step in the pipeline. Each script sets up the environment and runs the corresponding Python analysis script, making it easy to submit jobs to an HPC cluster using Sun Grid Engine (SGE) or compatible schedulers.
- **Output:** Shell scripts in the `qsub_scripts/` directory (one per analysis step), and log file paths in `logs/`.
- **Usage:**
  ```bash
  python generate_qsub_scripts.py
  # This creates qsub_scripts/<job>.sh for each main step
  # To submit a job, run:
  qsub qsub_scripts/<job>.sh
  # For example:
  qsub qsub_scripts/embed_fasta_transformer.sh
  ```
- **Customisation:**
  - Edit the `JOBS` list in `generate_qsub_scripts.py` to add, remove, or modify jobs and their arguments.
  - Edit the `QSUB_TEMPLATE` in the script to change resource requests (e.g., memory, runtime) or environment setup as needed for your cluster.
- **Note:**
  - The generated scripts assume you are using a conda environment named `canonical-fasta` and that your cluster uses SGE or a compatible system.
  - Log files for each job will be written to the `logs/` directory.

### 8. `embed_protein_fasta_transformer.py`
- **Purpose:** Embed all known protein sequences from Ensembl and all smORF protein sequences from `smorfs_v2.csv` using ProtTrans models (ProtBERT, ProtT5). Saves mean and max pooled embeddings, metadata, UMAP visualisation, and a nearest-neighbour network linking smORFs to their closest Ensembl proteins in embedding space.
- **Features:**
  - Downloads Ensembl protein FASTA if not present
  - Extracts unique protein sequences from Ensembl and smORFs
  - Embeds using selected ProtTrans model (`--type protbert` or `--type prott5`)
  - Saves embeddings and metadata to HDF5
  - Visualises UMAP (Ensembl vs smORF)
  - Finds k nearest Ensembl proteins for each smORF (network CSV)
  - **Test mode:** With `--test`, only 100 smORFs and 100 Ensembl proteins are embedded for quick validation
- **Input requirements:**
  - `data/Homo_sapiens.GRCh38.pep.all.fa.gz` (downloaded automatically if missing)
  - `data/smorfs_v2.csv` (must contain `releasev45_id` and `sequence_aa_MS` columns)
- **Outputs:**
  - `embeddings/protein_<model>_embeddings.h5`: Mean and max pooled embeddings, sequence metadata
  - `embeddings/protein_<model>_umap.png`: UMAP plot of mean embeddings, coloured by Ensembl/smORF
  - `embeddings/protein_<model>_smorf_ensembl_network.csv`: For each smORF, the k nearest Ensembl proteins in embedding space, with sequence and ID info
- **Usage:**
  ```bash
  python embed_protein_fasta_transformer.py --type protbert
  python embed_protein_fasta_transformer.py --type prott5
  # Add --test for a quick run on 100 smORFs and 100 Ensembl proteins
  python embed_protein_fasta_transformer.py --type protbert --test
  ```
- **Customisation:**
  - Use `--model` to specify a custom Hugging Face checkpoint
  - Use `--outdir` to change the output directory
  - The number of nearest neighbours (default 5) can be changed in the script

## Data Files

- **GTF/FASTA:** Downloaded automatically from Ensembl if not present, or can be pre-downloaded for offline use with `python setup_offline.py`.
- **smORFs:** Place `smorfs_v2.csv` in `data/`.

## Output Files

- **FASTA:**
  - `data/canonical_genes.fasta` — One transcript per gene (canonical/longest)
  - `data/all_transcripts.fasta` — All transcripts
- **Summary:**
  - `data/biotype_summary.txt` — Number of genes per biotype
- **Plots:**
  - `plots/` — Sequence statistics, classification, and UMAP plots
- **Embeddings:**
  - `embeddings/` — HDF5 files with sequence embeddings

## Notes
- All large/generated files in `data/`, `plots/`, and `embeddings/` are gitignored by default.
- For large-scale embedding, use a GPU for best performance.
- For questions or issues, see script docstrings or contact the maintainer.

## Canonical vs. All Transcripts FASTA: What is the Difference?

| Script                              | Output FASTA Contains                                                                 | Use Case / Rationale                                                                 |
|--------------------------------------|--------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------|
| `generate_canonical_fasta.py`        | **One transcript per gene**: the canonical transcript for protein-coding genes (as defined by Ensembl), and the longest transcript for non-coding genes | Use when you want a non-redundant set: one representative transcript per gene, suitable for gene-level analyses, summary statistics, or when you want to avoid transcript redundancy. |
| `generate_all_transcripts_fasta.py`  | **All annotated transcripts** from the GTF (including all isoforms for each gene)     | Use when you want to analyse every transcript isoform, e.g., for isoform-level expression, alternative splicing, or exhaustive transcriptome analysis. |

- **Canonical FASTA**: Smaller, non-redundant, one entry per gene. Good for summary statistics, gene-level plots, or when you want to avoid double-counting genes.
- **All Transcripts FASTA**: Larger, includes all isoforms. Good for isoform-level analysis, alternative splicing, or when you want to capture all transcript diversity.

---
