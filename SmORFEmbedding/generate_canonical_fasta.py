"""
generate_canonical_fasta.py
--------------------------
Downloads the latest Ensembl GTF, cDNA, and ncRNA files and creates a FASTA with one transcript per gene: canonical for protein-coding, longest for non-coding. Also writes a biotype summary.
"""
#!/usr/bin/env python3
import argparse
import os
import sys
import requests
import gzip
import re
from pathlib import Path
from collections import defaultdict, Counter
from typing import Dict, Tuple
from Bio import SeqIO

ENSEMBL_FTP = "https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/"
ENSEMBL_FASTA = "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cdna/"
ENSEMBL_NCRNA = "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/ncrna/"
GTF_FILENAME = "Homo_sapiens.GRCh38.110.gtf.gz"  # Will update dynamically
FASTA_FILENAME = "Homo_sapiens.GRCh38.cdna.all.fa.gz"  # Will update dynamically
DATA_DIR = "data"


def download_file(url: str, dest: str):
    """
    Download a file from a URL if it does not already exist locally.
    Args:
        url (str): URL to download from.
        dest (str): Local destination path.
    """
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    if os.path.exists(dest):
        print(f"File {dest} already exists. Skipping download.")
        return
    print(f"Downloading {url} ...")
    r = requests.get(url, stream=True)
    r.raise_for_status()
    with open(dest, 'wb') as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)
    print(f"Downloaded to {dest}")


def get_latest_ensembl_filenames() -> Tuple[str, str, str]:
    """
    Scrape the Ensembl FTP directory for the latest GTF, cDNA, and ncRNA filenames.
    Returns:
        tuple: (gtf_file, fasta_file, ncrna_file)
    """
    gtf_url = ENSEMBL_FTP
    fasta_url = ENSEMBL_FASTA
    ncrna_url = ENSEMBL_NCRNA
    gtf_resp = requests.get(gtf_url)
    fasta_resp = requests.get(fasta_url)
    ncrna_resp = requests.get(ncrna_url)
    gtf_file = None
    fasta_file = None
    ncrna_file = None
    # Match the main GTF file (not abinitio, chr, or patch)
    gtf_pattern = re.compile(r'Homo_sapiens\.GRCh38\.\d+\.gtf\.gz')
    fasta_pattern = re.compile(r'Homo_sapiens\.GRCh38\.cdna\.all\.fa\.gz')
    ncrna_pattern = re.compile(r'Homo_sapiens\.GRCh38\.ncrna\.fa\.gz')
    for line in gtf_resp.text.splitlines():
        match = gtf_pattern.search(line)
        if match and 'abinitio' not in match.group(0) and 'chr' not in match.group(0) and 'patch' not in match.group(0):
            gtf_file = match.group(0)
    for line in fasta_resp.text.splitlines():
        match = fasta_pattern.search(line)
        if match:
            fasta_file = match.group(0)
    for line in ncrna_resp.text.splitlines():
        match = ncrna_pattern.search(line)
        if match:
            ncrna_file = match.group(0)
    if not gtf_file or not fasta_file or not ncrna_file:
        print("Could not determine latest Ensembl GTF, cDNA, or ncRNA file.")
        sys.exit(1)
    return gtf_file, fasta_file, ncrna_file


def parse_gtf_for_canonical_and_longest(gtf_path: str, seqs: dict) -> dict:
    """
    For protein-coding genes, use canonical transcript. For non-coding genes, use the longest transcript present in seqs.
    Returns a dict: gene_id -> {transcript_id, gene_name, chrom, start, end, strand, biotype, ...}
    Args:
        gtf_path (str): Path to GTF file (gzipped).
        seqs (dict): transcript_id -> sequence.
    Returns:
        dict: gene_id -> transcript info
    """
    genes = {}
    noncoding_transcripts = {}
    with gzip.open(gtf_path, 'rt') as gtf:
        for line in gtf:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attr = fields
            attrs = {k: v.strip('"') for k, v in [x.strip().split(' ', 1) for x in attr.strip(';').split('; ') if ' ' in x]}
            if feature == 'transcript':
                gene_id = attrs.get('gene_id')
                transcript_id = attrs.get('transcript_id')
                gene_name = attrs.get('gene_name', gene_id)
                biotype = attrs.get('gene_biotype', attrs.get('gene_type', 'NA'))
                tags = attrs.get('tag', '')
                if biotype == 'protein_coding' and 'Ensembl_canonical' in tags:
                    genes[gene_id] = {
                        'transcript_id': transcript_id,
                        'gene_name': gene_name,
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'biotype': biotype
                    }
                elif biotype != 'protein_coding':
                    # Only consider transcripts present in seqs
                    tid_short = transcript_id.split('.')[0]
                    if tid_short in seqs:
                        length = len(seqs[tid_short])
                        if gene_id not in noncoding_transcripts or length > noncoding_transcripts[gene_id]['length']:
                            noncoding_transcripts[gene_id] = {
                                'transcript_id': transcript_id,
                                'gene_name': gene_name,
                                'chrom': chrom,
                                'start': start,
                                'end': end,
                                'strand': strand,
                                'biotype': biotype,
                                'length': length
                            }
    # Merge noncoding genes
    for gene_id, info in noncoding_transcripts.items():
        genes[gene_id] = {k: v for k, v in info.items() if k != 'length'}
    return genes


def parse_fasta(fasta_path: str) -> Dict[str, str]:
    """
    Parse a FASTA file and return a dict of transcript_id -> sequence.
    Args:
        fasta_path (str): Path to FASTA file (gzipped).
    Returns:
        dict: transcript_id -> sequence
    """
    seqs = {}
    with gzip.open(fasta_path, 'rt') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            transcript_id = record.id.split('.')[0]  # Remove version
            seqs[transcript_id] = str(record.seq)
    return seqs


def write_gene_fasta(genes: Dict[str, Dict], seqs: Dict[str, str], out_path: str):
    """
    Write one transcript per gene to a FASTA file with informative headers.
    Args:
        genes (dict): gene_id -> transcript info.
        seqs (dict): transcript_id -> sequence.
        out_path (str): Output FASTA file path.
    """
    with open(out_path, 'w') as out:
        for gene_id, info in genes.items():
            tid = info['transcript_id'].split('.')[0]
            seq = seqs.get(tid)
            if not seq:
                continue
            header = f">{info['gene_name']}|{gene_id} chrom={info['chrom']}:{info['start']}-{info['end']} strand={info['strand']} biotype={info['biotype']} transcript_id={info['transcript_id']}"
            out.write(header + '\n')
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + '\n')


def write_biotype_summary(genes: dict, out_path: str):
    """
    Write a summary of gene biotype counts to a text file and print to stdout.
    Args:
        genes (dict): gene_id -> transcript info.
        out_path (str): Output summary file path.
    """
    biotypes = [info['biotype'] for info in genes.values()]
    counts = Counter(biotypes)
    with open(out_path, 'w') as f:
        f.write('Biotype\tGeneCount\n')
        for biotype, count in counts.items():
            f.write(f'{biotype}\t{count}\n')
    print('\nSummary statistics:')
    for biotype, count in counts.items():
        print(f'{biotype}: {count}')


def main():
    """
    Main workflow: downloads Ensembl files, parses GTF and FASTA, writes canonical/longest transcript per gene to FASTA, and writes biotype summary.
    """
    parser = argparse.ArgumentParser(description="Generate canonical gene FASTA from Ensembl GTF and cDNA/ncRNA.")
    parser.add_argument('--gtf', type=str, help='Custom GTF file (gzipped)')
    parser.add_argument('--cdna', type=str, help='Custom cDNA FASTA file (gzipped)')
    parser.add_argument('--ncrna', type=str, help='Custom ncRNA FASTA file (gzipped)')
    parser.add_argument('-o', '--output', type=str, default=os.path.join(DATA_DIR, 'canonical_genes.fasta'), help='Output FASTA file')
    args = parser.parse_args()

    gtf_file, fasta_file, ncrna_file = get_latest_ensembl_filenames()
    gtf_path = args.gtf or os.path.join(DATA_DIR, gtf_file)
    fasta_path = args.cdna or os.path.join(DATA_DIR, fasta_file)
    ncrna_path = args.ncrna or os.path.join(DATA_DIR, ncrna_file)

    if not args.gtf:
        download_file(ENSEMBL_FTP + gtf_file, gtf_path)
    if not args.cdna:
        download_file(ENSEMBL_FASTA + fasta_file, fasta_path)
    if not args.ncrna:
        download_file(ENSEMBL_NCRNA + ncrna_file, ncrna_path)

    print("Parsing cDNA FASTA...")
    seqs_cdna = parse_fasta(fasta_path)
    print(f"Loaded {len(seqs_cdna)} cDNA transcript sequences.")
    print("Parsing ncRNA FASTA...")
    seqs_ncrna = parse_fasta(ncrna_path)
    print(f"Loaded {len(seqs_ncrna)} ncRNA transcript sequences.")
    # Merge transcript dicts (ncRNA can overwrite cDNA for same transcript_id, which is rare)
    seqs = {**seqs_cdna, **seqs_ncrna}
    print(f"Total unique transcript sequences: {len(seqs)}\n")

    print("Parsing GTF for canonical and longest non-coding transcripts...")
    genes = parse_gtf_for_canonical_and_longest(gtf_path, seqs)
    print(f"Found {len(genes)} genes (canonical protein-coding + longest non-coding).")

    print(f"Writing output to {args.output} ...")
    write_gene_fasta(genes, seqs, args.output)
    # Write and print summary statistics
    summary_path = os.path.join(DATA_DIR, 'biotype_summary.txt')
    write_biotype_summary(genes, summary_path)
    print(f"Biotype summary written to {summary_path}")
    print("Done.")

if __name__ == "__main__":
    main() 