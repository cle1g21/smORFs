import pytest
import tempfile
import os
from generate_canonical_fasta import parse_fasta, write_gene_fasta, write_biotype_summary

def test_parse_fasta(tmp_path):
    # Create a small gzipped FASTA file
    import gzip
    fasta_content = '>ENST0001\nATGC\n>ENST0002\nGCGC\n'
    fasta_path = tmp_path / 'test.fa.gz'
    with gzip.open(fasta_path, 'wt') as f:
        f.write(fasta_content)
    seqs = parse_fasta(str(fasta_path))
    assert seqs['ENST0001'] == 'ATGC'
    assert seqs['ENST0002'] == 'GCGC'

def test_write_gene_fasta(tmp_path):
    genes = {
        'gene1': {'transcript_id': 'ENST0001', 'gene_name': 'GENE1', 'chrom': '1', 'start': '1', 'end': '4', 'strand': '+', 'biotype': 'protein_coding'}
    }
    seqs = {'ENST0001': 'ATGC'}
    out_path = tmp_path / 'out.fa'
    write_gene_fasta(genes, seqs, str(out_path))
    with open(out_path) as f:
        lines = f.read().splitlines()
    assert lines[0].startswith('>GENE1|gene1')
    assert lines[1] == 'ATGC'

def test_write_biotype_summary(tmp_path):
    genes = {
        'gene1': {'biotype': 'protein_coding'},
        'gene2': {'biotype': 'lncRNA'}
    }
    out_path = tmp_path / 'summary.txt'
    write_biotype_summary(genes, str(out_path))
    with open(out_path) as f:
        content = f.read()
    assert 'protein_coding' in content
    assert 'lncRNA' in content 