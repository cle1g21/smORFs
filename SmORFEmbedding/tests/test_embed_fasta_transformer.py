import pytest
import numpy as np
import pandas as pd
from embed_fasta_transformer import fasta_to_df, scan_and_embed, MODEL_CONFIGS

class DummyTokenizer:
    def __call__(self, input_str, return_tensors=None):
        # Return a dummy tensor dict
        import torch
        return {'input_ids': torch.ones((1, 10), dtype=torch.long)}

class DummyModel:
    def to(self, device):
        return self
    def __call__(self, **inputs):
        import torch
        class Out:
            last_hidden_state = torch.ones((1, 10, 4))
        return Out()

def test_fasta_to_df(tmp_path):
    fasta_content = '>tx1 transcript_id=tx1 biotype=protein_coding\nATGC\n>tx2 transcript_id=tx2 biotype=lncRNA\nGCGC\n'
    fasta_path = tmp_path / 'test.fa'
    with open(fasta_path, 'w') as f:
        f.write(fasta_content)
    df = fasta_to_df(str(fasta_path))
    assert set(df.columns) == {'transcript_id', 'biotype', 'seq'}
    assert set(df['transcript_id']) == {'tx1', 'tx2'}
    assert set(df['biotype']) == {'protein_coding', 'lncRNA'}
    assert set(df['seq']) == {'ATGC', 'GCGC'}

def test_scan_and_embed():
    seq = 'ATGCATGCATGC'
    tokenizer = DummyTokenizer()
    model = DummyModel()
    device = 'cpu'
    model_type = 'dnabert'
    k = MODEL_CONFIGS[model_type]['k']
    max_len = 6
    mean_pool, max_pool = scan_and_embed(seq, tokenizer, model, device, model_type, k, max_len)
    assert mean_pool.shape == (4,)
    assert max_pool.shape == (4,)
    assert np.allclose(mean_pool, 1.0)
    assert np.allclose(max_pool, 1.0) 