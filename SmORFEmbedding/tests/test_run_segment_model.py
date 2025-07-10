import pytest
from run_segment_model import get_flanking_regions, predict_conservation

class DummyGenome:
    def __contains__(self, key):
        return key in ["chr1", "1"]
    def __getitem__(self, key):
        class DummySeq:
            @property
            def seq(self):
                return "A" * 10
            @property
            def reverse(self):
                return self
            @property
            def complement(self):
                return self
        return DummySeq()

def test_get_flanking_regions(monkeypatch):
    # Patch pyfaidx.Fasta to return DummyGenome
    import run_segment_model
    monkeypatch.setattr(run_segment_model.pyfaidx, "Fasta", lambda path: DummyGenome())
    # Patch gtfparse.read_gtf to return a minimal DataFrame
    import pandas as pd
    df = pd.DataFrame({
        'feature': ['transcript'],
        'transcript_id': ['tx1'],
        'seqname': ['1'],
        'start': [11],
        'end': [20],
        'strand': ['+']
    })
    monkeypatch.setattr(run_segment_model.gtfparse, "read_gtf", lambda path: df)
    result = get_flanking_regions("dummy.gtf", ["tx1"], "dummy.fa", upstream=5, downstream=5)
    assert "tx1_upstream" in result
    assert "tx1_downstream" in result
    assert result["tx1_upstream"] == "A" * 10
    assert result["tx1_downstream"] == "A" * 10

def test_predict_conservation():
    class DummyModel:
        def __call__(self, **inputs):
            class Out:
                logits = type('obj', (object,), {'mean': lambda self, dim: [0.5]})()
            return Out()
    class DummyTokenizer:
        def __call__(self, seqs, return_tensors, padding, truncation, max_length):
            return {"input_ids": [0]}
    result = predict_conservation(["ACGT"], DummyModel(), DummyTokenizer())
    assert result == [0.5] 