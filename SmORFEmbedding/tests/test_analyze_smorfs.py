import pytest
from analyze_smorfs import gc_content, sequence_entropy, has_kozak, kmer_freq

def test_gc_content():
    assert gc_content("GCGCGC") == 1.0
    assert gc_content("ATATAT") == 0.0
    assert gc_content("ATGC") == 0.5

def test_sequence_entropy():
    assert sequence_entropy("AAAA") == 0.0
    assert round(sequence_entropy("ATGC"), 2) == 2.0

def test_has_kozak():
    # Should return 1 if Kozak consensus is present
    assert has_kozak("GCCACCATGG") == 1
    # Should return 0 if not present
    assert has_kozak("ATGATGATG") == 0

def test_kmer_freq():
    freqs = kmer_freq("ATGCATGC", k=2)
    assert abs(freqs["AT"] - 2/7) < 1e-6
    assert abs(freqs["TG"] - 2/7) < 1e-6 