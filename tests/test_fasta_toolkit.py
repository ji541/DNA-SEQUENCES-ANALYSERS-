# Test suite for the DNA sequence toolkit

from Bio.Seq import Seq

from fasta_toolkit import clean_seq, find_orfs, sequence_stats, translate_dna


class DummyRecord:
    def __init__(self, seq, identifier="test"):
        self.seq = seq
        self.id = identifier
        self.description = "test record"


def test_clean_seq_uppercases_and_replaces_rna_u():
    rec = DummyRecord("atg uaa")
    cleaned = clean_seq(rec)
    assert str(cleaned) == "ATGUAA"


def test_sequence_stats_counts_gc_and_length():
    rec = DummyRecord("ATGGCCNNNA")
    stats = sequence_stats(rec)
    assert stats["length"] == 11
    assert stats["gc_percent"] == 50.0
    assert stats["ambiguous"] == 3


def test_find_orfs_returns_longest_orf_for_basic_example():
    seq = Seq("ATGAAATAG")
    orfs = find_orfs(seq, min_aa=1, table=1)
    assert len(orfs) == 1
    assert orfs[0]["strand"] == "+"
    assert orfs[0]["aa_length"] == 2
    assert orfs[0]["protein"] == "MK"


def test_translate_dna_handles_reverse_frame_and_to_stop():
    seq = Seq("ATGAAATAG")
    protein = translate_dna(seq, frame=1, table=1, to_stop=False)
    assert protein == "MK*"

    rev = translate_dna(seq, frame=-1, table=1, to_stop=False)
    assert rev
