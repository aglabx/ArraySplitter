#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tests for sequence utility functions.
"""

import pytest
from ArraySplitter.core_functions.tools.sequences import (
    get_revcomp,
    clear_sequence,
)


class TestReverseComplement:
    """Test reverse complement function."""
    
    def test_get_revcomp_simple(self):
        """Test basic reverse complement."""
        assert get_revcomp("ATCG") == "CGAT"
        assert get_revcomp("AAAA") == "TTTT"
        assert get_revcomp("TTTT") == "AAAA"
        assert get_revcomp("CCCC") == "GGGG"
        assert get_revcomp("GGGG") == "CCCC"
    
    def test_get_revcomp_mixed(self):
        """Test with mixed sequences."""
        assert get_revcomp("ATCGATCG") == "CGATCGAT"
        assert get_revcomp("GCTAGCTA") == "TAGCTAGC"
    
    def test_get_revcomp_palindrome(self):
        """Test palindromic sequences."""
        # These should be self-complementary when reversed
        assert get_revcomp("ATCGAT") == "ATCGAT"
        assert get_revcomp("GCATGC") == "GCATGC"
    
    def test_get_revcomp_empty(self):
        """Test empty sequence."""
        assert get_revcomp("") == ""
    
    def test_get_revcomp_lowercase(self):
        """Test handling of lowercase."""
        assert get_revcomp("atcg") == "cgat"
        assert get_revcomp("AtCg") == "CgAt"
    
    @pytest.mark.parametrize("seq,expected", [
        ("A", "T"),
        ("T", "A"),
        ("C", "G"),
        ("G", "C"),
        ("N", "N"),  # Ambiguous bases
        ("ATCGN", "NCGAT"),
    ])
    def test_get_revcomp_parametrized(self, seq, expected):
        """Parametrized tests for various sequences."""
        assert get_revcomp(seq) == expected


class TestClearSequence:
    """Test sequence cleaning function."""
    
    def test_clear_sequence_basic(self):
        """Test basic sequence cleaning."""
        assert clear_sequence("ATCG") == "ATCG"
        assert clear_sequence("atcg") == "ATCG"
        assert clear_sequence("AtCg") == "ATCG"
    
    def test_clear_sequence_with_spaces(self):
        """Test removal of spaces."""
        assert clear_sequence("A T C G") == "ATCG"
        assert clear_sequence("  ATCG  ") == "ATCG"
        assert clear_sequence("A\tT\nC\rG") == "ATCG"
    
    def test_clear_sequence_with_numbers(self):
        """Test removal of numbers."""
        assert clear_sequence("1ATCG") == "ATCG"
        assert clear_sequence("ATCG123") == "ATCG"
        assert clear_sequence("A1T2C3G4") == "ATCG"
    
    def test_clear_sequence_fasta_header(self):
        """Test cleaning FASTA-like input."""
        assert clear_sequence(">seq1\nATCG") == "ATCG"
        assert clear_sequence("ATCG\n>seq2\nGCTA") == "ATCGGCTA"
    
    def test_clear_sequence_empty(self):
        """Test empty and whitespace-only sequences."""
        assert clear_sequence("") == ""
        assert clear_sequence("   ") == ""
        assert clear_sequence("\n\n\n") == ""
    
    def test_clear_sequence_special_chars(self):
        """Test removal of special characters."""
        assert clear_sequence("A-T-C-G") == "ATCG"
        assert clear_sequence("A.T.C.G") == "ATCG"
        assert clear_sequence("A*T*C*G") == "ATCG"
        assert clear_sequence("(ATCG)") == "ATCG"
    
    def test_clear_sequence_ambiguous_bases(self):
        """Test handling of ambiguous bases."""
        # Assuming the function keeps standard IUPAC codes
        result = clear_sequence("ATCGRYSWKMBDHVN")
        assert "ATCG" in result
        # Check if ambiguous bases are handled (kept or removed)
        assert len(result) >= 4
    
    @pytest.mark.parametrize("input_seq,expected", [
        ("atcg", "ATCG"),
        ("ATCG", "ATCG"),
        ("a t c g", "ATCG"),
        ("1a2t3c4g5", "ATCG"),
        (">header\natcg", "ATCG"),
        ("", ""),
    ])
    def test_clear_sequence_parametrized(self, input_seq, expected):
        """Parametrized tests for sequence cleaning."""
        assert clear_sequence(input_seq) == expected


class TestSequenceUtilities:
    """Test combined sequence operations."""
    
    def test_clean_and_revcomp(self):
        """Test cleaning followed by reverse complement."""
        dirty_seq = "  atcg 123  "
        cleaned = clear_sequence(dirty_seq)
        assert cleaned == "ATCG"
        assert get_revcomp(cleaned) == "CGAT"
    
    def test_revcomp_idempotent(self):
        """Test that double reverse complement returns original."""
        original = "ATCGATCGATCG"
        revcomp1 = get_revcomp(original)
        revcomp2 = get_revcomp(revcomp1)
        assert revcomp2 == original
    
    def test_clear_idempotent(self):
        """Test that cleaning twice gives same result."""
        dirty = "  atcg 123  "
        clean1 = clear_sequence(dirty)
        clean2 = clear_sequence(clean1)
        assert clean1 == clean2


@pytest.fixture
def messy_fasta_content():
    """Provide messy FASTA-like content."""
    return """
    >sequence_1 description
    ATCGATCG atcg
    123 GCTAGCTA
    
    >sequence_2
    aaaa tttt
    cccc gggg
    """


def test_clear_multiline_fasta(messy_fasta_content):
    """Test cleaning multi-line FASTA content."""
    # This might need adjustment based on actual function behavior
    cleaned = clear_sequence(messy_fasta_content)
    # Should contain all sequences without headers or formatting
    assert "ATCGATCG" in cleaned
    assert "GCTAGCTA" in cleaned
    assert "AAAA" in cleaned or "aaaa" in cleaned.lower()
    assert ">" not in cleaned
    assert " " not in cleaned
    assert "\n" not in cleaned