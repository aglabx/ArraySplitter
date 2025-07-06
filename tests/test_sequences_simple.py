#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Simple tests for sequence utility functions.
"""

import pytest
from ArraySplitter.core_functions.tools.sequences import (
    get_revcomp,
    REVCOMP_DICTIONARY,
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
        ("ATCGN", "NCGAT"),
    ])
    def test_get_revcomp_parametrized(self, seq, expected):
        """Parametrized tests for various sequences."""
        assert get_revcomp(seq) == expected
    
    def test_revcomp_idempotent(self):
        """Test that double reverse complement returns original."""
        original = "ATCGATCGATCG"
        revcomp1 = get_revcomp(original)
        revcomp2 = get_revcomp(revcomp1)
        assert revcomp2 == original
    
    def test_revcomp_dictionary(self):
        """Test the REVCOMP_DICTIONARY constant."""
        # Check basic mappings
        assert REVCOMP_DICTIONARY['A'] == 'T'
        assert REVCOMP_DICTIONARY['T'] == 'A'
        assert REVCOMP_DICTIONARY['C'] == 'G'
        assert REVCOMP_DICTIONARY['G'] == 'C'
        
        # Check lowercase
        assert REVCOMP_DICTIONARY['a'] == 't'
        assert REVCOMP_DICTIONARY['t'] == 'a'
        
        # Check ambiguous codes
        assert REVCOMP_DICTIONARY['N'] == 'N'
        assert REVCOMP_DICTIONARY['n'] == 'n'


@pytest.fixture
def test_sequences():
    """Provide test sequences."""
    return {
        "simple": "ATCG",
        "repeat": "ATATATAT",
        "complex": "ATCGATCGATCGATCG",
        "ambiguous": "ATCGNNNNATCG",
    }


def test_revcomp_preserves_length(test_sequences):
    """Test that reverse complement preserves sequence length."""
    for name, seq in test_sequences.items():
        revcomp = get_revcomp(seq)
        assert len(revcomp) == len(seq), f"Length mismatch for {name}"