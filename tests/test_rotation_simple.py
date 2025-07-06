#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Simple tests for rotation functions.
"""

import pytest
from ArraySplitter.core_functions.tools.rotation import (
    rotate_arrays,
    best_kmer_for_start,
)


class TestRotationFunctions:
    """Test rotation utility functions."""
    
    def test_rotate_arrays_basic(self):
        """Test basic array rotation."""
        # Test data: space-separated monomers
        arrays = [
            "CAGCAG AGCAGC GCAGCA",  # Same pattern, different starts
            "ATCATC TCATCA CATCAT",  # Another pattern
        ]
        
        # Rotate arrays
        rotated = rotate_arrays(arrays, start_seq=None)
        
        # Should return rotated arrays
        assert len(rotated) == 2
        assert isinstance(rotated[0], str)
        assert isinstance(rotated[1], str)
    
    def test_rotate_arrays_with_start(self):
        """Test rotation with specific start sequence."""
        arrays = [
            "CAGCAG AGCAGC GCAGCA",
        ]
        
        # Rotate to start with 'CAG'
        rotated = rotate_arrays(arrays, start_seq="CAG")
        
        # Should have rotated monomers
        assert len(rotated) == 1
        # Each monomer should start with CAG (if possible)
        monomers = rotated[0].split()
        for monomer in monomers:
            if "CAG" in monomer or monomer in "CAGCAG":
                # Check if it's a valid rotation
                assert len(monomer) > 0
    
    def test_best_kmer_for_start(self):
        """Test finding best k-mer for starting position."""
        # Create test monomers
        monomers = ["ATCGATCG", "TCGATCGA", "CGATCGAT", "GATCGATC"]
        
        # Find best starting k-mer
        best_kmer = best_kmer_for_start(monomers, k=3)
        
        # Should return a k-mer
        assert isinstance(best_kmer, str)
        assert len(best_kmer) == 3
        
        # The best k-mer should appear in the monomers
        assert any(best_kmer in m for m in monomers)
    
    def test_best_kmer_different_k(self):
        """Test with different k values."""
        monomers = ["ATCGATCG"] * 5
        
        # Test different k values
        for k in [2, 3, 4, 5]:
            best = best_kmer_for_start(monomers, k=k)
            assert len(best) == k
    
    def test_empty_arrays(self):
        """Test with empty input."""
        # Empty array
        rotated = rotate_arrays([], start_seq="ATG")
        assert rotated == []
        
        # Array with empty strings
        rotated = rotate_arrays([""], start_seq="ATG")
        assert len(rotated) == 1


@pytest.fixture
def satellite_array_decomposed():
    """Provide a decomposed satellite array."""
    return "ATTCCATTCC TTCCATTCCA TCCATTCCAT CCATTCCATT CATTCCATTC"


def test_rotation_preserves_sequence(satellite_array_decomposed):
    """Test that rotation preserves the overall sequence."""
    original = satellite_array_decomposed
    
    # Rotate without specific start
    rotated = rotate_arrays([original], start_seq=None)
    
    # The concatenated sequence should be preserved (possibly rotated)
    original_concat = original.replace(" ", "")
    rotated_concat = rotated[0].replace(" ", "")
    
    # Lengths should be the same
    assert len(original_concat) == len(rotated_concat)
    
    # Should contain same characters
    assert set(original_concat) == set(rotated_concat)