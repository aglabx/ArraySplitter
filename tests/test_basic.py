#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Basic working tests for ArraySplitter core functionality.
"""

import pytest
from ArraySplitter.decompose import (
    get_top1_nucleotide,
    compute_cuts,
    decompose_array,
    get_fs_tree,
)
from ArraySplitter.core_functions.tools.sequences import get_revcomp
from ArraySplitter.core_functions.tools.fs_tree import WeightedValueHeap


class TestBasicFunctionality:
    """Test basic ArraySplitter functionality."""
    
    def test_get_top1_nucleotide(self):
        """Test finding most frequent nucleotide."""
        assert get_top1_nucleotide("AAATC") == 'A'
        assert get_top1_nucleotide("TTTTC") == 'T'
        assert get_top1_nucleotide("CGCGC") == 'C'
        assert get_top1_nucleotide("GGGAT") == 'G'
    
    def test_simple_decomposition(self):
        """Test decomposition of simple repeat."""
        array = "ATCATCATCATC"
        monomers, _, _, _, period = decompose_array(array, depth=5, cutoff=1)
        
        # Should identify ATC as monomer
        assert period == 3
        assert len(monomers) == 4
        assert all(m == "ATC" for m in monomers)
    
    def test_dinucleotide_repeat(self):
        """Test decomposition of dinucleotide repeat."""
        array = "AT" * 50
        monomers, _, _, _, period = decompose_array(array, depth=5, cutoff=1)
        
        assert period == 2
        assert len(monomers) == 50
        assert all(m == "AT" for m in monomers)
    
    def test_get_revcomp(self):
        """Test reverse complement."""
        assert get_revcomp("ATCG") == "CGAT"
        assert get_revcomp("") == ""
        assert get_revcomp("AAAA") == "TTTT"
    
    def test_weighted_heap(self):
        """Test the heap data structure."""
        heap = WeightedValueHeap()
        heap.insert(3, "three")
        heap.insert(1, "one")
        heap.insert(2, "two")
        
        assert heap.pop() == (1, "one")
        assert heap.pop() == (2, "two")
        assert heap.pop() == (3, "three")
    
    def test_fs_tree_construction(self):
        """Test fs_tree can be built."""
        sequence = "ATCATCATC"
        fs_tree = get_fs_tree(sequence, 'A', cutoff=1)
        
        # Should create a tree
        assert isinstance(fs_tree, dict)
        assert 0 in fs_tree  # Root node
    
    def test_compute_cuts_basic(self):
        """Test cut computation."""
        array = "CAGCAGCAGCAG"
        hints = [(3, "CAG", 4)]
        
        cut_seq, score, period = compute_cuts(array, hints)
        
        assert cut_seq == "CAG"
        assert period == 3
        assert score > 0


class TestEdgeCases:
    """Test edge cases."""
    
    def test_single_repeat(self):
        """Test with single repeat unit."""
        array = "ATCGATCG"
        monomers, _, _, _, period = decompose_array(array, depth=10, cutoff=1)
        
        assert len(monomers) == 1
        assert monomers[0] == array
    
    def test_empty_hints(self):
        """Test compute_cuts with no hints."""
        array = "ATCG"
        hints = []
        
        cut_seq, score, period = compute_cuts(array, hints)
        
        # Should return the whole array as period
        assert score == 0
        assert period == len(array)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])