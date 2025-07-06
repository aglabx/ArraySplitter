#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tests for the main decomposition algorithm.
"""

import pytest
from ArraySplitter.decompose import (
    get_top1_nucleotide,
    compute_cuts,
    decompose_array,
    decompose_array_iter1,
    decompose_array_iter2,
)


class TestNucleotideSelection:
    """Test nucleotide selection functions."""
    
    def test_get_top1_nucleotide_simple(self):
        """Test finding most frequent nucleotide."""
        assert get_top1_nucleotide("AAATC") == 'A'
        assert get_top1_nucleotide("TTTTC") == 'T'
        assert get_top1_nucleotide("CGCGC") == 'C'
        assert get_top1_nucleotide("GGGAT") == 'G'
    
    def test_get_top1_nucleotide_tie(self):
        """Test behavior with tied frequencies."""
        # With ties, it should return one of them consistently
        result = get_top1_nucleotide("AATTCCGG")
        assert result in ['A', 'T', 'C', 'G']
    
    def test_get_top1_nucleotide_empty(self):
        """Test with empty sequence."""
        with pytest.raises(Exception):
            get_top1_nucleotide("")


class TestCutSequences:
    """Test cut sequence identification."""
    
    def test_split_array_simple(self):
        """Test splitting array by cut sequence."""
        array = "ATCGATCGATCG"
        cut = "CG"
        
        parts = array.split(cut)
        assert parts == ['AT', 'AT', 'AT', '']
    
    def test_split_array_no_cut(self):
        """Test when cut sequence not present."""
        array = "ATATATATAT"
        cut = "CG"
        
        parts = array.split(cut)
        assert parts == ['ATATATATAT']
    
    def test_split_array_edge_cases(self):
        """Test edge cases for cutting."""
        # Cut at beginning
        assert "CGATAT".split("CG") == ['', 'ATAT']
        
        # Cut at end
        assert "ATATCG".split("CG") == ['ATAT', '']
        
        # Multiple consecutive cuts
        assert "CGCGAT".split("CG") == ['', '', 'AT']


class TestComputeCuts:
    """Test optimal cut computation."""
    
    def test_compute_cuts_perfect_repeat(self):
        """Test with perfect tandem repeat."""
        array = "CAGCAGCAGCAGCAG"
        hints = [
            (1, "C", 5),
            (2, "CA", 5),
            (3, "CAG", 5),
        ]
        
        best_hint, best_score, best_parts = compute_cuts(array, hints)
        
        assert best_hint == "CAG"
        assert best_score == 1.0  # Perfect score
        assert len(best_parts) == 5
    
    def test_compute_cuts_imperfect_repeat(self):
        """Test with imperfect repeat."""
        array = "CAGCAGCATCAGCAG"  # One mutation
        hints = [
            (3, "CAG", 4),
            (3, "CAT", 1),
        ]
        
        best_hint, best_score, best_parts = compute_cuts(array, hints)
        
        assert best_hint == "CAG"
        assert best_score < 1.0  # Not perfect due to mutation


class TestDecomposition:
    """Test the main decomposition functions."""
    
    def test_decompose_array_simple_repeat(self):
        """Test decomposition of simple repeat."""
        array = "ATCATCATCATC"
        monomers = decompose_array(array, depth=5, cutoff=1)
        
        # Should identify ATC as monomer
        assert len(monomers) == 4
        assert all(m == "ATC" for m in monomers)
    
    def test_decompose_array_with_incomplete_end(self):
        """Test with incomplete repeat at end."""
        array = "ATCATCATCAT"  # Missing last C
        monomers = decompose_array(array, depth=5, cutoff=1)
        
        # Should handle incomplete end
        assert len(monomers) >= 3
        assert monomers[0] == "ATC"
    
    @pytest.mark.parametrize("repeat_unit,count", [
        ("AG", 10),
        ("CAG", 8),
        ("ATTCC", 5),
    ])
    def test_decompose_array_various_repeats(self, repeat_unit, count):
        """Test with various repeat units."""
        array = repeat_unit * count
        monomers = decompose_array(array, depth=10, cutoff=1)
        
        assert len(monomers) == count
        assert all(m == repeat_unit for m in monomers)
    
    def test_decompose_array_iter1(self):
        """Test first iteration of decomposition."""
        array = "CAGCAGCAGCAG"
        best_cut = "CAG"
        best_parts = ["", "", "", "", ""]
        period = 3
        
        monomers = decompose_array_iter1(array, best_cut, best_parts, period)
        
        assert len(monomers) == 4
        assert all(m == "CAG" for m in monomers)
    
    def test_decompose_array_iter2_no_refinement_needed(self):
        """Test second iteration when no refinement needed."""
        monomers = ["CAG", "CAG", "CAG", "CAG"]
        most_common = "CAG"
        
        refined = decompose_array_iter2(monomers, most_common)
        
        # Should not change perfect monomers
        assert refined == monomers
    
    def test_decompose_array_iter2_with_long_monomer(self):
        """Test second iteration with concatenated monomers."""
        monomers = ["CAG", "CAG", "CAGCAG", "CAG"]  # One double monomer
        most_common = "CAG"
        
        refined = decompose_array_iter2(monomers, most_common)
        
        # Should split the double monomer
        assert len(refined) == 5
        assert refined[2] == "CAG"
        assert refined[3] == "CAG"


class TestEdgeCases:
    """Test edge cases and error conditions."""
    
    def test_empty_array(self):
        """Test with empty array."""
        with pytest.raises(Exception):
            decompose_array("", depth=5, cutoff=1)
    
    def test_single_nucleotide(self):
        """Test with single nucleotide."""
        monomers = decompose_array("A", depth=5, cutoff=1)
        assert monomers == ["A"]
    
    def test_no_repeats(self):
        """Test with non-repetitive sequence."""
        array = "ACGTACGTACGT"
        monomers = decompose_array(array, depth=5, cutoff=5)
        
        # High cutoff should return whole sequence
        assert len(monomers) == 1
    
    def test_very_long_repeat_unit(self):
        """Test with very long repeat unit."""
        repeat_unit = "ATCGATCGATCGATCGATCG"  # 20bp
        array = repeat_unit * 3
        
        monomers = decompose_array(array, depth=25, cutoff=1)
        assert len(monomers) == 3
        assert all(m == repeat_unit for m in monomers)


@pytest.fixture
def satellite_array():
    """Generate a satellite array with some variation."""
    base_monomer = "ATTCCATTCCATTCC"
    # Add some variations
    monomers = [base_monomer] * 8
    monomers[3] = "ATTCCATTCTATTCC"  # Single mutation
    monomers[6] = "ATTCCATTCCATTC"   # Deletion
    return "".join(monomers)


def test_realistic_satellite(satellite_array):
    """Test with realistic satellite array."""
    monomers = decompose_array(satellite_array, depth=20, cutoff=3)
    
    # Should find approximately correct number of monomers
    assert 6 <= len(monomers) <= 10
    
    # Most monomers should be similar
    base = "ATTCCATTCCATTCC"
    similar_count = sum(1 for m in monomers if len(m) >= 14 and len(m) <= 16)
    assert similar_count >= len(monomers) * 0.7  # At least 70% similar length