#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tests for monomer rotation functions.
"""

import pytest
from ArraySplitter.core_functions.tools.rotation import (
    rotate_sequence,
    find_rotation,
    get_rotations,
)


class TestRotateSequence:
    """Test sequence rotation function."""
    
    def test_rotate_sequence_basic(self):
        """Test basic rotation."""
        seq = "ATCG"
        assert rotate_sequence(seq, 0) == "ATCG"
        assert rotate_sequence(seq, 1) == "TCGA"
        assert rotate_sequence(seq, 2) == "CGAT"
        assert rotate_sequence(seq, 3) == "GATC"
        assert rotate_sequence(seq, 4) == "ATCG"  # Full rotation
    
    def test_rotate_sequence_negative(self):
        """Test negative rotation (rotate left)."""
        seq = "ATCG"
        assert rotate_sequence(seq, -1) == "GATC"
        assert rotate_sequence(seq, -2) == "CGAT"
        assert rotate_sequence(seq, -3) == "TCGA"
        assert rotate_sequence(seq, -4) == "ATCG"
    
    def test_rotate_sequence_large(self):
        """Test rotation with positions > sequence length."""
        seq = "ATCG"
        assert rotate_sequence(seq, 5) == "TCGA"  # 5 % 4 = 1
        assert rotate_sequence(seq, 8) == "ATCG"  # 8 % 4 = 0
        assert rotate_sequence(seq, 10) == "CGAT" # 10 % 4 = 2
    
    def test_rotate_sequence_empty(self):
        """Test rotation of empty sequence."""
        assert rotate_sequence("", 0) == ""
        assert rotate_sequence("", 5) == ""
    
    def test_rotate_sequence_single(self):
        """Test rotation of single character."""
        assert rotate_sequence("A", 0) == "A"
        assert rotate_sequence("A", 1) == "A"
        assert rotate_sequence("A", -1) == "A"
    
    @pytest.mark.parametrize("seq,pos,expected", [
        ("ATCG", 0, "ATCG"),
        ("ATCG", 1, "TCGA"),
        ("ATCG", 2, "CGAT"),
        ("ATCGATCG", 3, "GATCGATC"),
        ("ABC", -1, "CAB"),
    ])
    def test_rotate_sequence_parametrized(self, seq, pos, expected):
        """Parametrized rotation tests."""
        assert rotate_sequence(seq, pos) == expected


class TestGetRotations:
    """Test getting all rotations of a sequence."""
    
    def test_get_rotations_basic(self):
        """Test getting all rotations."""
        seq = "ABC"
        rotations = get_rotations(seq)
        
        assert len(rotations) == 3
        assert "ABC" in rotations
        assert "BCA" in rotations
        assert "CAB" in rotations
    
    def test_get_rotations_longer(self):
        """Test rotations of longer sequence."""
        seq = "ATCG"
        rotations = get_rotations(seq)
        
        assert len(rotations) == 4
        assert rotations[0] == "ATCG"
        assert rotations[1] == "TCGA"
        assert rotations[2] == "CGAT"
        assert rotations[3] == "GATC"
    
    def test_get_rotations_duplicate(self):
        """Test rotations with repeating pattern."""
        seq = "ATAT"
        rotations = get_rotations(seq)
        
        assert len(rotations) == 4  # Still 4 rotations
        # But only 2 unique patterns
        unique_rotations = set(rotations)
        assert len(unique_rotations) == 2
        assert "ATAT" in unique_rotations
        assert "TATA" in unique_rotations
    
    def test_get_rotations_all_same(self):
        """Test rotations of homopolymer."""
        seq = "AAAA"
        rotations = get_rotations(seq)
        
        assert len(rotations) == 4
        assert all(r == "AAAA" for r in rotations)
    
    def test_get_rotations_empty(self):
        """Test rotations of empty sequence."""
        rotations = get_rotations("")
        assert rotations == []


class TestFindRotation:
    """Test finding best rotation to match a pattern."""
    
    def test_find_rotation_exact_match(self):
        """Test finding rotation with exact match."""
        seq = "CGATCG"
        pattern = "ATCG"
        
        rotation_pos = find_rotation(seq, pattern)
        rotated = rotate_sequence(seq, rotation_pos)
        assert rotated.startswith(pattern)
    
    def test_find_rotation_no_match(self):
        """Test when no rotation matches."""
        seq = "AAAA"
        pattern = "TTT"
        
        # Should return best effort (e.g., 0)
        rotation_pos = find_rotation(seq, pattern)
        assert isinstance(rotation_pos, int)
    
    def test_find_rotation_partial_match(self):
        """Test finding best partial match."""
        seq = "ATCGATCG"
        pattern = "TCGA"
        
        rotation_pos = find_rotation(seq, pattern)
        rotated = rotate_sequence(seq, rotation_pos)
        assert rotated.startswith("TCGA")
    
    def test_find_rotation_self(self):
        """Test finding rotation of sequence to itself."""
        seq = "ATCG"
        pattern = "ATCG"
        
        rotation_pos = find_rotation(seq, pattern)
        assert rotation_pos == 0
    
    def test_find_rotation_longer_pattern(self):
        """Test with pattern longer than sequence."""
        seq = "ATC"
        pattern = "ATCGATCG"
        
        # Should handle gracefully
        rotation_pos = find_rotation(seq, pattern)
        assert isinstance(rotation_pos, int)
        assert 0 <= rotation_pos < len(seq)


class TestRotationScenarios:
    """Test realistic rotation scenarios."""
    
    def test_normalize_tandem_monomers(self):
        """Test normalizing monomers from tandem repeat."""
        # Monomers from same repeat but different starting positions
        monomers = [
            "CAGCAGCAG",
            "AGCAGCAGC",
            "GCAGCAGCA",
        ]
        
        # All should rotate to same canonical form
        pattern = "CAG"
        normalized = []
        for monomer in monomers:
            pos = find_rotation(monomer, pattern)
            normalized.append(rotate_sequence(monomer, pos))
        
        # All should start with CAG
        assert all(m.startswith("CAG") for m in normalized)
    
    def test_rotation_preserves_length(self):
        """Test that rotation preserves sequence length."""
        sequences = ["A", "AT", "ATCG", "ATCGATCGATCG"]
        
        for seq in sequences:
            for i in range(len(seq) + 2):
                rotated = rotate_sequence(seq, i)
                assert len(rotated) == len(seq)
    
    def test_rotation_circular_property(self):
        """Test circular property of rotation."""
        seq = "ATCGATCG"
        
        # Rotating by length should return original
        assert rotate_sequence(seq, len(seq)) == seq
        
        # Rotating by i then by j == rotating by i+j
        for i in range(len(seq)):
            for j in range(len(seq)):
                rot_i = rotate_sequence(seq, i)
                rot_ij = rotate_sequence(rot_i, j)
                rot_direct = rotate_sequence(seq, i + j)
                assert rot_ij == rot_direct


@pytest.fixture
def satellite_monomers():
    """Provide realistic satellite monomers."""
    return [
        "ATTCCATTCCATTCC",
        "TTCCATTCCATTCCA",
        "TCCATTCCATTCCAT",
        "CCATTCCATTCCATT",
    ]


def test_normalize_satellite_monomers(satellite_monomers):
    """Test normalizing real satellite monomers."""
    # Find most common starting pattern
    pattern = "ATTCC"
    
    normalized = []
    for monomer in satellite_monomers:
        pos = find_rotation(monomer, pattern)
        normalized.append(rotate_sequence(monomer, pos))
    
    # All should now start with ATTCC
    assert all(m.startswith(pattern) for m in normalized)