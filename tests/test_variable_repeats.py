#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tests for ArraySplitter with highly variable repeats.
These are challenging cases that test the algorithm's robustness.
"""

import pytest
from ArraySplitter.decompose import decompose_array, get_top1_nucleotide


class TestVariableRepeats:
    """Test decomposition of variable/degenerate repeats."""
    
    @pytest.fixture
    def g_rich_variable_repeat(self):
        """Highly variable G-rich repeat from real data."""
        return "GGGGAAAATGGGGGGAAAATGGGAAAAATGGGAGGAAATTGGGGGAAATGGGGAAAAAATGGGGGAAAATGGGGAAAATTTGGGAGAAAATGGGGGGAAATGGGCGGGAAATGGGGAGAAATTGGGGAGGAAATGGGGGGGAAATGGGGGAAATGGGGAGAAATATGGGAAATTTTGTAAGGAAATGGGGAAAATATGGGAAAAAATTGTGGGGATATGGGGAGGAGAATGGGGGAAATGTGGGGAAAATGGGGGAGAAATGGAAGAGAAATTGTGGGGAAATGGGGGGAAAATATAGGGAATTTGGGGGGAAATGGGAGAGATATTGTGGGGAGATGTGGGGGGAAATGGGGGAGAAATTGGGGGGAAATGGGGAGAAATTGGGGGAAAATTAGGGGAAAATGGGGGGAAATACGGGAAAAATTGTGGGGAAATGGAGAAAATGTGGGGAAAATTGTGGGGAAATGGGAGAGATATTGTGGGGAGATGTAGGGGGAAATGGGGAAAATGGGGGAGAAAATGGGGGTAAAATGAGGGGAAATGGGAGGAAAATTGGGGGGAAAATGGGGGGAAATTGGGGGGGAAATGGGGGGAAAATCTGGGAAAAAATGTGGGAAATTTGGGGGGGAAAGGGGGGGAATGTGGGGGGATTTTGGGGGAAATGGGGGGAAATGGGGGGAAATACGGGAAAAATTGTGGGGAAATGGGGAAAATGTGGGGAAAATTGTGGGGAAATGGGGAGAGGAATGGGAGAAATGTGGGAAAAATGGGGGGATGGGAGAGAAATTGTGGGGAGATGTGGGGGGAAATGGGGAGGAAATATGGGGGGGAAATGGGGGAAAAACGTGGGGAAATGGGGAAGGAATGAAGGGGAAAATGGAAAAATGGGGGGGGAATGTGGGAAAATGAGGGGAAACAGAGAAAATGGGGAGGAATTGGGGGGAAATCGGGGAGAAATTGAGGGAAAATGGGGGAAATTGGGGAGAAATGAGGGCAAACTGGGGGGAAACGGGGAAAATTTGGGAGAAATTAGTGGGGAAATGAGGGGATAATGGTGGAAAATGAGGGGAAAT"
    
    def test_g_rich_decomposition(self, g_rich_variable_repeat):
        """Test decomposition of G-rich variable repeat."""
        array = g_rich_variable_repeat
        
        # This should not crash
        try:
            monomers, counts, cut_seq, cut_score, period = decompose_array(
                array, depth=50, cutoff=None, verbose=False
            )
            
            # Basic sanity checks
            assert len(monomers) > 0
            assert period > 0
            assert cut_seq is not None
            
            # The sequence should be reconstructable
            reconstructed = "".join(monomers)
            assert reconstructed == array
            
            # Print results for analysis
            print(f"\nG-rich repeat decomposition:")
            print(f"  Period: {period}")
            print(f"  Number of monomers: {len(monomers)}")
            print(f"  Cut sequence: '{cut_seq}'")
            print(f"  Cut score: {cut_score:.3f}")
            
            # Check monomer variability
            unique_monomers = len(set(monomers))
            print(f"  Unique monomers: {unique_monomers}")
            print(f"  Variability: {unique_monomers/len(monomers)*100:.1f}%")
            
        except Exception as e:
            pytest.fail(f"Decomposition failed with error: {e}")
    
    def test_nucleotide_composition(self, g_rich_variable_repeat):
        """Test nucleotide frequency analysis."""
        array = g_rich_variable_repeat
        
        # Check most frequent nucleotide
        top_nucl = get_top1_nucleotide(array)
        assert top_nucl == 'G'  # Should be G for this sequence
        
        # Count nucleotides
        counts = {n: array.count(n) for n in 'ACGT'}
        total = sum(counts.values())
        
        print(f"\nNucleotide composition:")
        for n, c in sorted(counts.items(), key=lambda x: -x[1]):
            print(f"  {n}: {c} ({c/total*100:.1f}%)")
    
    def test_pattern_analysis(self, g_rich_variable_repeat):
        """Analyze common patterns in the repeat."""
        array = g_rich_variable_repeat
        
        # Look for common short patterns
        patterns = {}
        for length in [4, 5, 6, 8, 10, 12]:
            pattern_count = {}
            for i in range(len(array) - length + 1):
                pattern = array[i:i+length]
                pattern_count[pattern] = pattern_count.get(pattern, 0) + 1
            
            # Find most common pattern of this length
            if pattern_count:
                most_common = max(pattern_count.items(), key=lambda x: x[1])
                patterns[length] = most_common
        
        print(f"\nMost common patterns by length:")
        for length, (pattern, count) in sorted(patterns.items()):
            print(f"  Length {length}: '{pattern}' (appears {count} times)")
    
    @pytest.mark.parametrize("depth,cutoff", [
        (10, None),   # Auto cutoff
        (20, 5),      # Higher cutoff
        (50, 10),     # Deep search
        (100, 20),    # Very deep search
    ])
    def test_different_parameters(self, g_rich_variable_repeat, depth, cutoff):
        """Test decomposition with different parameters."""
        array = g_rich_variable_repeat
        
        monomers, _, cut_seq, score, period = decompose_array(
            array, depth=depth, cutoff=cutoff, verbose=False
        )
        
        # Should always reconstruct correctly
        assert "".join(monomers) == array
        
        print(f"\nParameters depth={depth}, cutoff={cutoff}:")
        print(f"  Period: {period}, Monomers: {len(monomers)}")
        print(f"  Cut: '{cut_seq[:20]}...' Score: {score:.3f}")


class TestVariableRepeatEdgeCases:
    """Test edge cases with variable repeats."""
    
    def test_high_g_content(self):
        """Test with extremely high G content."""
        # 90% G content
        array = "GGGGGGGGGGAAGGGGGGGGGGAAGGGGGGGGGG"
        
        monomers, _, _, _, period = decompose_array(array, depth=10)
        assert len(monomers) > 0
        assert "".join(monomers) == array
    
    def test_variable_length_monomers(self):
        """Test repeat with variable length units."""
        # Mix of similar but different length patterns
        array = "GGGGAAATGGGGAAAATGGGGAAAAATGGGGAAAAAATGGGGAAAAAAATGGGGAAAAAAAT"
        
        monomers, _, _, _, period = decompose_array(array, depth=20)
        assert len(monomers) > 0
        assert "".join(monomers) == array
        
        # Check monomer lengths
        lengths = [len(m) for m in monomers]
        print(f"\nMonomer lengths: {lengths}")
    
    def test_with_interruptions(self):
        """Test repeat with random interruptions."""
        base = "GGGGAAAAT"
        interrupted = base + "CCC" + base + "TTT" + base + "AAA" + base
        
        monomers, _, cut_seq, _, period = decompose_array(interrupted, depth=15)
        assert len(monomers) > 0
        assert "".join(monomers) == interrupted
        
        print(f"\nInterrupted repeat:")
        print(f"  Cut: '{cut_seq}', Period: {period}")
        print(f"  Monomers: {monomers}")


def test_error_handling():
    """Test that algorithm handles errors gracefully."""
    
    # Very short sequence
    try:
        decompose_array("GGG", depth=5)
    except Exception as e:
        pytest.fail(f"Failed on short sequence: {e}")
    
    # All same nucleotide
    try:
        decompose_array("G" * 100, depth=5)
    except Exception as e:
        pytest.fail(f"Failed on homopolymer: {e}")
    
    # Empty hints case
    try:
        decompose_array("ACGT", depth=1, cutoff=100)  # High cutoff = no hints
    except Exception as e:
        pytest.fail(f"Failed with no hints: {e}")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])  # -s to see print statements