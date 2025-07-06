#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test to verify that decomposition can be reconstructed correctly.
This tests the fix for the reconstruction bug.
"""

import pytest
from ArraySplitter.decompose import decompose_array


class TestReconstruction:
    """Test that decomposition can be perfectly reconstructed."""
    
    def test_simple_reconstruction(self):
        """Test reconstruction of simple repeats."""
        test_cases = [
            "ATCATCATCATC",
            "CAGCAGCAGCAG",
            "AT" * 20,
            "GGGGAAATGGGGAAAT",
        ]
        
        for array in test_cases:
            monomers, _, cut_seq, _, period = decompose_array(array, depth=20, cutoff=1)
            
            # Try to reconstruct
            # Note: Current implementation loses cut sequences during split
            reconstructed = "".join(monomers)
            
            print(f"\nArray: {array[:30]}...")
            print(f"Cut sequence: '{cut_seq}'")
            print(f"Period: {period}")
            print(f"Monomers: {len(monomers)}")
            print(f"Reconstructed length: {len(reconstructed)} vs Original: {len(array)}")
            
            if reconstructed != array:
                print(f"MISMATCH!")
                print(f"Missing characters: {len(array) - len(reconstructed)}")
                
                # Find where they differ
                for i, (a, b) in enumerate(zip(array, reconstructed)):
                    if a != b:
                        print(f"First difference at position {i}: '{a}' vs '{b}'")
                        break
    
    def test_variable_repeat_reconstruction(self):
        """Test reconstruction of the problematic G-rich repeat."""
        array = "GGGGAAAATGGGGGGAAAATGGGAAAAATGGGAGGAAATTGGGGGAAATGGGGAAAAAATGGGGGAAAATGGGGAAAATTTGGGAGAAAATGGGGGGAAATGGGCGGGAAATGGGGAGAAATTGGGGAGGAAATGGGGGGGAAATGGGGGAAATGGGGAGAAATATGGGAAATTTTGTAAGGAAATGGGGAAAATATGGGAAAAAATTGTGGGGATATGGGGAGGAGAATGGGGGAAATGTGGGGAAAATGGGGGAGAAATGGAAGAGAAATTGTGGGGAAATGGGGGGAAAATATAGGGAATTTGGGGGGAAATGGGAGAGATATTGTGGGGAGATGTGGGGGGAAATGGGGGAGAAATTGGGGGGAAATGGGGAGAAATTGGGGGAAAATTAGGGGAAAATGGGGGGAAATACGGGAAAAATTGTGGGGAAATGGAGAAAATGTGGGGAAAATTGTGGGGAAATGGGAGAGATATTGTGGGGAGATGTAGGGGGAAATGGGGAAAATGGGGGAGAAAATGGGGGTAAAATGAGGGGAAATGGGAGGAAAATTGGGGGGAAAATGGGGGGAAATTGGGGGGGAAATGGGGGGAAAATCTGGGAAAAAATGTGGGAAATTTGGGGGGGAAAGGGGGGGAATGTGGGGGGATTTTGGGGGAAATGGGGGGAAATGGGGGGAAATACGGGAAAAATTGTGGGGAAATGGGGAAAATGTGGGGAAAATTGTGGGGAAATGGGGAGAGGAATGGGAGAAATGTGGGAAAAATGGGGGGATGGGAGAGAAATTGTGGGGAGATGTGGGGGGAAATGGGGAGGAAATATGGGGGGGAAATGGGGGAAAAACGTGGGGAAATGGGGAAGGAATGAAGGGGAAAATGGAAAAATGGGGGGGGAATGTGGGAAAATGAGGGGAAACAGAGAAAATGGGGAGGAATTGGGGGGAAATCGGGGAGAAATTGAGGGAAAATGGGGGAAATTGGGGAGAAATGAGGGCAAACTGGGGGGAAACGGGGAAAATTTGGGAGAAATTAGTGGGGAAATGAGGGGATAATGGTGGAAAATGAGGGGAAAT"
        
        # Try different parameters to see if any work
        params = [
            (10, None),
            (20, 5),
            (50, 10),
            (100, 20),
        ]
        
        for depth, cutoff in params:
            monomers, _, cut_seq, score, period = decompose_array(
                array, depth=depth, cutoff=cutoff
            )
            
            reconstructed = "".join(monomers)
            
            print(f"\nParameters: depth={depth}, cutoff={cutoff}")
            print(f"Cut: '{cut_seq}', Period: {period}, Score: {score:.3f}")
            print(f"Reconstruction match: {reconstructed == array}")
            
            if reconstructed != array:
                print(f"Length difference: {len(array) - len(reconstructed)}")


def test_decompose_array_iter1_reconstruction():
    """Test the iter1 function directly to understand the bug."""
    from ArraySplitter.decompose import decompose_array_iter1
    
    # Simple test case
    array = "ATCATCATCATC"
    cut_seq = "ATC"
    period = 3
    
    decomposition, counts = decompose_array_iter1(
        array, cut_seq, period, verbose=False
    )
    
    print(f"\nDirect iter1 test:")
    print(f"Array: {array}")
    print(f"Cut: {cut_seq}")
    print(f"Decomposition: {decomposition}")
    print(f"Reconstructed: {''.join(decomposition)}")
    print(f"Match: {''.join(decomposition) == array}")
    
    # Check what split gives us
    parts = array.split(cut_seq)
    print(f"\nSplit result: {parts}")
    print(f"Split loses the cut sequence!")


def analyze_split_behavior():
    """Analyze Python's split behavior to understand the bug."""
    
    print("\nPython split() behavior analysis:")
    
    test = "ATCATCATCATC"
    cut = "ATC"
    
    parts = test.split(cut)
    print(f"String: {test}")
    print(f"Cut: {cut}")
    print(f"Parts: {parts}")
    print(f"Reconstructed with cut: {cut.join(parts)}")
    print(f"Match: {cut.join(parts) == test}")
    
    # Edge cases
    print("\nEdge cases:")
    
    # Starts with cut
    test2 = "ATCGATCGATCG"
    parts2 = test2.split(cut)
    print(f"Starts with cut: {test2} -> {parts2}")
    
    # Ends with cut  
    test3 = "GATCGATCGATC"
    parts3 = test3.split(cut)
    print(f"Ends with cut: {test3} -> {parts3}")
    
    # Multiple consecutive cuts
    test4 = "ATCATCGATCATC"
    parts4 = test4.split(cut)
    print(f"Consecutive cuts: {test4} -> {parts4}")


if __name__ == "__main__":
    # Run analysis
    analyze_split_behavior()
    
    # Run tests
    test_decompose_array_iter1_reconstruction()
    
    # Run pytest
    pytest.main([__file__, "-v", "-s"])