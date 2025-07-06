#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tests to ensure independent processing of sequences in multi-FASTA files.
"""

import pytest
import tempfile
from pathlib import Path
from ArraySplitter.decompose import decompose_array, main
from ArraySplitter.core_functions.io.fasta_reader import sc_iter_fasta_file


class TestIndependentProcessing:
    """Test that sequences are processed independently in multi-FASTA."""
    
    def test_cutoff_independence(self):
        """Test that cutoff is set independently for each sequence."""
        # Create sequences of different sizes
        small_array = "CAG" * 100  # 300bp - should use cutoff=3
        large_array = "ATT" * 40000  # 120kb - should use cutoff=250
        
        # Process them in different orders
        results1 = []
        results2 = []
        
        # Order 1: small then large
        r1_small = decompose_array(small_array, depth=10, cutoff=None, verbose=False)
        r1_large = decompose_array(large_array, depth=10, cutoff=None, verbose=False)
        results1.append((r1_small, r1_large))
        
        # Order 2: large then small  
        r2_large = decompose_array(large_array, depth=10, cutoff=None, verbose=False)
        r2_small = decompose_array(small_array, depth=10, cutoff=None, verbose=False)
        results2.append((r2_small, r2_large))
        
        # Results should be the same regardless of order
        # Compare decompositions (first element of tuple)
        assert r1_small[0] == r2_small[0], "Small array results differ based on processing order"
        assert r1_large[0] == r2_large[0], "Large array results differ based on processing order"
    
    def test_multifasta_independence(self, temp_dir):
        """Test processing of multi-FASTA file."""
        # Create a multi-FASTA with different repeat types
        fasta_content = """>array1_dinucleotide
AT" * 50 + """
>array2_trinucleotide
CAG""" + "CAG" * 30 + """
>array3_complex
ATTCCATTCCATTCC""" + "ATTCCATTCCATTCC" * 10 + """
>array4_dinucleotide_again
GC""" + "GC" * 50
        
        # Write test file
        test_file = temp_dir / "test_multi.fa"
        test_file.write_text(fasta_content)
        
        # Process the file
        output_prefix = str(temp_dir / "output")
        main(str(test_file), output_prefix, "fasta", 1)
        
        # Read results
        output_file = Path(f"{output_prefix}.decomposed.fasta")
        assert output_file.exists()
        
        results = {}
        for header, seq in sc_iter_fasta_file(str(output_file)):
            results[header.split()[0]] = seq
        
        # Check that each array was decomposed correctly
        assert "array1_dinucleotide" in results
        assert "array2_trinucleotide" in results
        assert "array3_complex" in results
        assert "array4_dinucleotide_again" in results
        
        # Check decomposition quality
        # Array1: AT repeats
        monomers1 = results["array1_dinucleotide"].split()
        assert all(m == "AT" for m in monomers1), f"Unexpected monomers in array1: {set(monomers1)}"
        
        # Array2: CAG repeats
        monomers2 = results["array2_trinucleotide"].split()
        assert all(m == "CAG" for m in monomers2), f"Unexpected monomers in array2: {set(monomers2)}"
        
        # Array4: GC repeats (should be independent of array1 AT repeats)
        monomers4 = results["array4_dinucleotide_again"].split()
        assert all(m == "GC" for m in monomers4), f"Unexpected monomers in array4: {set(monomers4)}"
    
    def test_no_cross_contamination(self):
        """Test that patterns from one sequence don't affect another."""
        # Create two arrays with different patterns
        array1 = "ATCGATCG" * 20  # 8bp repeat
        array2 = "CAGCAGCAGCAG" * 15  # 3bp repeat
        
        # Process array1 first
        result1 = decompose_array(array1, depth=20, cutoff=3, verbose=False)
        monomers1 = result1[0]
        
        # Process array2 
        result2 = decompose_array(array2, depth=20, cutoff=3, verbose=False)
        monomers2 = result2[0]
        
        # Check no contamination
        # Array1 should have 8bp monomers
        assert all(len(m) == 8 for m in monomers1), f"Array1 monomer lengths: {[len(m) for m in monomers1]}"
        assert all("ATCGATCG" in m or m in "ATCGATCG" for m in monomers1)
        
        # Array2 should have 3bp monomers
        assert all(len(m) == 3 for m in monomers2), f"Array2 monomer lengths: {[len(m) for m in monomers2]}"
        assert all(m == "CAG" for m in monomers2)
    
    def test_parameter_isolation(self):
        """Test that each sequence gets appropriate parameters."""
        # Different sized arrays should get different cutoffs
        sizes_and_expected_cutoffs = [
            (500, 3),      # < 10kb
            (15000, 10),   # > 10kb
            (150000, 250), # > 100kb  
            (1500000, 1000), # > 1MB
        ]
        
        for size, expected_cutoff in sizes_and_expected_cutoffs:
            # Create array of specific size
            repeat_unit = "ATTCC"
            repeat_count = size // len(repeat_unit)
            array = repeat_unit * repeat_count
            
            # We can't directly test cutoff used, but we can verify
            # that results are consistent with expected behavior
            result = decompose_array(array, depth=10, cutoff=None, verbose=False)
            monomers = result[0]
            
            # Should find the repeat unit
            assert len(monomers) > 0
            # Most monomers should be the repeat unit
            correct_monomers = sum(1 for m in monomers if m == repeat_unit)
            assert correct_monomers >= len(monomers) * 0.8


@pytest.mark.integration 
class TestMultiFastaIntegration:
    """Integration tests for multi-FASTA processing."""
    
    def test_mixed_complexity_file(self, temp_dir):
        """Test file with mixed complexity arrays."""
        # Create complex test case
        fasta_content = """>perfect_short
""" + "AG" * 100 + """
>perfect_long
""" + "CAG" * 10000 + """
>imperfect_with_errors
""" + "ATTCC" * 50 + "ATTCT" + "ATTCC" * 49 + """
>complex_nested
""" + ("ATCGATCG" + "ATCGATCT") * 25 + """
>very_long
""" + "GCGCGC" * 200000
        
        test_file = temp_dir / "complex_test.fa"
        test_file.write_text(fasta_content)
        
        output_prefix = str(temp_dir / "complex_output")
        main(str(test_file), output_prefix, "fasta", 1)
        
        output_file = Path(f"{output_prefix}.decomposed.fasta")
        results = {}
        for header, seq in sc_iter_fasta_file(str(output_file)):
            name = header.split()[0]
            period = int(header.split()[1])
            results[name] = (seq.split(), period)
        
        # Verify each was processed correctly
        assert results["perfect_short"][1] == 2  # AG repeat
        assert results["perfect_long"][1] == 3   # CAG repeat
        assert results["imperfect_with_errors"][1] == 5  # ATTCC repeat
        assert results["very_long"][1] in [2, 3, 6]  # GC, GCG, or GCGCGC