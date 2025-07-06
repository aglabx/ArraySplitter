#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test to verify the bug fix for independent sequence processing.
This test specifically checks that cutoff values don't leak between sequences.
"""

import pytest
import tempfile
from pathlib import Path
from ArraySplitter.decompose import main
from ArraySplitter.core_functions.io.fasta_reader import sc_iter_fasta_file


def test_cutoff_leak_bug_fixed(temp_dir):
    """
    Test that demonstrates the fix for cutoff value leaking between sequences.
    
    Before the fix: If a large sequence (>1MB) was processed first, its high cutoff (1000)
    would be used for all subsequent sequences, even small ones.
    
    After the fix: Each sequence gets its own appropriate cutoff based on its size.
    """
    # Create a multi-FASTA with a large sequence followed by small sequences
    # The large sequence will trigger cutoff=1000
    large_sequence = "ATCG" * 300000  # 1.2MB - should use cutoff=1000
    
    # These small sequences should use cutoff=3, not 1000
    small_sequence1 = "CAG" * 100     # 300bp - should use cutoff=3
    small_sequence2 = "AT" * 200      # 400bp - should use cutoff=3
    
    fasta_content = f""">large_array
{large_sequence}
>small_array1
{small_sequence1}
>small_array2  
{small_sequence2}
"""
    
    # Write test file
    test_file = temp_dir / "test_cutoff_bug.fa"
    test_file.write_text(fasta_content)
    
    # Process the file
    output_prefix = str(temp_dir / "output")
    main(str(test_file), output_prefix, "fasta", 1)
    
    # Read results
    output_file = Path(f"{output_prefix}.decomposed.fasta")
    results = {}
    for header, seq in sc_iter_fasta_file(str(output_file)):
        name = header.split()[0]
        monomers = seq.split()
        results[name] = monomers
    
    # Check that small sequences were properly decomposed
    # If the bug existed, high cutoff would prevent finding short repeats
    
    # Small array 1: CAG repeats
    assert len(results["small_array1"]) == 100, \
        f"Expected 100 CAG monomers, got {len(results['small_array1'])}"
    assert all(m == "CAG" for m in results["small_array1"]), \
        f"Expected all CAG monomers, got {set(results['small_array1'])}"
    
    # Small array 2: AT repeats  
    assert len(results["small_array2"]) == 200, \
        f"Expected 200 AT monomers, got {len(results['small_array2'])}"
    assert all(m == "AT" for m in results["small_array2"]), \
        f"Expected all AT monomers, got {set(results['small_array2'])}"
    
    print("✓ Bug fix verified: Each sequence uses appropriate cutoff")
    print(f"  - Large sequence ({len(large_sequence)}bp) → {len(results['large_array'])} monomers")
    print(f"  - Small sequence 1 ({len(small_sequence1)}bp) → {len(results['small_array1'])} monomers")
    print(f"  - Small sequence 2 ({len(small_sequence2)}bp) → {len(results['small_array2'])} monomers")


def test_reverse_order_consistency(temp_dir):
    """
    Test that processing order doesn't affect results.
    Process same sequences in different orders to ensure independence.
    """
    # Test data
    sequences = {
        "small": "GC" * 50,        # 100bp
        "medium": "ATTCC" * 3000,  # 15kb  
        "large": "CAG" * 50000,    # 150kb
    }
    
    # Test both orders
    orders = [
        ["small", "medium", "large"],
        ["large", "medium", "small"],
    ]
    
    all_results = []
    
    for i, order in enumerate(orders):
        # Create FASTA in specific order
        fasta_content = ""
        for name in order:
            fasta_content += f">{name}\n{sequences[name]}\n"
        
        # Process
        test_file = temp_dir / f"test_order_{i}.fa"
        test_file.write_text(fasta_content)
        
        output_prefix = str(temp_dir / f"output_order_{i}")
        main(str(test_file), output_prefix, "fasta", 1)
        
        # Collect results
        output_file = Path(f"{output_prefix}.decomposed.fasta")
        results = {}
        for header, seq in sc_iter_fasta_file(str(output_file)):
            name = header.split()[0]
            results[name] = seq.split()
        
        all_results.append(results)
    
    # Compare results from different orders
    for seq_name in sequences:
        result1 = all_results[0][seq_name]
        result2 = all_results[1][seq_name]
        
        assert result1 == result2, \
            f"Results differ for {seq_name} based on processing order:\n" \
            f"Order 1: {len(result1)} monomers\n" \
            f"Order 2: {len(result2)} monomers"
    
    print("✓ Processing order independence verified")


if __name__ == "__main__":
    # Allow running this test directly to verify the bug fix
    with tempfile.TemporaryDirectory() as tmpdir:
        temp_path = Path(tmpdir)
        test_cutoff_leak_bug_fixed(temp_path)
        test_reverse_order_consistency(temp_path)