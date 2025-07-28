#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Improved rotation algorithm that:
1. Properly orients sequences to canonical form (A>T, C>G)
2. Uses cut sequences as conservative regions for rotation
"""

import sys
sys.path.insert(0, 'src')

from collections import Counter
from ArraySplitter.core_functions.tools.sequences import get_revcomp


def get_canonical_orientation(sequence):
    """
    Determine canonical orientation where A>T and C>G.
    Returns True if sequence is already canonical, False if needs reversal.
    """
    a_count = sequence.count('A')
    t_count = sequence.count('T')
    c_count = sequence.count('C')
    g_count = sequence.count('G')
    
    # Primary criterion: A > T
    if a_count != t_count:
        return a_count > t_count
    
    # Secondary criterion: C > G (when A == T)
    return c_count > g_count


def orient_sequences_canonically(arrays):
    """
    Orient all sequences to canonical form.
    Returns list of (oriented_sequence, was_reversed) tuples.
    """
    oriented = []
    
    for array in arrays:
        # Join monomers to analyze full sequence
        full_seq = "".join(array.split())
        
        if get_canonical_orientation(full_seq):
            # Already canonical
            oriented.append((array, False))
        else:
            # Need reverse complement
            rev_array = get_revcomp(array)
            oriented.append((rev_array, True))
    
    return oriented


def find_conserved_regions(arrays, k=8):
    """
    Find most conserved k-mers across arrays.
    These are good candidates for rotation anchors.
    """
    # Count k-mers across all arrays
    kmer_counts = Counter()
    kmer_positions = {}  # Track where each k-mer appears
    
    for array_idx, array in enumerate(arrays):
        monomers = array.split()
        full_seq = "".join(monomers)
        
        # Find all k-mers
        seen_in_array = set()
        for i in range(len(full_seq) - k + 1):
            kmer = full_seq[i:i+k]
            
            # Only count once per array
            if kmer not in seen_in_array:
                kmer_counts[kmer] += 1
                seen_in_array.add(kmer)
                
                if kmer not in kmer_positions:
                    kmer_positions[kmer] = []
                kmer_positions[kmer].append((array_idx, i))
    
    # Find k-mers that appear in most arrays
    total_arrays = len(arrays)
    conserved_kmers = []
    
    for kmer, count in kmer_counts.most_common(50):
        # Calculate conservation score
        presence_ratio = count / total_arrays
        
        # Check if it appears roughly once per monomer
        appearances = []
        for array_idx, _ in kmer_positions[kmer]:
            array = arrays[array_idx]
            monomers = array.split()
            monomer_count = sum(1 for m in monomers if kmer in m)
            appearances.append(monomer_count / len(monomers))
        
        avg_appearances = sum(appearances) / len(appearances) if appearances else 0
        
        # Good candidate: present in many arrays, ~1 per monomer
        if presence_ratio > 0.8 and 0.8 < avg_appearances < 1.2:
            conserved_kmers.append((kmer, presence_ratio, avg_appearances))
    
    return conserved_kmers


def rotate_with_cut_sequences(arrays, cut_sequences):
    """
    Use cut sequences from decomposition as rotation anchors.
    Cut sequences are by definition conserved regions.
    """
    # First, orient all sequences canonically
    oriented_data = orient_sequences_canonically(arrays)
    oriented_arrays = [seq for seq, _ in oriented_data]
    
    # Find which cut sequence is most common
    cut_counter = Counter()
    for array in oriented_arrays:
        for cut in cut_sequences:
            if cut in array:
                cut_counter[cut] += array.count(cut)
    
    if not cut_counter:
        print("No cut sequences found in arrays")
        return oriented_arrays
    
    # Use most common cut as rotation anchor
    best_cut, _ = cut_counter.most_common(1)[0]
    print(f"Using cut sequence '{best_cut}' as rotation anchor")
    
    # Rotate each array
    rotated_arrays = []
    for array in oriented_arrays:
        monomers = array.split()
        rotated_monomers = []
        
        for monomer in monomers:
            if best_cut in monomer:
                # Find first occurrence of cut
                pos = monomer.find(best_cut)
                if pos > 0:
                    # Rotate so cut is at start
                    rotated = monomer[pos:] + monomer[:pos]
                    rotated_monomers.append(rotated)
                else:
                    # Already starts with cut
                    rotated_monomers.append(monomer)
            else:
                # No cut in this monomer (flank?)
                rotated_monomers.append(monomer)
        
        rotated_arrays.append(" ".join(rotated_monomers))
    
    return rotated_arrays


def improved_rotate_arrays(arrays, cut_sequences=None, starting_kmer=None):
    """
    Improved rotation that:
    1. Orients to canonical form (A>T, C>G)
    2. Uses cut sequences or conserved regions for rotation
    """
    if cut_sequences:
        # Prefer cut sequences as they're proven conserved regions
        return rotate_with_cut_sequences(arrays, cut_sequences)
    
    # Orient canonically first
    oriented_data = orient_sequences_canonically(arrays)
    oriented_arrays = [seq for seq, _ in oriented_data]
    
    if starting_kmer:
        # Use provided k-mer
        best_kmer = starting_kmer
    else:
        # Find conserved regions
        conserved = find_conserved_regions(oriented_arrays)
        if conserved:
            best_kmer = conserved[0][0]  # Most conserved
            print(f"Found conserved k-mer: '{best_kmer}'")
        else:
            print("No conserved regions found")
            return oriented_arrays
    
    # Rotate using best k-mer
    rotated_arrays = []
    for array in oriented_arrays:
        monomers = array.split()
        rotated_monomers = []
        
        for monomer in monomers:
            if best_kmer in monomer:
                pos = monomer.find(best_kmer)
                if pos > 0:
                    rotated = monomer[pos:] + monomer[:pos]
                    rotated_monomers.append(rotated)
                else:
                    rotated_monomers.append(monomer)
            else:
                rotated_monomers.append(monomer)
        
        rotated_arrays.append(" ".join(rotated_monomers))
    
    return rotated_arrays


# Test the improved rotation
if __name__ == "__main__":
    # Test canonical orientation
    test_seqs = [
        "AAATTTCCCGGG",  # A=3, T=3, C=3, G=3 -> C>G decides
        "AAATTCCCGGG",   # A=3, T=2, C=3, G=3 -> A>T decides
        "ATATATCGCGCG",  # A=3, T=3, C=3, G=3 -> equal
    ]
    
    print("Testing canonical orientation:")
    for seq in test_seqs:
        is_canonical = get_canonical_orientation(seq)
        print(f"{seq}: {'canonical' if is_canonical else 'needs reversal'}")
        if not is_canonical:
            print(f"  -> {get_revcomp(seq)}")
    
    # Test with arrays containing cut sequences
    print("\n\nTesting with cut sequences:")
    test_arrays = [
        "ATGATGATGATG ATGATGATGATG",
        "TGATGATGATGA TGATGATGATGA",  # Rotated version
        "CATCATCATCAT CATCATCATCAT",  # Reverse complement
    ]
    
    cut_sequences = ["ATG", "CAT"]  # ATG and its reverse complement
    
    rotated = improved_rotate_arrays(test_arrays, cut_sequences=cut_sequences)
    print("\nRotated arrays:")
    for orig, rot in zip(test_arrays, rotated):
        print(f"Original: {orig}")
        print(f"Rotated:  {rot}")
        print()