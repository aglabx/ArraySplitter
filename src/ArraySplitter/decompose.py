


import argparse

import os
from collections import Counter
import re
from tqdm import tqdm
from statistics import mean, stdev

import editdistance as ed

from .core_functions.io.fasta_reader import \
    sc_iter_fasta_file
from .core_functions.io.satellome_reader import \
    sc_iter_satellome_file
from .core_functions.io.trf_reader import sc_iter_trf_file
from .core_functions.tools.fs_tree import \
    iter_fs_tree_from_sequence
from .core_functions.tools.sequences import get_revcomp
from .core_functions.tools.anchor_graph import AnchorGraphDecomposer


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


def rotate_monomers_to_cut(decomposition, cut_sequence):
    """
    Rotate monomers so they start with the cut sequence.
    Returns rotated monomers.
    """
    rotated = []
    
    for monomer in decomposition:
        if monomer.startswith(cut_sequence):
            # Already starts with cut
            rotated.append(monomer)
        elif cut_sequence in monomer:
            # Find cut and rotate
            pos = monomer.find(cut_sequence)
            rotated_monomer = monomer[pos:] + monomer[:pos]
            rotated.append(rotated_monomer)
        else:
            # No cut (flank), keep as is
            rotated.append(monomer)
    
    return rotated


def get_top1_nucleotide(array):
    ### Step 1. Find the most frequent nucleotide (TODO: check all nucleotides and find with the best final score
    c = Counter()
    for n in "ACTG":
        c[n] = array.count(n)
        # print(n, array.count(n))
    return c.most_common(1)[0][0]


def get_fs_tree(array, top1_nucleotide, cutoff):
    ### Step 2. Build fs_tree (TODO:  optimize it for long sequences)
    names_ = [i for i in range(len(array)) if array[i] == top1_nucleotide]
    positions_ = names_[::]
    # print(f"Starting positions: {len(positions_)}")
    return iter_fs_tree_from_sequence(
        array, top1_nucleotide, names_, positions_, cutoff
    )


def is_self_repeating(pattern):
    """Check if a pattern is composed of repeated smaller units."""
    n = len(pattern)
    for sub_len in range(1, n // 2 + 1):
        if n % sub_len == 0:
            sub_pattern = pattern[:sub_len]
            if pattern == sub_pattern * (n // sub_len):
                return sub_pattern
    return None


def iterate_hints(array, fs_tree, depth):
    ### Step 3. Find a list of hints (hint is the sequenece for array cutoff)
    ### Modified to stop chains when self-repeating patterns are detected

    current_length = 0
    buffer = []
    found_patterns = {}  # Track patterns by their minimal unit
    
    for L, names, positions in fs_tree:
        if L != current_length:
            if buffer:
                max_n = 0
                found_seq = None
                for start, end, N in buffer:
                    if N > max_n:
                        max_n = N
                        found_seq = array[start : end + 1]
                
                # Check if this is a self-repeating pattern
                minimal_unit = is_self_repeating(found_seq)
                
                if minimal_unit:
                    # This is self-repeating, yield the minimal unit instead
                    # But only if we haven't yielded it before
                    if minimal_unit not in found_patterns:
                        # Find the frequency of the minimal unit
                        min_len = len(minimal_unit)
                        min_count = array.count(minimal_unit)
                        yield min_len, minimal_unit, min_count
                        found_patterns[minimal_unit] = True
                    # Don't yield the longer self-repeating pattern
                else:
                    # Not self-repeating, yield as normal
                    yield current_length, found_seq, max_n
                    found_patterns[found_seq] = True
            
            buffer = []
            current_length = L
            if current_length > depth:
                break
                
        start = names[0]
        end = positions[0]
        N = len(names)
        buffer.append((start, end, N))
    
    if buffer:
        max_n = 0
        found_seq = None
        for start, end, N in buffer:
            if N > max_n:
                max_n = N
                found_seq = array[start : end + 1]
        
        minimal_unit = is_self_repeating(found_seq)
        if minimal_unit and minimal_unit not in found_patterns:
            min_len = len(minimal_unit)
            min_count = array.count(minimal_unit)
            yield min_len, minimal_unit, min_count
        elif not minimal_unit:
            yield current_length, found_seq, max_n


def gcd(a, b):
    """Greatest common divisor."""
    while b:
        a, b = b, a % b
    return a


def find_gcd_of_list(numbers):
    """Find GCD of a list of numbers."""
    if not numbers:
        return 0
    result = numbers[0]
    for num in numbers[1:]:
        result = gcd(result, num)
        if result == 1:
            return 1
    return result


def compute_cuts(array, hints, score_threshold=0.05, fragmentation_threshold=0.5):
    ### Step 4. Find optimal cutoff with improved criteria
    
    candidates = []
    
    # Calculate metrics for each hint
    for L, cut_sequence, N in hints:
        parts = array.split(cut_sequence)
        periods = []
        non_empty_periods = []  # Track periods from non-empty parts
        
        for i, part in enumerate(parts):
            # Handle edge cases for first/last parts
            if i < len(parts) - 1 or part:  # All except possibly empty last
                period = len(part) + len(cut_sequence)
                periods.append(period)
                # Track non-empty parts separately
                if len(part) > 0:
                    non_empty_periods.append(period)
        
        if not periods:
            continue
        
        # Determine if this is a perfect/near-perfect repeat
        empty_ratio = (len(periods) - len(non_empty_periods)) / len(periods) if periods else 0
        
        if empty_ratio >= 0.8:  # 80% or more empty parts = perfect/near-perfect repeat
            # This is a perfect or near-perfect repeat
            # Use the cut sequence length as the period
            period_counts = Counter([len(cut_sequence)])
            mode_period = len(cut_sequence)
            mode_count = len(periods)
            total_segments = len(periods)
        elif non_empty_periods:
            # Use only non-empty parts for period calculation
            period_counts = Counter(non_empty_periods)
            mode_period, mode_count = period_counts.most_common(1)[0]
            total_segments = len(non_empty_periods)
        else:
            # Fallback: use all periods
            period_counts = Counter(periods)
            mode_period, mode_count = period_counts.most_common(1)[0]
            total_segments = len(periods)
        
        # Base score (uniformity)
        base_score = mode_count / total_segments
        
        # Fragmentation penalty
        short_threshold = mode_period * fragmentation_threshold
        short_fragments = sum(1 for p in periods if p < short_threshold)
        fragmentation = short_fragments / total_segments
        
        # Check for periodicity/divisibility
        unique_periods = list(period_counts.keys())
        period_gcd = find_gcd_of_list(unique_periods) if len(unique_periods) > 1 else mode_period
        
        candidates.append({
            'cut': cut_sequence,
            'mode_period': mode_period,
            'base_score': base_score,
            'fragmentation': fragmentation,
            'period_gcd': period_gcd,
            'period_distribution': period_counts,
            'num_segments': total_segments,
            'num_parts': len(parts),  # Total number of parts from split
            'empty_ratio': empty_ratio
        })
    
    if not candidates:
        return array, 0, len(array)
    
    # Group candidates by score
    candidates.sort(key=lambda x: x['base_score'], reverse=True)
    best_base_score = candidates[0]['base_score']
    
    # Get all candidates within threshold
    similar_candidates = [
        c for c in candidates 
        if c['base_score'] >= best_base_score - score_threshold
    ]
    
    # Check for fundamental period (GCD > 1)
    fundamental_candidates = []
    for c in similar_candidates:
        if c['period_gcd'] > 1 and c['period_gcd'] < c['mode_period']:
            # Check if most periods are multiples of GCD
            multiples = sum(1 for p in c['period_distribution'] 
                          if p % c['period_gcd'] == 0)
            if multiples >= c['num_segments'] * 0.8:  # 80% are multiples
                c['fundamental_period'] = c['period_gcd']
                fundamental_candidates.append(c)
    
    # If we found fundamental periods, use those
    if fundamental_candidates:
        best = min(fundamental_candidates, key=lambda x: x['fundamental_period'])
        return best['cut'], best['base_score'], best['fundamental_period']
    
    # Otherwise, penalize fragmentation and choose minimal period
    for c in similar_candidates:
        c['adjusted_score'] = c['base_score'] * (1 - c['fragmentation'] * 0.5)
    
    # Sort by adjusted score, then by number of segments (fewer is better), then by period (smaller is better)
    similar_candidates.sort(key=lambda x: (-x['adjusted_score'], x['num_segments'], x['mode_period']))
    
    best = similar_candidates[0]
    return best['cut'], best['base_score'], best['mode_period']


### Step 5a. Try to cut long monomers to expected
def refine_repeat_even(repeat, best_period):
    # Protection against destructive splitting when period=1
    if best_period <= 1 and len(repeat) > 1:
        # Don't split multi-character monomers into single nucleotides
        yield repeat
        return
    
    if len(repeat) % best_period == 0:
        start = 0
        for _ in range(len(repeat) // best_period):
            yield repeat[start : start + best_period]
            start += best_period
    else:
        yield repeat


def detect_and_fix_ab_pattern(decomposition, cut_seq, verbose=False):
    """
    Detect A-B alternating pattern and merge if it improves length uniformity.

    A-B pattern occurs when the cut sequence appears at TWO positions within each
    true monomer, creating an alternating pattern of short (A) and long (B) fragments.

    Decision is based on CV improvement, NOT sequence similarity.
    Different parts of the same monomer CAN have dissimilar sequences.

    Returns:
        (new_decomposition, was_merged)
    """
    if len(decomposition) < 4:
        return decomposition, False

    # Get lengths of monomers that start with cut sequence
    monomer_lengths = []
    monomer_indices = []

    for i, m in enumerate(decomposition):
        if m.startswith(cut_seq):
            monomer_lengths.append(len(m))
            monomer_indices.append(i)

    if len(monomer_lengths) < 4:
        return decomposition, False

    # Check for bimodal distribution
    length_counts = Counter(monomer_lengths)

    # Need at least 2 distinct length groups with significant counts
    if len(length_counts) < 2:
        return decomposition, False

    # Find the two most common length modes
    top_modes = length_counts.most_common(2)
    mode_a_len, mode_a_count = top_modes[0]
    mode_b_len, mode_b_count = top_modes[1]

    # Both modes should have significant representation (at least 15% each)
    total_monomers = len(monomer_lengths)
    if mode_a_count < total_monomers * 0.15 or mode_b_count < total_monomers * 0.15:
        return decomposition, False

    # Check that the two modes are significantly different (not just noise)
    # They should differ by at least 50% of the smaller one
    smaller_len = min(mode_a_len, mode_b_len)
    larger_len = max(mode_a_len, mode_b_len)

    if larger_len < smaller_len * 1.5:
        return decomposition, False

    if verbose:
        print(f"  Detected potential A-B pattern: {mode_a_len}bp ({mode_a_count}x) vs {mode_b_len}bp ({mode_b_count}x)")

    # Calculate current CV
    current_mean = mean(monomer_lengths)
    current_std = stdev(monomer_lengths) if len(monomer_lengths) > 1 else 0
    current_cv = current_std / current_mean if current_mean > 0 else 0

    # Check if adjacent pairs sum to consistent values
    # This is the key indicator of A-B pattern
    merged_lengths = []
    for i in range(0, len(monomer_indices) - 1, 2):
        idx1 = monomer_indices[i]
        idx2 = monomer_indices[i + 1]

        # Only merge adjacent fragments in the original decomposition
        if idx2 == idx1 + 1:
            merged_len = len(decomposition[idx1]) + len(decomposition[idx2])
            merged_lengths.append(merged_len)

    if len(merged_lengths) < 3:
        return decomposition, False

    # Calculate CV of merged lengths
    merged_mean = mean(merged_lengths)
    merged_std = stdev(merged_lengths) if len(merged_lengths) > 1 else 0
    merged_cv = merged_std / merged_mean if merged_mean > 0 else 0

    if verbose:
        print(f"  Current CV: {current_cv:.3f} (mean={current_mean:.0f})")
        print(f"  Merged CV:  {merged_cv:.3f} (mean={merged_mean:.0f})")

    # Merge if CV improves significantly (at least 50% reduction)
    if merged_cv >= current_cv * 0.5:
        if verbose:
            print(f"  CV improvement not significant enough, skipping merge")
        return decomposition, False

    if verbose:
        print(f"  ✓ Merging A-B pattern (CV improved by {(1 - merged_cv/current_cv)*100:.0f}%)")

    # Perform the merge
    original_sequence = "".join(decomposition)
    new_decomposition = []
    i = 0

    while i < len(decomposition):
        current = decomposition[i]

        # Check if this and next are both monomers that should be merged
        if (i < len(decomposition) - 1 and
            current.startswith(cut_seq) and
            decomposition[i + 1].startswith(cut_seq)):

            # Check if their lengths match the A-B pattern
            curr_len = len(current)
            next_len = len(decomposition[i + 1])

            # Check if one is close to mode_a and other is close to mode_b
            curr_is_a = abs(curr_len - smaller_len) / smaller_len < 0.2
            curr_is_b = abs(curr_len - larger_len) / larger_len < 0.2
            next_is_a = abs(next_len - smaller_len) / smaller_len < 0.2
            next_is_b = abs(next_len - larger_len) / larger_len < 0.2

            # Merge if one is A and other is B (in either order)
            if (curr_is_a and next_is_b) or (curr_is_b and next_is_a):
                merged = current + decomposition[i + 1]
                new_decomposition.append(merged)
                i += 2
                continue

        new_decomposition.append(current)
        i += 1

    # Verify reconstruction
    reconstructed = "".join(new_decomposition)
    if reconstructed != original_sequence:
        if verbose:
            print(f"  ERROR: Reconstruction failed after A-B merge, reverting")
        return decomposition, False

    return new_decomposition, True


def optimize_monomer_lengths(decomposition, cut_seq, verbose=True, array_id=None):
    """
    Post-processing optimization to merge short frequent monomers with adjacent longer ones.
    Goal: minimize variance of monomer lengths.

    This runs AFTER the main decomposition is complete.
    """
    from collections import Counter
    import editdistance as ed

    if len(decomposition) < 3:
        return decomposition
    
    # Calculate original sequence for verification
    original_sequence = "".join(decomposition)
    original_length = len(original_sequence)
    
    # Skip if first fragment doesn't start with cut (it's a flank)
    start_idx = 0
    if decomposition[0] and not decomposition[0].startswith(cut_seq):
        start_idx = 1
    
    # Get monomer lengths (excluding flanks)
    monomer_info = []  # List of (index, length) tuples
    for i in range(start_idx, len(decomposition)):
        # Skip likely right flank
        if i == len(decomposition) - 1:
            avg_len = sum(len(d) for d in decomposition[start_idx:i]) / max(1, i - start_idx)
            if len(decomposition[i]) < avg_len * 0.7:
                continue
        monomer_info.append((i, len(decomposition[i])))
    
    if len(monomer_info) < 2:
        return decomposition
    
    # Count length frequencies
    length_counts = Counter(info[1] for info in monomer_info)
    all_lengths = [info[1] for info in monomer_info]
    
    # Calculate initial variance
    initial_mean = mean(all_lengths)
    initial_variance = sum((x - initial_mean) ** 2 for x in all_lengths) / len(all_lengths)
    
    if verbose:
        print(f"Initial variance: {initial_variance:.1f} (mean={initial_mean:.1f})")
        print(f"Length distribution: {dict(sorted(length_counts.items())[:5])}...")
    
    # Find frequently occurring short monomers (at least 2 occurrences or 15% of monomers)
    min_frequency = max(2, int(len(monomer_info) * 0.15))
    short_lengths = [length for length, count in length_counts.items() 
                     if count >= min_frequency and length < initial_mean * 0.5]
    
    if not short_lengths:
        return decomposition
    
    if verbose:
        print(f"Frequent short lengths: {short_lengths}")
    
    # Try merging short monomers with adjacent ones
    working_decomposition = decomposition.copy()
    merge_occurred = True
    iteration = 0
    
    while merge_occurred and iteration < 1:  # Only one iteration to avoid over-merging
        merge_occurred = False
        iteration += 1
        new_decomposition = []
        i = 0
        
        while i < len(working_decomposition):
            if i < len(working_decomposition) - 1:
                current = working_decomposition[i]
                next_frag = working_decomposition[i + 1]
                
                # Check if current is a short frequent monomer or similar to one
                current_is_short = False
                if current.startswith(cut_seq) and next_frag.startswith(cut_seq):
                    if len(current) in short_lengths:
                        current_is_short = True
                    else:
                        # Check if it's within 5% of any frequent short length
                        for short_len in short_lengths:
                            if abs(len(current) - short_len) / short_len < 0.05:
                                current_is_short = True
                                break
                
                if current_is_short:
                    # Check if we're dealing with alternating different sequences (A B A B pattern)
                    # Compare the short and long fragments to see if they're different types
                    seq_short = current[len(cut_seq):min(len(current), len(cut_seq)+50)]
                    seq_long = next_frag[len(cut_seq):min(len(next_frag), len(cut_seq)+50)]
                    
                    if len(seq_short) > 0 and len(seq_long) > 0:
                        # Calculate similarity between short and long
                        similarity = 1 - (ed.eval(seq_short, seq_long) / max(len(seq_short), len(seq_long)))
                        
                        # If short and long are dissimilar, check for A-B pattern
                        if similarity < 0.8:  # Less than 80% similar means different types
                            # Count how many short and long fragments we have
                            short_count = 0
                            long_count = 0
                            
                            for j, frag in enumerate(working_decomposition):
                                if frag.startswith(cut_seq):
                                    frag_len = len(frag)
                                    # Check if it's a short fragment (similar length to current)
                                    if abs(frag_len - len(current)) / len(current) < 0.2:
                                        short_count += 1
                                    # Check if it's a long fragment (similar length to next)
                                    elif abs(frag_len - len(next_frag)) / len(next_frag) < 0.2:
                                        long_count += 1
                            
                            # If we have multiple instances of both short and long, it's likely A-B
                            if short_count >= 3 and long_count >= 3:
                                if verbose:
                                    print(f"  Detected A-B alternating pattern (dissimilar sequences), skipping merge of {len(current)}+{len(next_frag)}")
                                # Skip this merge
                                new_decomposition.append(current)
                                i += 1
                                continue
                    
                    # Calculate what variance would be after merge
                    test_lengths = []
                    for j, frag in enumerate(working_decomposition):
                        if j == i:
                            test_lengths.append(len(current) + len(next_frag))
                        elif j == i + 1:
                            continue
                        else:
                            # Only count monomers, not flanks
                            if frag.startswith(cut_seq) or (j > 0 and j < len(working_decomposition) - 1):
                                test_lengths.append(len(frag))
                    
                    if test_lengths:
                        test_mean = mean(test_lengths)
                        test_variance = sum((x - test_mean) ** 2 for x in test_lengths) / len(test_lengths)
                        
                        # For alternating patterns, also check coefficient of variation
                        test_cv = (test_variance ** 0.5) / test_mean if test_mean > 0 else 0
                        initial_cv = (initial_variance ** 0.5) / initial_mean if initial_mean > 0 else 0
                        
                        if verbose:
                            print(f"  Testing merge {len(current)}+{len(next_frag)}: CV {initial_cv:.3f} -> {test_cv:.3f}")
                        
                        # Accept merge if:
                        # 1. It reduces CV significantly, OR
                        # 2. We're merging a short fragment with a much longer one (likely overcutting)
                        length_ratio = len(next_frag) / len(current) if len(current) > 0 else 1
                        
                        if test_cv < initial_cv * 0.98 or length_ratio > 3:
                            merged = current + next_frag
                            new_decomposition.append(merged)
                            i += 2
                            merge_occurred = True
                            if verbose:
                                print(f"  ✓ Merged!")
                            continue
                
                # Try merging next with current if next is short
                else:
                    next_is_short = False
                    if next_frag.startswith(cut_seq) and current.startswith(cut_seq):
                        if len(next_frag) in short_lengths:
                            next_is_short = True
                        else:
                            # Check if it's within 5% of any frequent short length
                            for short_len in short_lengths:
                                if abs(len(next_frag) - short_len) / short_len < 0.05:
                                    next_is_short = True
                                    break
                    
                    if next_is_short:
                        # Don't merge if current is already long (avoid merging already merged monomers)
                        if len(current) > initial_mean * 1.5:
                            # Skip - current is already a merged monomer
                            pass
                        else:
                            # Calculate variance after merge
                            test_lengths = []
                            for j, frag in enumerate(working_decomposition):
                                if j == i:
                                    test_lengths.append(len(current) + len(next_frag))
                                elif j == i + 1:
                                    continue
                                else:
                                    if frag.startswith(cut_seq) or (j > 0 and j < len(working_decomposition) - 1):
                                        test_lengths.append(len(frag))
                    
                            if test_lengths:
                                test_mean = mean(test_lengths)
                                test_variance = sum((x - test_mean) ** 2 for x in test_lengths) / len(test_lengths)
                                
                                # Check coefficient of variation
                                test_cv = (test_variance ** 0.5) / test_mean if test_mean > 0 else 0
                                initial_cv = (initial_variance ** 0.5) / initial_mean if initial_mean > 0 else 0
                                
                                if verbose:
                                    print(f"  Testing merge {len(current)}+{len(next_frag)}: CV {initial_cv:.3f} -> {test_cv:.3f}")
                                
                                # Accept merge for short+long patterns
                                length_ratio = len(current) / len(next_frag) if len(next_frag) > 0 else 1
                                
                                if test_cv < initial_cv * 0.98 or length_ratio > 3:
                                    merged = current + next_frag
                                    new_decomposition.append(merged)
                                    i += 2
                                    merge_occurred = True
                                    if verbose:
                                        print(f"  ✓ Merged!")
                                    continue
            
            # No merge, keep fragment
            new_decomposition.append(working_decomposition[i])
            i += 1
        
        working_decomposition = new_decomposition
        
        # Update variance for next iteration
        current_lengths = [len(f) for f in working_decomposition 
                          if f.startswith(cut_seq) or working_decomposition.index(f) > 0]
        if current_lengths:
            current_mean = mean(current_lengths)
            initial_variance = sum((x - current_mean) ** 2 for x in current_lengths) / len(current_lengths)
            
            # Also update initial_mean for ratio calculations
            initial_mean = current_mean
    
    # Final verification - ensure perfect reconstruction
    final_sequence = "".join(working_decomposition)
    if final_sequence != original_sequence:
        print(f"ERROR: Sequence changed during merging! {original_length} != {len(final_sequence)}")
        print(f"Reverting to original decomposition")
        return decomposition
    
    if verbose and working_decomposition != decomposition:
        final_lengths = [len(f) for f in working_decomposition if f.startswith(cut_seq)]
        if final_lengths:
            final_mean = mean(final_lengths)
            final_variance = sum((x - final_mean) ** 2 for x in final_lengths) / len(final_lengths)
            print(f"Optimization complete: variance {sum((x - mean(all_lengths)) ** 2 for x in all_lengths) / len(all_lengths):.1f} -> {final_variance:.1f}")
    
    return working_decomposition


def find_mutant_anchor(sequence, anchor, expected_pos, window=50, max_dist=2):
    """
    Find anchor with mutations (including indels) near expected position.

    Args:
        sequence: The monomer sequence to search in
        anchor: The original anchor sequence
        expected_pos: Expected position of anchor (e.g., median monomer length)
        window: Search window around expected position
        max_dist: Maximum edit distance to accept

    Returns:
        Position of mutant anchor, or None if not found
    """
    anchor_len = len(anchor)
    start = max(0, expected_pos - window)
    end = min(len(sequence), expected_pos + window + anchor_len)

    best_pos = None
    best_dist = max_dist + 1
    best_pos_dist = float('inf')  # Distance from expected position

    for pos in range(start, end):
        # Check different lengths due to possible indels (±2bp)
        for length in range(anchor_len - 2, anchor_len + 3):
            if pos + length > len(sequence):
                continue
            candidate = sequence[pos:pos + length]
            dist = ed.eval(anchor, candidate)
            pos_dist = abs(pos - expected_pos)

            # Prefer: 1) lower edit distance, 2) closer to expected position
            if dist < best_dist or (dist == best_dist and pos_dist < best_pos_dist):
                best_dist = dist
                best_pos = pos
                best_pos_dist = pos_dist

    if best_dist <= max_dist:
        return best_pos
    return None


def split_long_monomers(decomposition, cut_seq, verbose=False, array_id=None):
    """
    Post-processing: split monomers that are 2x, 3x, etc. longer than expected.

    This happens when anchor has mutation and we missed a cut point.
    We search for mutant anchor near expected position and split there.
    Iterates until no more splits can be made.

    Args:
        decomposition: List of monomer sequences
        cut_seq: The cut/anchor sequence
        verbose: Print debug info
        array_id: Array identifier for logging

    Returns:
        New decomposition with long monomers split
    """
    if len(decomposition) < 3:
        return decomposition

    # Save original for verification
    original_sequence = "".join(decomposition)

    # Iterate until no more splits can be made
    max_iterations = 10  # Safety limit
    total_splits = 0
    current_decomposition = decomposition
    arr_info = f"[{array_id}] " if array_id else ""

    for iteration in range(max_iterations):
        # Calculate expected monomer length (median)
        # For short repeats, use ALL monomers (after split_duplicate_halves, monomers
        # may not start with original cut_seq anymore, e.g., 6bp pieces from 12bp split)
        lengths = [len(mono) for mono in current_decomposition]

        if not lengths:
            break

        # Use median to avoid outliers affecting expected length
        sorted_lengths = sorted(lengths)
        median_length = sorted_lengths[len(sorted_lengths) // 2]

        if verbose and iteration == 0:
            print(f"{arr_info}Split long monomers: median={median_length}bp, cut={cut_seq}")

        result = []
        splits_made = 0

        # For short repeats, find the most common monomer to use as reference
        if median_length < 23:
            mono_counts = Counter(current_decomposition)
            most_common_mono = mono_counts.most_common(1)[0][0] if mono_counts else None
        else:
            most_common_mono = None

        for mono_idx, mono in enumerate(current_decomposition):
            # Skip flanks (don't start with cut or most common monomer for short repeats)
            if median_length < 23:
                # For short repeats, skip if it's the most common monomer (already correct size)
                if mono == most_common_mono:
                    result.append(mono)
                    continue
            else:
                # For longer repeats, skip flanks that don't start with cut
                if not mono.startswith(cut_seq):
                    result.append(mono)
                    continue

            # Check if monomer is longer than expected
            ratio = len(mono) / median_length

            if ratio < 1.3:
                # Normal length monomer
                result.append(mono)
                continue

            # Try to split - either 2x/3x (mutant anchor) or 1.3-1.7x (exact anchor with tail)
            n_expected = max(2, round(ratio))

            if verbose:
                print(f"{arr_info}  #{mono_idx}: Long monomer {len(mono)}bp ({ratio:.1f}x), trying to split into {n_expected} parts")

            # Try to find mutant anchors and split
            parts = [mono]

            for split_num in range(1, n_expected):
                # Expected position of next anchor
                expected_pos = split_num * median_length

                # Search in the last (unsplit) part
                current_part = parts[-1]

                # Adjust expected_pos relative to current part
                offset = sum(len(p) for p in parts[:-1])
                rel_expected_pos = expected_pos - offset

                # Skip if expected position is too close to edges
                # For short repeats, allow position at exactly the expected spot
                min_pos = len(cut_seq) // 2  # Allow some flexibility for short repeats
                max_pos = len(current_part) - len(cut_seq) // 2
                if rel_expected_pos < min_pos or rel_expected_pos > max_pos:
                    continue

                # Find split position
                split_pos = None

                if median_length < 23:
                    # For short repeats, just split at median_length position
                    # No need to search for anchor - just use the expected position
                    if rel_expected_pos > 0 and rel_expected_pos < len(current_part):
                        split_pos = rel_expected_pos
                elif ratio >= 1.7:
                    # For 2x, 3x: search for mutant anchor with edit distance
                    search_window = max(5, int(median_length * 0.15))
                    split_pos = find_mutant_anchor(
                        current_part,
                        cut_seq,
                        rel_expected_pos,
                        window=search_window,
                        max_dist=2
                    )
                else:
                    # For 1.3-1.7x: split at median_length if monomer starts with exact anchor
                    # This handles "exact anchor + short tail" cases
                    if current_part.startswith(cut_seq) and rel_expected_pos >= len(cut_seq):
                        split_pos = rel_expected_pos

                if split_pos is not None and split_pos > 0:
                    # Split at found position
                    part1 = current_part[:split_pos]
                    part2 = current_part[split_pos:]

                    # Check variance criterion: split only if it reduces total deviation
                    # Without split: |len(current_part) - median|
                    # With split: |len(part1) - median| + |len(part2) - median|
                    dev_no_split = abs(len(current_part) - median_length)
                    dev_with_split = abs(len(part1) - median_length) + abs(len(part2) - median_length)

                    if dev_with_split >= dev_no_split:
                        if verbose:
                            print(f"{arr_info}    #{mono_idx}: Skip split at {split_pos}: deviation {dev_no_split} -> {dev_with_split} (no improvement)")
                        continue

                    # Replace last part with two new parts
                    parts[-1] = part1
                    parts.append(part2)
                    splits_made += 1

                    if verbose:
                        print(f"{arr_info}    #{mono_idx}: Split at pos {split_pos}: {len(part1)}bp + {len(part2)}bp (dev {dev_no_split} -> {dev_with_split})")

            result.extend(parts)

        # Update for next iteration
        total_splits += splits_made
        current_decomposition = result

        if splits_made == 0:
            # No more splits possible
            break

        if verbose:
            print(f"{arr_info}  Iteration {iteration + 1}: {splits_made} splits, continuing...")

    # Verify reconstruction
    final_sequence = "".join(current_decomposition)
    if final_sequence != original_sequence:
        print(f"ERROR: Sequence changed during split_long_monomers!")
        return decomposition

    if verbose and total_splits > 0:
        print(f"{arr_info}  Total: {total_splits} splits in {iteration + 1} iteration(s)")

    return current_decomposition


def split_duplicate_halves(decomposition, cut_seq, verbose=False, array_id=None):
    """
    Post-processing: split monomers that consist of two identical halves.

    This catches cases like CCTAACCCTAAC which is actually CCTAAC + CCTAAC
    (two rotated telomere units) that weren't caught by split_long_monomers
    because CCTAAC is not within edit distance of TAACCC.

    Iterates until no more splits can be made.

    Args:
        decomposition: List of monomer sequences
        cut_seq: The cut/anchor sequence
        verbose: Print debug info
        array_id: Array identifier for logging

    Returns:
        New decomposition with duplicate-half monomers split
    """
    if len(decomposition) < 3:
        return decomposition

    # Calculate expected monomer length (median)
    lengths = [len(mono) for mono in decomposition if mono.startswith(cut_seq)]
    if not lengths:
        return decomposition

    median_length = sorted(lengths)[len(lengths) // 2]

    # Only apply for short repeats where this pattern is common
    if median_length > 24:
        return decomposition

    arr_info = f"[{array_id}] " if array_id else ""

    # Save original for verification
    original_sequence = "".join(decomposition)

    # Iterate until no more splits
    max_iterations = 10
    total_splits = 0
    current_decomposition = decomposition

    for iteration in range(max_iterations):
        result = []
        splits_made = 0

        for mono_idx, mono in enumerate(current_decomposition):
            mono_len = len(mono)

            # For short repeats, check ANY even-length monomer for duplicate halves
            # This catches cases where median is 12bp but true unit is 6bp
            if mono_len % 2 != 0:
                result.append(mono)
                continue

            half_len = mono_len // 2
            first_half = mono[:half_len]
            second_half = mono[half_len:]

            # Check if both halves are identical or nearly identical
            if first_half == second_half:
                # Perfect duplicate - split!
                result.append(first_half)
                result.append(second_half)
                splits_made += 1

                if verbose:
                    print(f"{arr_info}  #{mono_idx}: Split duplicate halves: {mono_len}bp -> {half_len}bp + {half_len}bp ({first_half})")
            else:
                # Check edit distance for near-duplicates (allow 1 mismatch for sequences >= 8bp)
                dist = ed.eval(first_half, second_half)
                if dist <= 1 and half_len >= 4:
                    result.append(first_half)
                    result.append(second_half)
                    splits_made += 1

                    if verbose:
                        print(f"{arr_info}  #{mono_idx}: Split near-duplicate halves (dist={dist}): {mono_len}bp -> {half_len}bp + {half_len}bp")
                else:
                    result.append(mono)

        total_splits += splits_made
        current_decomposition = result

        if splits_made == 0:
            break

        if verbose:
            print(f"{arr_info}  Duplicate halves iteration {iteration + 1}: {splits_made} splits")

    # Verify reconstruction
    final_sequence = "".join(current_decomposition)
    if final_sequence != original_sequence:
        print(f"ERROR: Sequence changed during split_duplicate_halves!")
        return decomposition

    if verbose and total_splits > 0:
        print(f"{arr_info}  Total duplicate-half splits: {total_splits}")

    return current_decomposition


def split_by_popular_prefix(decomposition, cut_seq, verbose=False, array_id=None):
    """
    Post-processing: split monomers that start with the most popular monomer.

    For short repeats, if a monomer starts with the most common monomer sequence
    but is longer, split it into (popular_monomer) + (remainder).

    Example: If popular monomer is CCTAAC (6bp), then CCTAACCCTA (10bp)
    becomes CCTAAC (6bp) + CCTA (4bp).

    Args:
        decomposition: List of monomer sequences
        cut_seq: The cut/anchor sequence
        verbose: Print debug info
        array_id: Array identifier for logging

    Returns:
        New decomposition with prefix-matched monomers split
    """
    if len(decomposition) < 3:
        return decomposition

    # Find the most common monomer
    mono_counts = Counter(decomposition)
    if not mono_counts:
        return decomposition

    popular_mono, popular_count = mono_counts.most_common(1)[0]
    popular_len = len(popular_mono)

    # Only apply for short repeats where this pattern is common
    if popular_len > 20:
        return decomposition

    # Need significant majority to use as reference
    if popular_count < len(decomposition) * 0.3:
        return decomposition

    arr_info = f"[{array_id}] " if array_id else ""

    # Save original for verification
    original_sequence = "".join(decomposition)

    result = []
    splits_made = 0

    for mono_idx, mono in enumerate(decomposition):
        mono_len = len(mono)

        # Skip if it's the popular monomer or shorter
        if mono_len <= popular_len:
            result.append(mono)
            continue

        # Check if it starts with the popular monomer
        if mono.startswith(popular_mono):
            remainder = mono[popular_len:]
            remainder_len = len(remainder)

            # Variance check: only split if it reduces total deviation
            # Before: |mono_len - popular_len|
            # After: |popular_len - popular_len| + |remainder_len - popular_len| = 0 + |remainder_len - popular_len|
            dev_before = abs(mono_len - popular_len)
            dev_after = abs(remainder_len - popular_len)

            if dev_after < dev_before:
                # Split improves variance
                result.append(popular_mono)
                result.append(remainder)
                splits_made += 1

                if verbose:
                    print(f"{arr_info}  #{mono_idx}: Split by prefix: {mono_len}bp -> {popular_len}bp + {remainder_len}bp (dev {dev_before} -> {dev_after})")
            else:
                # Split would make variance worse, keep as is
                result.append(mono)
                if verbose:
                    print(f"{arr_info}  #{mono_idx}: Skip prefix split: {mono_len}bp (dev {dev_before} -> {dev_after}, no improvement)")
        else:
            result.append(mono)

    # Verify reconstruction
    final_sequence = "".join(result)
    if final_sequence != original_sequence:
        print(f"ERROR: Sequence changed during split_by_popular_prefix!")
        return decomposition

    if verbose and splits_made > 0:
        print(f"{arr_info}  Split {splits_made} monomers by popular prefix ({popular_mono})")

    return result


def decompose_array_iter1(array, best_cut_seq, best_period, verbose=True, array_id=None):
    """
    Decompose array using the cut sequence, ensuring perfect reconstruction.
    Cut sequence is the START of each monomer (except the first fragment).
    """
    repeats2count = Counter()
    decomposition = []
    
    if not best_cut_seq or best_cut_seq not in array:
        # No cut sequence or not found, return whole array
        decomposition.append(array)
        repeats2count[array] = 1
        return decomposition, repeats2count
    
    # Split by cut sequence
    parts = array.split(best_cut_seq)
    
    # Build monomers: cut + following part
    # First part is a special case (flank)
    if parts[0]:
        # First fragment (before first cut) - this is a flank
        decomposition.append(parts[0])
        repeats2count[parts[0]] += 1
        if verbose:
            print(f"Flank: {len(parts[0])}bp")
    
    # Process all other parts: cut + part = monomer
    for i in range(1, len(parts)):
        monomer = best_cut_seq + parts[i]
        decomposition.append(monomer)
        repeats2count[monomer] += 1
        if verbose:
            print(f"Monomer {i}: {len(monomer)}bp (cut {len(best_cut_seq)}bp + part {len(parts[i])}bp)")
    
    # Don't merge here - we'll do optimization at the very end of the pipeline
    
    # Verify reconstruction
    reconstructed = "".join(decomposition)
    if reconstructed != array:
        print(f"WARNING: Reconstruction mismatch! {len(array)} != {len(reconstructed)}")
        if array_id:
            print(f"  Sequence ID: {array_id}")
        # Print array identifier for debugging
        array_preview = array[:50] + "..." if len(array) > 50 else array
        print(f"  Array preview: {array_preview}")
        print(f"  Cut sequence: '{best_cut_seq}'")
    
    return decomposition, repeats2count


### Step 5b. Try to cut long monomers to expected
def refine_repeat_odd(repeat, best_period, most_common_monomer, verbose=False):
    if len(repeat) / best_period > 1.3:
        n = len(most_common_monomer)
        optimal_cut = 0
        best_ed = n

        begin_positions = [i for i in range(min(len(repeat) - n + 1, 5))]
        end_positions = [i for i in range(max(0, len(repeat) - n + 1 - 5), len(repeat) - n + 1)]

        for i in begin_positions+end_positions:
            rep_b = repeat[i : i + n]
            dist = ed.eval(most_common_monomer, rep_b)
            if dist < best_ed:
                best_ed = dist
                optimal_cut = i
                if verbose:
                    print(
                        "Optimal cut",
                        best_ed,
                        optimal_cut,
                        len(repeat[:optimal_cut]),
                        len(repeat[optimal_cut:]),
                    )
        if best_ed < n / 2:
            if optimal_cut == 0:
                optimal_cut += n
            a = repeat[:optimal_cut]
            b = repeat[optimal_cut:]
            if min(len(a), len(b)) < 0:  # n/3:
                yield repeat
            else:
                if a:
                    yield a
                if b:
                    yield b
        else:
            yield repeat
    else:
        yield repeat


def decompose_array_iter2(decomposition, best_period, repeats2count_ref, verbose=True):
    repeats2count = Counter()
    refined_decomposition = []
    most_common_monomer = None
    for monomer, tf in repeats2count_ref.most_common(1000):
        if len(monomer) == best_period:
            most_common_monomer = monomer
            break
    if not most_common_monomer:
        # No monomer with exact period length found
        # Try to find the most common monomer of any length
        if repeats2count_ref:
            most_common_monomer = repeats2count_ref.most_common(1)[0][0]
        else:
            # Use first monomer if nothing else available
            most_common_monomer = decomposition[0] if decomposition else ""
        
        if verbose:
            print(f"No monomer of length {best_period} found, using '{most_common_monomer}' (len={len(most_common_monomer)})")
    for repeat in decomposition:
        if verbose:
            print("Repeat under consideration", len(repeat), repeat)
        for repeat in refine_repeat_odd(
            repeat, best_period, most_common_monomer, verbose=verbose
        ):
            if verbose:
                print("Added:", len(repeat), repeat)
            repeats2count[repeat] += 1
            refined_decomposition.append(repeat)
    return (
        refined_decomposition,
        repeats2count,
        len(refined_decomposition) != len(decomposition),
    )


def print_monomers(decomposition, repeats2count, best_period):
    start2tf = Counter()
    for monomer in decomposition:
        start2tf[monomer[:5]] += 1
    print(start2tf)

    most_common_monomer = None
    for monomer, tf in repeats2count.most_common(1000):
        if len(monomer) == best_period:
            most_common_monomer = monomer
            break
    assert most_common_monomer
    for repeat in decomposition:
        print(
            len(repeat),
            start2tf[repeat[:5]],
            repeat,
            ed.eval(repeat, most_common_monomer),
        )


def print_pause_clean(decomposition, repeats2count, best_period):
    print_monomers(decomposition, repeats2count, best_period)
    input("?")


#   clear_output(wait=True)


def decompose_array_with_cuts(array, cut_sequences, verbose=False, array_id=None):
    """
    Decompose array using predefined cut sequences.
    Tries each cut and selects the best one based on scoring.
    """
    if not cut_sequences:
        raise ValueError("No cut sequences provided")
    
    # Try each cut sequence
    best_result = None
    best_score = -1
    
    for cut_seq in cut_sequences:
        if cut_seq not in array:
            continue
            
        # Split by this cut
        parts = array.split(cut_seq)
        
        # Calculate score (same logic as compute_cuts)
        periods = []
        for i, part in enumerate(parts):
            if i < len(parts) - 1 or part:
                period = len(part) + len(cut_seq)
                periods.append(period)
        
        if not periods:
            continue
            
        # Calculate score
        period_counts = Counter(periods)
        mode_period, mode_count = period_counts.most_common(1)[0]
        total_segments = len(periods)
        base_score = mode_count / total_segments
        
        # Fragmentation penalty
        short_threshold = mode_period * 0.5
        short_fragments = sum(1 for p in periods if p < short_threshold)
        fragmentation = short_fragments / total_segments
        adjusted_score = base_score * (1 - fragmentation * 0.5)
        
        if adjusted_score > best_score:
            best_score = adjusted_score
            best_result = (cut_seq, base_score, mode_period, len(parts))
    
    if best_result is None:
        # No cuts found, return whole array
        return [array], Counter({array: 1}), "", 0, len(array)
    
    cut_seq, score, period, num_parts = best_result
    
    # Decompose with best cut
    decomposition, repeats2count = decompose_array_iter1(
        array, cut_seq, period, verbose=verbose, array_id=array_id
    )
    
    return decomposition, repeats2count, cut_seq, score, period


def decompose_array(array, depth=500, cutoff=None, verbose=False, array_id=None):
    ### Step 0. Set cutoff based on array size if not provided
    if cutoff is None:
        if len(array) > 1_000_000:
            cutoff = 1000
        elif len(array) > 100_000:
            cutoff = 250
        elif len(array) > 10_000:
            cutoff = 10
        else:
            cutoff = 3
    
    ### Step 1-3. Get hints from all nucleotides instead of just the most frequent
    all_hints = []
    hint_sources = {}  # Track which nucleotide generated each hint
    
    for nucleotide in "ACTG":
        # Get positions of this nucleotide
        positions = [i for i in range(len(array)) if array[i] == nucleotide]
        
        if len(positions) <= cutoff:
            continue
            
        # Build fs_tree for this nucleotide
        fs_tree = get_fs_tree(array, nucleotide, cutoff=cutoff)
        
        # Get hints from this fs_tree
        for hint in iterate_hints(array, fs_tree, depth):
            all_hints.append(hint)
            hint_key = (hint[0], hint[1])  # (length, sequence)
            if hint_key not in hint_sources or hint[2] > hint_sources[hint_key][0]:
                hint_sources[hint_key] = (hint[2], nucleotide)
    
    # Remove duplicates, keeping the one with highest frequency
    unique_hints = {}
    for length, sequence, freq in all_hints:
        key = (length, sequence)
        if key not in unique_hints or freq > unique_hints[key][2]:
            unique_hints[key] = (length, sequence, freq)
    
    hints = list(unique_hints.values())
    
    ### Step 4. PRIMARY: Try anchor graph decomposition first
    ### Anchor graph selects the cycle with minimum CV (best length uniformity)
    graph_result = decompose_with_anchor_graph(array, hints, verbose=verbose)

    if graph_result is not None:
        decomposition, best_cut_seq, best_period, cv = graph_result
        if verbose:
            print(f"  Using anchor graph decomposition: {len(decomposition)} monomers, period={best_period}, CV={cv:.3f}")

        # Calculate best_cut_score as fraction of monomers starting with cut
        monomers_with_cut = sum(1 for m in decomposition if m.startswith(best_cut_seq))
        best_cut_score = monomers_with_cut / len(decomposition) if decomposition else 0
        repeats2count = Counter(decomposition)

        return decomposition, repeats2count, best_cut_seq, best_cut_score, best_period

    ### Step 5. FALLBACK: Use FS-tree decomposition if anchor graph fails
    if verbose:
        print(f"  Anchor graph failed, falling back to FS-tree decomposition")

    best_cut_seq, best_cut_score, best_period = compute_cuts(array, hints)

    decomposition, repeats2count = decompose_array_iter1(
        array, best_cut_seq, best_period, verbose=verbose, array_id=array_id
    )

    return decomposition, repeats2count, best_cut_seq, best_cut_score, best_period


def get_candidates_from_hints(array, hints):
    """
    Convert hints to candidates list with scores (same logic as compute_cuts but return all).
    Used for anchor graph decomposition.
    """
    candidates = []

    for L, cut_sequence, N in hints:
        parts = array.split(cut_sequence)
        periods = []
        non_empty_periods = []

        for i, part in enumerate(parts):
            if i < len(parts) - 1 or part:
                period = len(part) + len(cut_sequence)
                periods.append(period)
                if len(part) > 0:
                    non_empty_periods.append(period)

        if not periods:
            continue

        empty_ratio = (len(periods) - len(non_empty_periods)) / len(periods) if periods else 0

        if empty_ratio >= 0.8:
            period_counts = Counter([len(cut_sequence)])
            mode_period = len(cut_sequence)
            mode_count = len(periods)
            total_segments = len(periods)
        elif non_empty_periods:
            period_counts = Counter(non_empty_periods)
            mode_period, mode_count = period_counts.most_common(1)[0]
            total_segments = len(non_empty_periods)
        else:
            period_counts = Counter(periods)
            mode_period, mode_count = period_counts.most_common(1)[0]
            total_segments = len(periods)

        base_score = mode_count / total_segments
        short_threshold = mode_period * 0.5
        short_fragments = sum(1 for p in periods if p < short_threshold)
        fragmentation = short_fragments / total_segments
        adjusted_score = base_score * (1 - fragmentation * 0.5)

        candidates.append({
            'cut': cut_sequence,
            'length': L,
            'frequency': N,
            'mode_period': mode_period,
            'base_score': base_score,
            'adjusted_score': adjusted_score,
            'fragmentation': fragmentation,
            'num_segments': total_segments,
        })

    return candidates


def get_all_hints_for_graph(array, depth=100, cutoff=3):
    """
    Get hints from FS-tree for all nucleotides with small cutoff.
    Used for anchor graph analysis to find longer/rarer anchors.
    """
    all_hints = []

    for nucleotide in "ACTG":
        positions = [i for i in range(len(array)) if array[i] == nucleotide]

        if len(positions) <= cutoff:
            continue

        fs_tree = get_fs_tree(array, nucleotide, cutoff=cutoff)

        for hint in iterate_hints(array, fs_tree, depth):
            all_hints.append(hint)

    # Remove duplicates, keeping highest frequency
    unique = {}
    for length, anchor, freq in all_hints:
        key = (length, anchor)
        if key not in unique or freq > unique[key][2]:
            unique[key] = (length, anchor, freq)

    return list(unique.values())


def decompose_with_anchor_graph(array, hints, verbose=False):
    """
    Primary decomposition method using anchor graph.

    Anchor graph is built from FS-tree hints and selects the cycle
    with minimum CV (best length uniformity).

    Returns:
        (decomposition, cut_seq, period, cv) or None if failed
    """
    # For large sequences, get more hints with smaller cutoff
    if len(array) > 10000:
        graph_hints = get_all_hints_for_graph(array, depth=100, cutoff=3)
        if len(graph_hints) > len(hints):
            hints = graph_hints

    # Build candidates from hints
    candidates = get_candidates_from_hints(array, hints)

    if not candidates:
        return None

    # Build anchor graph
    decomposer = AnchorGraphDecomposer()
    decomposer.build_from_candidates(array, candidates, top_k=15, verbose=verbose)

    stats = decomposer.get_stats()
    cycle = stats['cycle']
    graph_period = stats['estimated_monomer_length']

    if not cycle or graph_period < 10:
        return None

    # Decompose using graph
    graph_decomposition = decomposer.decompose(verbose=verbose)

    # Verify reconstruction
    reconstructed = "".join(graph_decomposition)
    if reconstructed != array:
        if verbose:
            print(f"  Graph reconstruction failed")
        return None

    # Calculate CV of resulting decomposition
    cut_seq = cycle[0]
    lengths = [len(m) for m in graph_decomposition if m.startswith(cut_seq)]

    if len(lengths) < 2:
        return None

    mean_len = sum(lengths) / len(lengths)
    variance = sum((x - mean_len) ** 2 for x in lengths) / len(lengths)
    cv = (variance ** 0.5) / mean_len if mean_len > 0 else float('inf')

    if verbose:
        print(f"  Graph result: {len(graph_decomposition)} monomers, period={graph_period:.0f}, CV={cv:.3f}")

    return graph_decomposition, cut_seq, int(graph_period), cv


def check_and_fix_overcutting(array, decomposition, best_cut_seq, best_period, hints, verbose=False):
    """
    DEPRECATED: This function is kept for backwards compatibility.
    The anchor graph is now the primary decomposition method.
    """
    # This is now a no-op since graph is applied first
    return decomposition, best_cut_seq, best_period, False


def get_array_generator(input_file, format):
    '''Get array generator by format.'''
    if format == "fasta":
        return sc_iter_fasta_file(input_file)
    if format == "trf":
        return sc_iter_trf_file(input_file)
    if format == "satellome":
        return sc_iter_satellome_file(input_file)
    
    print(f"Unknown format: {format}")
    exit(1)
    

def main(input_file, output_prefix, format, threads, predefined_cuts=None, depth=100, verbose=False):
    """Main function."""

    sequences = get_array_generator(input_file, format)
    total = 0
    for _ in sequences:
        total += 1
    sequences = get_array_generator(input_file, format)

    print(f"Start processing")
    if predefined_cuts:
        print(f"Using predefined cuts: {', '.join(predefined_cuts)}")
    else:
        print(f"Will discover cuts automatically (depth={depth})")

    if output_prefix.endswith(".fasta"):
        print("Remove .fasta from output prefix")
        output_prefix = output_prefix[:-6]
    elif output_prefix.endswith(".fa"):
        print("Remove .fa from output prefix")
        output_prefix = output_prefix[:-3]

    output_file = f"{output_prefix}.decomposed.fasta"
    detail_file = f"{output_prefix}.monomers.tsv"
    lengths_file = f"{output_prefix}.lengths"
    print(f"Output file: {output_file}")
    print(f"Detail file: {detail_file}")
    print(f"Lengths file: {lengths_file}")
    
    # Open all output files
    with open(output_file, "w") as fw, open(detail_file, "w") as fw_detail, open(lengths_file, "w") as fw_lengths:
        # Write header for detail file
        fw_detail.write("sequence_id\torientation\tindex\ttype\tlength\tis_flank\tsequence\n")
        
        for header, array in tqdm(sequences, total=total):
            # Check canonical orientation
            is_canonical = get_canonical_orientation(array)
            was_reversed = False
            
            if not is_canonical:
                # Need reverse complement
                array = get_revcomp(array)
                was_reversed = True
            
            # print(len(array), end=" ")
            # Use predefined cuts or discover automatically
            if predefined_cuts:
                (
                    decomposition,
                    repeats2count,
                    best_cut_seq,
                    best_cut_score,
                    best_period,
                ) = decompose_array_with_cuts(array, predefined_cuts, verbose=verbose, array_id=header)
            else:
                # cutoff will be set automatically based on array size
                (
                    decomposition,
                    repeats2count,
                    best_cut_seq,
                    best_cut_score,
                    best_period,
                ) = decompose_array(array, depth=depth, cutoff=None, verbose=verbose, array_id=header)

            # Verify initial decomposition before any modifications
            reconstructed = "".join(decomposition)
            if reconstructed != array:
                print(f"FATAL ERROR: Initial decomposition does not match original sequence!")
                print(f"  Array ID: {header}")
                print(f"  Original length: {len(array)}")
                print(f"  Reconstructed length: {len(reconstructed)}")
                print(f"  Cut sequence: {best_cut_seq}")
                for i in range(min(len(array), len(reconstructed))):
                    if array[i] != reconstructed[i]:
                        print(f"  First mismatch at position {i}")
                        print(f"    Original: ...{array[max(0,i-10):i+20]}...")
                        print(f"    Reconstructed: ...{reconstructed[max(0,i-10):i+20]}...")
                        break
                raise ValueError(f"Initial decomposition verification failed for {header}")

            # Note: We do NOT rotate monomers internally - that would change the sequence.
            # "Rotation" means choosing a different cut point, not circular shifting within monomers.
            # decomposition = rotate_monomers_to_cut(decomposition, best_cut_seq)  # REMOVED - was breaking sequence

            # Detect and fix A-B alternating pattern (before other optimizations)
            decomposition, ab_merged = detect_and_fix_ab_pattern(decomposition, best_cut_seq, verbose=verbose)
            if ab_merged:
                # Recount monomers after A-B merge
                repeats2count = Counter(decomposition)

            # Apply post-processing optimization to merge short frequent monomers
            decomposition = optimize_monomer_lengths(decomposition, best_cut_seq, verbose=verbose, array_id=header)

            # Split monomers that are duplicate halves (e.g., CCTAACCCTAAC -> CCTAAC + CCTAAC)
            # Do this FIRST to establish correct median length before splitting long monomers
            decomposition = split_duplicate_halves(decomposition, best_cut_seq, verbose=verbose, array_id=header)

            # Split long monomers where anchor has mutation (x2, x3 longer than expected)
            # Now median is correct (e.g., 6bp instead of 12bp), so 11bp, 17bp etc. will be split
            decomposition = split_long_monomers(decomposition, best_cut_seq, verbose=verbose, array_id=header)

            # Final polish: split monomers that start with popular monomer
            # Example: CCTAACCCTA -> CCTAAC + CCTA if CCTAAC is the most common
            decomposition = split_by_popular_prefix(decomposition, best_cut_seq, verbose=verbose, array_id=header)

            # Final verification: decomposition must reconstruct original sequence exactly
            final_seq = "".join(decomposition)
            if final_seq != array:
                print(f"FATAL ERROR: Decomposition does not match original sequence!")
                print(f"  Array ID: {header}")
                print(f"  Original length: {len(array)}, Reconstructed length: {len(final_seq)}")
                for i in range(min(len(array), len(final_seq))):
                    if array[i] != final_seq[i]:
                        print(f"  First mismatch at position {i}")
                        print(f"    Original: ...{array[max(0,i-10):i+20]}...")
                        print(f"    Reconstructed: ...{final_seq[max(0,i-10):i+20]}...")
                        break
                raise ValueError(f"Decomposition verification failed for {header}")

            # Calculate statistics for internal monomers (excluding flanks)
            internal_monomers = []
            all_monomer_lengths = []
            
            # First pass: collect all potential monomer lengths
            for i, monomer in enumerate(decomposition):
                if monomer.startswith(best_cut_seq):
                    all_monomer_lengths.append(len(monomer))
            
            # Calculate average to better identify flanks
            if len(all_monomer_lengths) > 1:
                avg_monomer_len = sum(all_monomer_lengths) / len(all_monomer_lengths)
                flank_threshold = avg_monomer_len * 0.7  # 70% of average
            else:
                flank_threshold = best_period * 0.5
            
            # Second pass: identify true internal monomers
            for i, monomer in enumerate(decomposition):
                if monomer.startswith(best_cut_seq):
                    # Check if it's the last piece and too short (right flank)
                    if i == len(decomposition) - 1 and len(monomer) < flank_threshold:
                        continue  # Skip right flank
                    internal_monomers.append(monomer)
            
            if internal_monomers:
                internal_lengths = [len(m) for m in internal_monomers]
                min_len = min(internal_lengths)
                max_len = max(internal_lengths)
                avg_len = sum(internal_lengths) / len(internal_lengths)
                orientation = "rev" if was_reversed else "fwd"
                header_info = f"{header} cut={best_cut_seq} orientation={orientation} n_monomers={len(internal_monomers)} range={min_len}-{max_len} avg={avg_len:.1f}"
            else:
                orientation = "rev" if was_reversed else "fwd"
                header_info = f"{header} cut={best_cut_seq} orientation={orientation} n_monomers=0"
            
            fw.write(f">{header_info}\n")
            fw.write(" ".join(decomposition) + "\n")
            
            # Write lengths file
            fw_lengths.write(f">{header_info}\n")
            lengths = [str(len(m)) for m in decomposition]
            fw_lengths.write(" ".join(lengths) + "\n")
            
            # Write detailed monomer information
            for i, monomer in enumerate(decomposition):
                # LEFT_FLANK: only first fragment if it doesn't start with cut
                if i == 0 and not monomer.startswith(best_cut_seq):
                    piece_type = "LEFT_FLANK"
                    is_flank = "TRUE"
                # RIGHT_FLANK: only last fragment if too short
                elif i == len(decomposition) - 1 and len(monomer) < flank_threshold:
                    piece_type = "RIGHT_FLANK"
                    is_flank = "TRUE"
                else:
                    # Everything else is MONOMER (including mutant monomers without exact cut)
                    piece_type = "MONOMER"
                    is_flank = "FALSE"
                
                orientation = "rev" if was_reversed else "fwd"
                fw_detail.write(f"{header}\t{orientation}\t{i}\t{piece_type}\t{len(monomer)}\t{is_flank}\t{monomer}\n")


def run_it():
    parser = argparse.ArgumentParser(
        description="De novo decomposition of satellite DNA arrays into monomers"
    )
    parser.add_argument("-i", "--input", help="Input file", required=True)
    parser.add_argument(
        "--format",
        help="Input format: fasta, trf [fasta]",
        required=False,
        default="fasta",
    )
    parser.add_argument("-o", "--output", help="Output prefix", required=True)
    parser.add_argument(
        "-t", "--threads", help="Number of threads", required=False, default=4
    )
    parser.add_argument(
        "-c", "--cuts", 
        help="Comma-separated list of predefined cut sequences (e.g., ATG,ATGATG). If provided, skips cut discovery.", 
        required=False, 
        default=None
    )
    parser.add_argument(
        "-d", "--depth", 
        help="Depth for hint discovery (default: 100)", 
        required=False, 
        type=int,
        default=100
    )
    parser.add_argument(
        "-v", "--verbose", 
        help="Verbose output", 
        action="store_true"
    )
    args = parser.parse_args()

    input_file = args.input
    output_prefix = args.output
    format = args.format
    threads = int(args.threads)
    predefined_cuts = args.cuts.split(',') if args.cuts else None
    depth = args.depth
    verbose = args.verbose

    if not os.path.isfile(input_file):
        print(f"File {input_file} not found")
        exit(1)

    main(input_file, output_prefix, format, threads, predefined_cuts, depth, verbose)

if __name__ == "__main__":
    run_it()