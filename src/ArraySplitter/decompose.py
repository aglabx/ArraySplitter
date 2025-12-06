


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
    
    ### Step 4. Find the optimal cut sequence and the best period
    ### Defined as the maximal fraction of the cut sequence to the total cut sequence
    best_cut_seq, best_cut_score, best_period = compute_cuts(array, hints)

    ### Step 5. Cut the array
    ### The first iteration finds monomer frequencies
    # print("Firset iteration")
    decomposition, repeats2count = decompose_array_iter1(
        array, best_cut_seq, best_period, verbose=verbose, array_id=array_id
    )

    # assert "".join(decomposition) == array

    # DISABLED: The second iteration breaks the cut structure
    # It tries to refine monomers but creates pieces that don't start with cut
    # TODO: Fix decompose_array_iter2 to preserve cut structure

    # changed = True
    # while changed:
    #     # print("Firset iteration", len(decomposition))
    #     decomposition, repeats2count, changed = decompose_array_iter2(
    #         decomposition, best_period, repeats2count, verbose=verbose
    #     )
    #     # assert "".join(decomposition) == array

    ### Step 6. Check for overcutting using anchor graph
    decomposition, best_cut_seq, best_period, was_fixed = check_and_fix_overcutting(
        array, decomposition, best_cut_seq, best_period, hints, verbose=verbose
    )

    if was_fixed:
        # Recount monomers after graph-based decomposition
        repeats2count = Counter(decomposition)

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


def check_and_fix_overcutting(array, decomposition, best_cut_seq, best_period, hints, verbose=False):
    """
    Check if current decomposition shows signs of overcutting using anchor graph.

    Overcutting indicators:
    1. Multiple anchor cycle (>1 conserved parts per monomer)
    2. Graph period significantly larger than FS-tree period

    Returns:
        (decomposition, cut_seq, period, was_fixed)
    """
    if len(decomposition) < 3:
        return decomposition, best_cut_seq, best_period, False

    # For large sequences, get more hints with smaller cutoff for better anchor detection
    if len(array) > 10000:
        graph_hints = get_all_hints_for_graph(array, depth=100, cutoff=3)
        if len(graph_hints) > len(hints):
            hints = graph_hints

    # Build anchor graph from hints
    candidates = get_candidates_from_hints(array, hints)

    if not candidates:
        return decomposition, best_cut_seq, best_period, False

    decomposer = AnchorGraphDecomposer()
    decomposer.build_from_candidates(array, candidates, top_k=10, verbose=False)

    stats = decomposer.get_stats()
    graph_period = stats['estimated_monomer_length']
    cycle = stats['cycle']

    # Check for overcutting indicators
    is_overcutting = False
    reason = ""

    # Indicator 1: Multiple anchors in cycle (different conserved parts)
    if len(cycle) > 1:
        is_overcutting = True
        reason = f"multi-anchor cycle ({len(cycle)} anchors)"

    # Indicator 2: Graph period significantly larger (>2x)
    elif graph_period > best_period * 2 and graph_period > 50:
        is_overcutting = True
        reason = f"period mismatch (graph={graph_period:.0f} vs fstree={best_period})"

    if not is_overcutting:
        return decomposition, best_cut_seq, best_period, False

    # Try graph-based decomposition
    graph_decomposition = decomposer.decompose(verbose=False)

    # Verify reconstruction
    reconstructed = "".join(graph_decomposition)
    if reconstructed != array:
        if verbose:
            print(f"  Anchor graph reconstruction failed, keeping FS-tree result")
        return decomposition, best_cut_seq, best_period, False

    # Check if graph decomposition is better (fewer monomers with reasonable sizes)
    if len(graph_decomposition) < len(decomposition) * 0.9:
        if verbose:
            print(f"  Overcutting detected ({reason})")
            print(f"  FS-tree: {len(decomposition)} monomers, period={best_period}")
            print(f"  Graph:   {len(graph_decomposition)} monomers, period={graph_period:.0f}")
            print(f"  Using anchor graph decomposition")

        # Use first anchor in cycle as cut sequence
        new_cut_seq = cycle[0] if cycle else best_cut_seq
        return graph_decomposition, new_cut_seq, int(graph_period), True

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
            
            # Rotate monomers to start with cut
            decomposition = rotate_monomers_to_cut(decomposition, best_cut_seq)
            
            # Apply post-processing optimization to merge short frequent monomers
            decomposition = optimize_monomer_lengths(decomposition, best_cut_seq, verbose=verbose, array_id=header)

            # print("best period:", best_period, "len:", len(decomposition))
            # print_pause_clean(decomposition, repeats2count, best_period)

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
                if not monomer.startswith(best_cut_seq):
                    piece_type = "LEFT_FLANK"
                    is_flank = "TRUE"
                elif i == len(decomposition) - 1 and len(monomer) < flank_threshold:
                    piece_type = "RIGHT_FLANK"
                    is_flank = "TRUE"
                else:
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