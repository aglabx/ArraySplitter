


import argparse

import os
from collections import Counter
import re
from tqdm import tqdm

import editdistance as ed

from .core_functions.io.fasta_reader import \
    sc_iter_fasta_file
from .core_functions.io.satellome_reader import \
    sc_iter_satellome_file
from .core_functions.io.trf_reader import sc_iter_trf_file
from .core_functions.tools.fs_tree import \
    iter_fs_tree_from_sequence
from .core_functions.tools.sequences import get_revcomp


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

    ### TODO: The third iteration tries to glue short dangling monomers to the nearest monomer

    return decomposition, repeats2count, best_cut_seq, best_cut_score, best_period


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