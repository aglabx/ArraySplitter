//! Post-processing heuristics for decomposition refinement
//!
//! These functions refine the initial decomposition to handle:
//! - A-B alternating patterns (two anchors per monomer)
//! - Short frequent monomers that should be merged
//! - Long monomers that should be split
//! - Duplicate halves (e.g., ABAB -> AB + AB)
//! - Popular prefix patterns

use std::collections::HashMap;
use triple_accel::levenshtein::levenshtein_exp;

/// Calculate mean of a slice
fn mean(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.iter().sum::<f64>() / values.len() as f64
}

/// Calculate standard deviation of a slice
fn stdev(values: &[f64]) -> f64 {
    if values.len() < 2 {
        return 0.0;
    }
    let m = mean(values);
    let variance = values.iter().map(|x| (x - m).powi(2)).sum::<f64>() / values.len() as f64;
    variance.sqrt()
}

/// Calculate coefficient of variation
fn cv(values: &[f64]) -> f64 {
    let m = mean(values);
    if m == 0.0 {
        return 0.0;
    }
    stdev(values) / m
}

/// Count frequencies in a slice
fn count_frequencies<T: std::hash::Hash + Eq + Clone>(items: &[T]) -> HashMap<T, usize> {
    let mut counts = HashMap::new();
    for item in items {
        *counts.entry(item.clone()).or_insert(0) += 1;
    }
    counts
}

/// Edit distance between two strings
fn edit_distance(s1: &str, s2: &str) -> usize {
    levenshtein_exp(s1.as_bytes(), s2.as_bytes()) as usize
}

/// Detect and fix A-B alternating pattern
///
/// A-B pattern occurs when the cut sequence appears at TWO positions within each
/// true monomer, creating an alternating pattern of short (A) and long (B) fragments.
///
/// Decision is based on CV improvement, NOT sequence similarity.
pub fn detect_and_fix_ab_pattern(
    decomposition: &[String],
    cut_seq: &str,
    verbose: bool,
) -> (Vec<String>, bool) {
    if decomposition.len() < 4 {
        return (decomposition.to_vec(), false);
    }

    // Get lengths of monomers that start with cut sequence
    let mut monomer_lengths: Vec<usize> = Vec::new();
    let mut monomer_indices: Vec<usize> = Vec::new();

    for (i, m) in decomposition.iter().enumerate() {
        if m.starts_with(cut_seq) {
            monomer_lengths.push(m.len());
            monomer_indices.push(i);
        }
    }

    if monomer_lengths.len() < 4 {
        return (decomposition.to_vec(), false);
    }

    // Check for bimodal distribution
    let length_counts = count_frequencies(&monomer_lengths);

    if length_counts.len() < 2 {
        return (decomposition.to_vec(), false);
    }

    // Find the two most common length modes
    let mut counts_vec: Vec<(usize, usize)> = length_counts.into_iter().collect();
    counts_vec.sort_by(|a, b| b.1.cmp(&a.1));

    let (mode_a_len, mode_a_count) = counts_vec[0];
    let (mode_b_len, mode_b_count) = counts_vec[1];

    // Both modes should have significant representation (at least 15% each)
    let total_monomers = monomer_lengths.len();
    if mode_a_count < total_monomers * 15 / 100 || mode_b_count < total_monomers * 15 / 100 {
        return (decomposition.to_vec(), false);
    }

    // Check that the two modes are significantly different
    let smaller_len = mode_a_len.min(mode_b_len);
    let larger_len = mode_a_len.max(mode_b_len);

    if larger_len < smaller_len * 3 / 2 {
        return (decomposition.to_vec(), false);
    }

    if verbose {
        eprintln!(
            "  Detected potential A-B pattern: {}bp ({}x) vs {}bp ({}x)",
            mode_a_len, mode_a_count, mode_b_len, mode_b_count
        );
    }

    // Calculate current CV
    let lengths_f64: Vec<f64> = monomer_lengths.iter().map(|&x| x as f64).collect();
    let current_cv = cv(&lengths_f64);

    // Check if adjacent pairs sum to consistent values
    let mut merged_lengths: Vec<f64> = Vec::new();
    let mut i = 0;
    while i < monomer_indices.len() - 1 {
        let idx1 = monomer_indices[i];
        let idx2 = monomer_indices[i + 1];

        // Only merge adjacent fragments
        if idx2 == idx1 + 1 {
            let merged_len = decomposition[idx1].len() + decomposition[idx2].len();
            merged_lengths.push(merged_len as f64);
        }
        i += 2;
    }

    if merged_lengths.len() < 3 {
        return (decomposition.to_vec(), false);
    }

    // Calculate CV of merged lengths
    let merged_cv = cv(&merged_lengths);

    if verbose {
        eprintln!("  Current CV: {:.3} (mean={:.0})", current_cv, mean(&lengths_f64));
        eprintln!("  Merged CV:  {:.3} (mean={:.0})", merged_cv, mean(&merged_lengths));
    }

    // Merge if CV improves significantly (at least 50% reduction)
    if merged_cv >= current_cv * 0.5 {
        if verbose {
            eprintln!("  CV improvement not significant enough, skipping merge");
        }
        return (decomposition.to_vec(), false);
    }

    if verbose {
        eprintln!(
            "  ✓ Merging A-B pattern (CV improved by {:.0}%)",
            (1.0 - merged_cv / current_cv) * 100.0
        );
    }

    // Perform the merge
    let original_sequence: String = decomposition.join("");
    let mut new_decomposition: Vec<String> = Vec::new();
    let mut i = 0;

    while i < decomposition.len() {
        let current = &decomposition[i];

        if i < decomposition.len() - 1
            && current.starts_with(cut_seq)
            && decomposition[i + 1].starts_with(cut_seq)
        {
            let curr_len = current.len();
            let next_len = decomposition[i + 1].len();

            let curr_is_a = (curr_len as i64 - smaller_len as i64).abs() < (smaller_len as i64 / 5);
            let curr_is_b = (curr_len as i64 - larger_len as i64).abs() < (larger_len as i64 / 5);
            let next_is_a = (next_len as i64 - smaller_len as i64).abs() < (smaller_len as i64 / 5);
            let next_is_b = (next_len as i64 - larger_len as i64).abs() < (larger_len as i64 / 5);

            if (curr_is_a && next_is_b) || (curr_is_b && next_is_a) {
                let merged = format!("{}{}", current, decomposition[i + 1]);
                new_decomposition.push(merged);
                i += 2;
                continue;
            }
        }

        new_decomposition.push(current.clone());
        i += 1;
    }

    // Verify reconstruction
    let reconstructed: String = new_decomposition.join("");
    if reconstructed != original_sequence {
        if verbose {
            eprintln!("  ERROR: Reconstruction failed after A-B merge, reverting");
        }
        return (decomposition.to_vec(), false);
    }

    (new_decomposition, true)
}

/// Optimize monomer lengths by merging short frequent monomers
///
/// Post-processing optimization to merge short frequent monomers with adjacent longer ones.
/// Goal: minimize variance of monomer lengths.
pub fn optimize_monomer_lengths(
    decomposition: &[String],
    cut_seq: &str,
    verbose: bool,
) -> Vec<String> {
    if decomposition.len() < 3 {
        return decomposition.to_vec();
    }

    let original_sequence: String = decomposition.join("");

    // Skip if first fragment doesn't start with cut (it's a flank)
    let start_idx = if !decomposition[0].is_empty() && !decomposition[0].starts_with(cut_seq) {
        1
    } else {
        0
    };

    // Get monomer lengths
    let mut monomer_info: Vec<(usize, usize)> = Vec::new();
    for i in start_idx..decomposition.len() {
        if i == decomposition.len() - 1 {
            let avg_len: f64 = decomposition[start_idx..i]
                .iter()
                .map(|d| d.len() as f64)
                .sum::<f64>()
                / (i - start_idx).max(1) as f64;
            if (decomposition[i].len() as f64) < avg_len * 0.7 {
                continue;
            }
        }
        monomer_info.push((i, decomposition[i].len()));
    }

    if monomer_info.len() < 2 {
        return decomposition.to_vec();
    }

    // Count length frequencies
    let all_lengths: Vec<usize> = monomer_info.iter().map(|(_, len)| *len).collect();
    let length_counts = count_frequencies(&all_lengths);

    // Calculate initial variance
    let lengths_f64: Vec<f64> = all_lengths.iter().map(|&x| x as f64).collect();
    let initial_mean = mean(&lengths_f64);
    let initial_cv = cv(&lengths_f64);

    if verbose {
        eprintln!("Initial variance: CV={:.3} (mean={:.1})", initial_cv, initial_mean);
    }

    // Find frequently occurring short monomers
    let min_frequency = (monomer_info.len() * 15 / 100).max(2);
    let short_lengths: Vec<usize> = length_counts
        .iter()
        .filter(|(len, count)| **count >= min_frequency && (**len as f64) < initial_mean * 0.5)
        .map(|(len, _)| *len)
        .collect();

    if short_lengths.is_empty() {
        return decomposition.to_vec();
    }

    if verbose {
        eprintln!("Frequent short lengths: {:?}", short_lengths);
    }

    // Try merging short monomers with adjacent ones
    let mut working_decomposition = decomposition.to_vec();
    let mut new_decomposition: Vec<String> = Vec::new();
    let mut i = 0;

    while i < working_decomposition.len() {
        if i < working_decomposition.len() - 1 {
            let current = &working_decomposition[i];
            let next_frag = &working_decomposition[i + 1];

            let mut current_is_short = false;
            if current.starts_with(cut_seq) && next_frag.starts_with(cut_seq) {
                if short_lengths.contains(&current.len()) {
                    current_is_short = true;
                } else {
                    for &short_len in &short_lengths {
                        if (current.len() as i64 - short_len as i64).abs() < (short_len as i64 / 20) {
                            current_is_short = true;
                            break;
                        }
                    }
                }

                if current_is_short {
                    // Calculate variance after merge
                    let mut test_lengths: Vec<f64> = Vec::new();
                    for (j, frag) in working_decomposition.iter().enumerate() {
                        if j == i {
                            test_lengths.push((current.len() + next_frag.len()) as f64);
                        } else if j == i + 1 {
                            continue;
                        } else if frag.starts_with(cut_seq)
                            || (j > 0 && j < working_decomposition.len() - 1)
                        {
                            test_lengths.push(frag.len() as f64);
                        }
                    }

                    if !test_lengths.is_empty() {
                        let test_cv = cv(&test_lengths);
                        let length_ratio =
                            next_frag.len() as f64 / current.len().max(1) as f64;

                        if verbose {
                            eprintln!(
                                "  Testing merge {}+{}: CV {:.3} -> {:.3}",
                                current.len(),
                                next_frag.len(),
                                initial_cv,
                                test_cv
                            );
                        }

                        if test_cv < initial_cv * 0.98 || length_ratio > 3.0 {
                            let merged = format!("{}{}", current, next_frag);
                            new_decomposition.push(merged);
                            i += 2;
                            if verbose {
                                eprintln!("  ✓ Merged!");
                            }
                            continue;
                        }
                    }
                }
            }
        }

        new_decomposition.push(working_decomposition[i].clone());
        i += 1;
    }

    // Verify reconstruction
    let final_sequence: String = new_decomposition.join("");
    if final_sequence != original_sequence {
        eprintln!("ERROR: Sequence changed during merging! Reverting.");
        return decomposition.to_vec();
    }

    new_decomposition
}

/// Find anchor with mutations near expected position
pub fn find_mutant_anchor(
    sequence: &str,
    anchor: &str,
    expected_pos: usize,
    window: usize,
    max_dist: usize,
) -> Option<usize> {
    let anchor_len = anchor.len();
    let start = expected_pos.saturating_sub(window);
    let end = (expected_pos + window + anchor_len).min(sequence.len());

    let mut best_pos: Option<usize> = None;
    let mut best_dist = max_dist + 1;
    let mut best_pos_dist = usize::MAX;

    for pos in start..end {
        // Check different lengths due to possible indels
        for length in (anchor_len.saturating_sub(2))..=(anchor_len + 2) {
            if pos + length > sequence.len() {
                continue;
            }
            let candidate = &sequence[pos..pos + length];
            let dist = edit_distance(anchor, candidate);
            let pos_dist = (pos as i64 - expected_pos as i64).unsigned_abs() as usize;

            if dist < best_dist || (dist == best_dist && pos_dist < best_pos_dist) {
                best_dist = dist;
                best_pos = Some(pos);
                best_pos_dist = pos_dist;
            }
        }
    }

    if best_dist <= max_dist {
        best_pos
    } else {
        None
    }
}

/// Split monomers that are 2x, 3x longer than expected
///
/// This happens when anchor has mutation and we missed a cut point.
/// We search for mutant anchor near expected position and split there.
pub fn split_long_monomers(
    decomposition: &[String],
    cut_seq: &str,
    verbose: bool,
) -> Vec<String> {
    if decomposition.len() < 3 {
        return decomposition.to_vec();
    }

    let original_sequence: String = decomposition.join("");
    let mut current_decomposition = decomposition.to_vec();
    let mut total_splits = 0;

    for _iteration in 0..10 {
        // Calculate median length
        let lengths: Vec<usize> = current_decomposition.iter().map(|m| m.len()).collect();
        if lengths.is_empty() {
            break;
        }

        let mut sorted_lengths = lengths.clone();
        sorted_lengths.sort();
        let median_length = sorted_lengths[sorted_lengths.len() / 2];

        if verbose && _iteration == 0 {
            eprintln!(
                "Split long monomers: median={}bp, cut={}",
                median_length, cut_seq
            );
        }

        let mut result: Vec<String> = Vec::new();
        let mut splits_made = 0;

        // For short repeats, find the most common monomer
        let mono_counts = count_frequencies(&current_decomposition);
        let most_common_mono = if median_length < 23 {
            mono_counts
                .iter()
                .max_by_key(|(_, count)| *count)
                .map(|(mono, _)| mono.clone())
        } else {
            None
        };

        for (mono_idx, mono) in current_decomposition.iter().enumerate() {
            // Skip flanks
            if median_length < 23 {
                if let Some(ref common) = most_common_mono {
                    if mono == common {
                        result.push(mono.clone());
                        continue;
                    }
                }
            } else if !mono.starts_with(cut_seq) {
                result.push(mono.clone());
                continue;
            }

            let ratio = mono.len() as f64 / median_length as f64;

            if ratio < 1.3 {
                result.push(mono.clone());
                continue;
            }

            let n_expected = (ratio.round() as usize).max(2);

            if verbose {
                eprintln!(
                    "  #{}: Long monomer {}bp ({:.1}x), trying to split into {} parts",
                    mono_idx,
                    mono.len(),
                    ratio,
                    n_expected
                );
            }

            let mut parts = vec![mono.clone()];

            for split_num in 1..n_expected {
                let expected_pos = split_num * median_length;
                let current_part = parts.last().unwrap().clone();
                let offset: usize = parts[..parts.len() - 1].iter().map(|p| p.len()).sum();
                let rel_expected_pos = expected_pos.saturating_sub(offset);

                let min_pos = cut_seq.len() / 2;
                let max_pos = current_part.len().saturating_sub(cut_seq.len() / 2);

                if rel_expected_pos < min_pos || rel_expected_pos > max_pos {
                    continue;
                }

                let split_pos = if median_length < 23 {
                    if rel_expected_pos > 0 && rel_expected_pos < current_part.len() {
                        Some(rel_expected_pos)
                    } else {
                        None
                    }
                } else if ratio >= 1.7 {
                    let search_window = (median_length as f64 * 0.15).max(5.0) as usize;
                    find_mutant_anchor(&current_part, cut_seq, rel_expected_pos, search_window, 2)
                } else if current_part.starts_with(cut_seq) && rel_expected_pos >= cut_seq.len() {
                    Some(rel_expected_pos)
                } else {
                    None
                };

                if let Some(pos) = split_pos {
                    if pos > 0 {
                        let part1 = current_part[..pos].to_string();
                        let part2 = current_part[pos..].to_string();

                        // Check variance criterion
                        let dev_no_split =
                            (current_part.len() as i64 - median_length as i64).abs() as usize;
                        let dev_with_split = (part1.len() as i64 - median_length as i64).abs()
                            as usize
                            + (part2.len() as i64 - median_length as i64).abs() as usize;

                        if dev_with_split >= dev_no_split {
                            if verbose {
                                eprintln!(
                                    "    #{}: Skip split at {}: deviation {} -> {} (no improvement)",
                                    mono_idx, pos, dev_no_split, dev_with_split
                                );
                            }
                            continue;
                        }

                        *parts.last_mut().unwrap() = part1.clone();
                        parts.push(part2.clone());
                        splits_made += 1;

                        if verbose {
                            eprintln!(
                                "    #{}: Split at pos {}: {}bp + {}bp (dev {} -> {})",
                                mono_idx,
                                pos,
                                part1.len(),
                                part2.len(),
                                dev_no_split,
                                dev_with_split
                            );
                        }
                    }
                }
            }

            result.extend(parts);
        }

        total_splits += splits_made;
        current_decomposition = result;

        if splits_made == 0 {
            break;
        }

        if verbose {
            eprintln!("  Iteration {}: {} splits, continuing...", _iteration + 1, splits_made);
        }
    }

    // Verify reconstruction
    let final_sequence: String = current_decomposition.join("");
    if final_sequence != original_sequence {
        eprintln!("ERROR: Sequence changed during split_long_monomers!");
        return decomposition.to_vec();
    }

    if verbose && total_splits > 0 {
        eprintln!("  Total: {} splits", total_splits);
    }

    current_decomposition
}

/// Split monomers that consist of two identical halves
///
/// This catches cases like CCTAACCCTAAC which is actually CCTAAC + CCTAAC
pub fn split_duplicate_halves(
    decomposition: &[String],
    cut_seq: &str,
    verbose: bool,
) -> Vec<String> {
    if decomposition.len() < 3 {
        return decomposition.to_vec();
    }

    // Calculate median length
    let lengths: Vec<usize> = decomposition
        .iter()
        .filter(|m| m.starts_with(cut_seq))
        .map(|m| m.len())
        .collect();

    if lengths.is_empty() {
        return decomposition.to_vec();
    }

    let mut sorted_lengths = lengths.clone();
    sorted_lengths.sort();
    let median_length = sorted_lengths[sorted_lengths.len() / 2];

    // Only apply for short repeats
    if median_length > 24 {
        return decomposition.to_vec();
    }

    let original_sequence: String = decomposition.join("");
    let mut current_decomposition = decomposition.to_vec();
    let mut total_splits = 0;

    for _iteration in 0..10 {
        let mut result: Vec<String> = Vec::new();
        let mut splits_made = 0;

        for mono in &current_decomposition {
            let mono_len = mono.len();

            if mono_len % 2 != 0 {
                result.push(mono.clone());
                continue;
            }

            let half_len = mono_len / 2;
            let first_half = &mono[..half_len];
            let second_half = &mono[half_len..];

            if first_half == second_half {
                result.push(first_half.to_string());
                result.push(second_half.to_string());
                splits_made += 1;

                if verbose {
                    eprintln!(
                        "  Split duplicate halves: {}bp -> {}bp + {}bp ({})",
                        mono_len, half_len, half_len, first_half
                    );
                }
            } else {
                let dist = edit_distance(first_half, second_half);
                if dist <= 1 && half_len >= 4 {
                    result.push(first_half.to_string());
                    result.push(second_half.to_string());
                    splits_made += 1;

                    if verbose {
                        eprintln!(
                            "  Split near-duplicate halves (dist={}): {}bp -> {}bp + {}bp",
                            dist, mono_len, half_len, half_len
                        );
                    }
                } else {
                    result.push(mono.clone());
                }
            }
        }

        total_splits += splits_made;
        current_decomposition = result;

        if splits_made == 0 {
            break;
        }

        if verbose {
            eprintln!("  Duplicate halves iteration {}: {} splits", _iteration + 1, splits_made);
        }
    }

    // Verify reconstruction
    let final_sequence: String = current_decomposition.join("");
    if final_sequence != original_sequence {
        eprintln!("ERROR: Sequence changed during split_duplicate_halves!");
        return decomposition.to_vec();
    }

    if verbose && total_splits > 0 {
        eprintln!("  Total duplicate-half splits: {}", total_splits);
    }

    current_decomposition
}

/// Split monomers that start with the most popular monomer
///
/// For short repeats, if a monomer starts with the most common monomer sequence
/// but is longer, split it into (popular_monomer) + (remainder).
pub fn split_by_popular_prefix(
    decomposition: &[String],
    _cut_seq: &str,
    verbose: bool,
) -> Vec<String> {
    if decomposition.len() < 3 {
        return decomposition.to_vec();
    }

    // Find the most common monomer
    let mono_counts = count_frequencies(decomposition);
    if mono_counts.is_empty() {
        return decomposition.to_vec();
    }

    let (popular_mono, popular_count) = mono_counts
        .iter()
        .max_by_key(|(_, count)| *count)
        .map(|(mono, count)| (mono.clone(), *count))
        .unwrap();

    let popular_len = popular_mono.len();

    // Only apply for short repeats
    if popular_len > 20 {
        return decomposition.to_vec();
    }

    // Need significant majority
    if popular_count < decomposition.len() * 30 / 100 {
        return decomposition.to_vec();
    }

    let original_sequence: String = decomposition.join("");
    let mut result: Vec<String> = Vec::new();
    let mut splits_made = 0;

    for mono in decomposition {
        let mono_len = mono.len();

        if mono_len <= popular_len {
            result.push(mono.clone());
            continue;
        }

        if mono.starts_with(&popular_mono) {
            let remainder = &mono[popular_len..];
            let remainder_len = remainder.len();

            // Variance check
            let dev_before = (mono_len as i64 - popular_len as i64).abs() as usize;
            let dev_after = (remainder_len as i64 - popular_len as i64).abs() as usize;

            if dev_after < dev_before {
                result.push(popular_mono.clone());
                result.push(remainder.to_string());
                splits_made += 1;

                if verbose {
                    eprintln!(
                        "  Split by prefix: {}bp -> {}bp + {}bp (dev {} -> {})",
                        mono_len, popular_len, remainder_len, dev_before, dev_after
                    );
                }
            } else {
                result.push(mono.clone());
                if verbose {
                    eprintln!(
                        "  Skip prefix split: {}bp (dev {} -> {}, no improvement)",
                        mono_len, dev_before, dev_after
                    );
                }
            }
        } else {
            result.push(mono.clone());
        }
    }

    // Verify reconstruction
    let final_sequence: String = result.join("");
    if final_sequence != original_sequence {
        eprintln!("ERROR: Sequence changed during split_by_popular_prefix!");
        return decomposition.to_vec();
    }

    if verbose && splits_made > 0 {
        eprintln!(
            "  Split {} monomers by popular prefix ({})",
            splits_made, popular_mono
        );
    }

    result
}

/// Apply all post-processing heuristics to a decomposition
pub fn apply_all_heuristics(
    decomposition: &[String],
    cut_seq: &str,
    verbose: bool,
) -> Vec<String> {
    let mut result = decomposition.to_vec();

    // 1. Detect and fix A-B alternating pattern
    let (new_result, _ab_merged) = detect_and_fix_ab_pattern(&result, cut_seq, verbose);
    result = new_result;

    // 2. Optimize monomer lengths (merge short frequent monomers)
    result = optimize_monomer_lengths(&result, cut_seq, verbose);

    // 3. Split duplicate halves (e.g., CCTAACCCTAAC -> CCTAAC + CCTAAC)
    result = split_duplicate_halves(&result, cut_seq, verbose);

    // 4. Split long monomers where anchor has mutation
    result = split_long_monomers(&result, cut_seq, verbose);

    // 5. Split monomers that start with popular monomer
    result = split_by_popular_prefix(&result, cut_seq, verbose);

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_split_duplicate_halves() {
        // Function requires at least 3 elements
        let decomposition = vec![
            "ACGTACGT".to_string(), // Should split into ACGT + ACGT
            "ACGT".to_string(),     // Should not split
            "ACGT".to_string(),     // Keep as is
        ];

        let result = split_duplicate_halves(&decomposition, "ACGT", false);

        assert_eq!(result.len(), 4);
        assert_eq!(result[0], "ACGT");
        assert_eq!(result[1], "ACGT");
        assert_eq!(result[2], "ACGT");
        assert_eq!(result[3], "ACGT");
    }

    #[test]
    fn test_find_mutant_anchor() {
        let sequence = "ACGTACGTACGT";
        let anchor = "ACGT";

        // Should find at expected position
        let pos = find_mutant_anchor(sequence, anchor, 4, 10, 2);
        assert!(pos.is_some());
        assert_eq!(pos.unwrap(), 4);

        // Should find with 1 mismatch
        let sequence_mut = "ACGTACCTACGT"; // ACGT -> ACCT at pos 4
        let pos = find_mutant_anchor(sequence_mut, anchor, 4, 10, 2);
        assert!(pos.is_some());
    }

    #[test]
    fn test_cv_calculation() {
        let values = vec![100.0, 100.0, 100.0];
        assert_eq!(cv(&values), 0.0);

        let values = vec![100.0, 200.0];
        let c = cv(&values);
        assert!(c > 0.0);
    }
}
