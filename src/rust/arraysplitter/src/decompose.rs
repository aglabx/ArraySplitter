//! Core decomposition algorithm for satellite DNA arrays
//!
//! Implements the main decomposition pipeline that:
//! 1. Builds frequency suffix tree to find candidate cut sequences
//! 2. Uses anchor graph to select optimal cuts
//! 3. Decomposes arrays into monomers

use std::collections::HashMap;

use crate::sequence::find_gcd_of_list;
use crate::fs_tree::FsTree;
use crate::anchor_graph::AnchorGraphDecomposer;

/// Result of array decomposition
#[derive(Debug, Clone)]
pub struct Decomposition {
    /// List of monomer sequences
    pub monomers: Vec<String>,
    /// The cut sequence used
    pub cut_sequence: String,
    /// Score of the decomposition
    pub score: f64,
    /// Estimated period (monomer length)
    pub period: usize,
    /// Whether the array was reverse-complemented
    pub was_reversed: bool,
    /// Coefficient of variation
    pub cv: f64,
}

/// Count frequencies of items
fn count_frequencies<T: std::hash::Hash + Eq + Clone>(items: &[T]) -> HashMap<T, usize> {
    let mut counts = HashMap::new();
    for item in items {
        *counts.entry(item.clone()).or_insert(0) += 1;
    }
    counts
}

/// Get canonical orientation of a sequence
///
/// Returns True if sequence is already canonical (A > T, or A == T and C >= G)
/// When counts are equal, we consider it canonical to avoid unnecessary reversal.
pub fn is_canonical_orientation(sequence: &str) -> bool {
    let a_count = sequence.chars().filter(|&c| c == 'A' || c == 'a').count();
    let t_count = sequence.chars().filter(|&c| c == 'T' || c == 't').count();
    let c_count = sequence.chars().filter(|&c| c == 'C' || c == 'c').count();
    let g_count = sequence.chars().filter(|&c| c == 'G' || c == 'g').count();

    if a_count != t_count {
        return a_count > t_count;
    }
    // When A == T, use C >= G (treat tie as canonical)
    c_count >= g_count
}

/// Rotate monomers so they start with the cut sequence
pub fn rotate_monomers_to_cut(decomposition: &[String], cut_sequence: &str) -> Vec<String> {
    decomposition
        .iter()
        .map(|monomer| {
            if monomer.starts_with(cut_sequence) {
                monomer.clone()
            } else if let Some(pos) = monomer.find(cut_sequence) {
                format!("{}{}", &monomer[pos..], &monomer[..pos])
            } else {
                monomer.clone()
            }
        })
        .collect()
}

/// Get the most frequent nucleotide in a sequence
pub fn get_top1_nucleotide(array: &str) -> char {
    let mut counts: HashMap<char, usize> = HashMap::new();
    for n in "ACTG".chars() {
        counts.insert(n, array.chars().filter(|&c| c == n).count());
    }

    counts
        .into_iter()
        .max_by_key(|(_, count)| *count)
        .map(|(c, _)| c)
        .unwrap_or('A')
}

/// Candidate for cut sequence selection
#[derive(Debug, Clone)]
pub struct CutCandidate {
    pub cut: String,
    pub mode_period: usize,
    pub base_score: f64,
    pub adjusted_score: f64,
    pub fragmentation: f64,
    pub period_gcd: usize,
    pub num_segments: usize,
    pub num_parts: usize,
    pub empty_ratio: f64,
}

/// Compute best cut sequence from hints
pub fn compute_cuts(
    array: &str,
    hints: &[(usize, String, usize)],
    score_threshold: f64,
    fragmentation_threshold: f64,
) -> (String, f64, usize) {
    let mut candidates: Vec<CutCandidate> = Vec::new();

    for (_, cut_sequence, _) in hints {
        let parts: Vec<&str> = array.split(cut_sequence.as_str()).collect();
        let mut periods: Vec<usize> = Vec::new();
        let mut non_empty_periods: Vec<usize> = Vec::new();

        for (i, part) in parts.iter().enumerate() {
            if i < parts.len() - 1 || !part.is_empty() {
                let period = part.len() + cut_sequence.len();
                periods.push(period);
                if !part.is_empty() {
                    non_empty_periods.push(period);
                }
            }
        }

        if periods.is_empty() {
            continue;
        }

        let empty_ratio = (periods.len() - non_empty_periods.len()) as f64 / periods.len() as f64;

        let (period_counts, mode_period, mode_count, total_segments) = if empty_ratio >= 0.8 {
            let mut counts = HashMap::new();
            counts.insert(cut_sequence.len(), periods.len());
            (counts, cut_sequence.len(), periods.len(), periods.len())
        } else if !non_empty_periods.is_empty() {
            let counts = count_frequencies(&non_empty_periods);
            let (mode, count) = counts
                .iter()
                .max_by_key(|(_, count)| *count)
                .map(|(k, v)| (*k, *v))
                .unwrap_or((0, 0));
            (counts, mode, count, non_empty_periods.len())
        } else {
            let counts = count_frequencies(&periods);
            let (mode, count) = counts
                .iter()
                .max_by_key(|(_, count)| *count)
                .map(|(k, v)| (*k, *v))
                .unwrap_or((0, 0));
            (counts, mode, count, periods.len())
        };

        let base_score = mode_count as f64 / total_segments as f64;

        // Fragmentation penalty
        let short_threshold = (mode_period as f64 * fragmentation_threshold) as usize;
        let short_fragments = periods.iter().filter(|&&p| p < short_threshold).count();
        let fragmentation = short_fragments as f64 / total_segments as f64;

        // Check for periodicity/divisibility
        let unique_periods: Vec<usize> = period_counts.keys().copied().collect();
        let period_gcd = if unique_periods.len() > 1 {
            find_gcd_of_list(&unique_periods)
        } else {
            mode_period
        };

        candidates.push(CutCandidate {
            cut: cut_sequence.clone(),
            mode_period,
            base_score,
            adjusted_score: base_score * (1.0 - fragmentation * 0.5),
            fragmentation,
            period_gcd,
            num_segments: total_segments,
            num_parts: parts.len(),
            empty_ratio,
        });
    }

    if candidates.is_empty() {
        return (String::new(), 0.0, array.len());
    }

    // Sort by score
    candidates.sort_by(|a, b| {
        b.base_score
            .partial_cmp(&a.base_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let best_base_score = candidates[0].base_score;

    // Get all candidates within threshold
    let similar_candidates: Vec<&CutCandidate> = candidates
        .iter()
        .filter(|c| c.base_score >= best_base_score - score_threshold)
        .collect();

    // Check for fundamental period (GCD > 1)
    let mut fundamental_candidates: Vec<&CutCandidate> = Vec::new();
    for c in &similar_candidates {
        if c.period_gcd > 1 && c.period_gcd < c.mode_period {
            fundamental_candidates.push(c);
        }
    }

    if !fundamental_candidates.is_empty() {
        let best = fundamental_candidates
            .iter()
            .min_by_key(|c| c.period_gcd)
            .unwrap();
        return (best.cut.clone(), best.base_score, best.period_gcd);
    }

    // Sort by adjusted score, then by period
    let mut similar: Vec<CutCandidate> = similar_candidates.iter().map(|c| (*c).clone()).collect();
    similar.sort_by(|a, b| {
        let score_cmp = b
            .adjusted_score
            .partial_cmp(&a.adjusted_score)
            .unwrap_or(std::cmp::Ordering::Equal);
        if score_cmp == std::cmp::Ordering::Equal {
            a.mode_period.cmp(&b.mode_period)
        } else {
            score_cmp
        }
    });

    let best = &similar[0];
    (best.cut.clone(), best.base_score, best.mode_period)
}

/// Decompose array using the cut sequence (iteration 1)
pub fn decompose_array_iter1(
    array: &str,
    best_cut_seq: &str,
    _best_period: usize,
    verbose: bool,
) -> (Vec<String>, HashMap<String, usize>) {
    let mut repeats2count: HashMap<String, usize> = HashMap::new();
    let mut decomposition: Vec<String> = Vec::new();

    if best_cut_seq.is_empty() || !array.contains(best_cut_seq) {
        decomposition.push(array.to_string());
        *repeats2count.entry(array.to_string()).or_insert(0) += 1;
        return (decomposition, repeats2count);
    }

    // Split by cut sequence
    let parts: Vec<&str> = array.split(best_cut_seq).collect();

    // First part is a flank
    if !parts[0].is_empty() {
        decomposition.push(parts[0].to_string());
        *repeats2count.entry(parts[0].to_string()).or_insert(0) += 1;
        if verbose {
            eprintln!("Flank: {}bp", parts[0].len());
        }
    }

    // Process all other parts: cut + part = monomer
    for i in 1..parts.len() {
        let monomer = format!("{}{}", best_cut_seq, parts[i]);
        *repeats2count.entry(monomer.clone()).or_insert(0) += 1;
        decomposition.push(monomer.clone());
        if verbose {
            eprintln!(
                "Monomer {}: {}bp (cut {}bp + part {}bp)",
                i,
                monomer.len(),
                best_cut_seq.len(),
                parts[i].len()
            );
        }
    }

    // Verify reconstruction
    let reconstructed: String = decomposition.join("");
    if reconstructed != array {
        eprintln!(
            "WARNING: Reconstruction mismatch! {} != {}",
            array.len(),
            reconstructed.len()
        );
    }

    (decomposition, repeats2count)
}

/// Get candidates from hints for anchor graph
fn get_candidates_from_hints(array: &str, hints: &[(usize, String, usize)]) -> Vec<CutCandidate> {
    let mut candidates: Vec<CutCandidate> = Vec::new();

    for (l, cut_sequence, _n) in hints {
        let parts: Vec<&str> = array.split(cut_sequence.as_str()).collect();
        let mut periods: Vec<usize> = Vec::new();
        let mut non_empty_periods: Vec<usize> = Vec::new();

        for (i, part) in parts.iter().enumerate() {
            if i < parts.len() - 1 || !part.is_empty() {
                let period = part.len() + cut_sequence.len();
                periods.push(period);
                if !part.is_empty() {
                    non_empty_periods.push(period);
                }
            }
        }

        if periods.is_empty() {
            continue;
        }

        let empty_ratio = (periods.len() - non_empty_periods.len()) as f64 / periods.len() as f64;

        let (period_counts, mode_period, mode_count, total_segments) = if empty_ratio >= 0.8 {
            let mut counts = HashMap::new();
            counts.insert(cut_sequence.len(), periods.len());
            (counts, cut_sequence.len(), periods.len(), periods.len())
        } else if !non_empty_periods.is_empty() {
            let counts = count_frequencies(&non_empty_periods);
            let (mode, count) = counts
                .iter()
                .max_by_key(|(_, count)| *count)
                .map(|(k, v)| (*k, *v))
                .unwrap_or((0, 0));
            (counts, mode, count, non_empty_periods.len())
        } else {
            let counts = count_frequencies(&periods);
            let (mode, count) = counts
                .iter()
                .max_by_key(|(_, count)| *count)
                .map(|(k, v)| (*k, *v))
                .unwrap_or((0, 0));
            (counts, mode, count, periods.len())
        };

        let base_score = mode_count as f64 / total_segments as f64;
        let short_threshold = (mode_period as f64 * 0.5) as usize;
        let short_fragments = periods.iter().filter(|&&p| p < short_threshold).count();
        let fragmentation = short_fragments as f64 / total_segments as f64;
        let adjusted_score = base_score * (1.0 - fragmentation * 0.5);

        candidates.push(CutCandidate {
            cut: cut_sequence.clone(),
            mode_period,
            base_score,
            adjusted_score,
            fragmentation,
            period_gcd: mode_period,
            num_segments: total_segments,
            num_parts: parts.len(),
            empty_ratio,
        });
    }

    candidates
}

/// Decompose array using anchor graph
fn decompose_with_anchor_graph(
    array: &str,
    hints: &[(usize, String, usize)],
    verbose: bool,
) -> Option<(Vec<String>, String, usize, f64)> {
    // Build candidates from hints
    let candidates = get_candidates_from_hints(array, hints);

    if candidates.is_empty() {
        return None;
    }

    // Get anchor strings from candidates (max length 11)
    let anchors: Vec<String> = candidates
        .iter()
        .filter(|c| c.cut.len() <= 11)
        .map(|c| c.cut.clone())
        .take(15)
        .collect();

    if anchors.is_empty() {
        return None;
    }

    // Build anchor graph
    let mut decomposer = AnchorGraphDecomposer::new();
    decomposer.build_from_candidates(array, &anchors, verbose);

    // Decompose
    let result = decomposer.decompose(verbose)?;

    // Calculate CV
    let cut_seq = result.cut_sequence.clone();
    let lengths: Vec<usize> = result
        .monomers
        .iter()
        .filter(|m| m.starts_with(&cut_seq))
        .map(|m| m.len())
        .collect();

    if lengths.len() < 2 {
        return None;
    }

    let mean_len: f64 = lengths.iter().sum::<usize>() as f64 / lengths.len() as f64;
    let variance: f64 =
        lengths.iter().map(|&x| (x as f64 - mean_len).powi(2)).sum::<f64>() / lengths.len() as f64;
    let cv = if mean_len > 0.0 {
        variance.sqrt() / mean_len
    } else {
        f64::INFINITY
    };

    // Verify reconstruction
    let reconstructed: String = result.monomers.join("");
    if reconstructed != array {
        if verbose {
            eprintln!("  Graph reconstruction failed");
        }
        return None;
    }

    Some((result.monomers, cut_seq, result.period, cv))
}

/// Main decomposition function
///
/// Decomposes a satellite DNA array into monomers using:
/// 1. Frequency suffix tree to find candidate cut sequences
/// 2. Anchor graph to select optimal cut
/// 3. Post-processing heuristics to refine decomposition
pub fn decompose_array(
    array: &str,
    depth: usize,
    cutoff: Option<usize>,
    verbose: bool,
) -> Decomposition {
    // Set cutoff based on array size if not provided
    let cutoff = cutoff.unwrap_or_else(|| {
        if array.len() > 1_000_000 {
            1000
        } else if array.len() > 100_000 {
            250
        } else if array.len() > 10_000 {
            10
        } else {
            3
        }
    });

    // Collect hints from all nucleotides
    let mut all_hints: Vec<(usize, String, usize)> = Vec::new();

    for nucleotide in "ACTG".chars() {
        let positions: Vec<usize> = array
            .char_indices()
            .filter(|(_, c)| *c == nucleotide)
            .map(|(i, _)| i)
            .collect();

        if positions.len() <= cutoff {
            continue;
        }

        // Build fs_tree for this nucleotide
        let fs_tree = FsTree::new(array, 3, depth.min(100), cutoff);

        // Get hints from this fs_tree
        for hint in fs_tree.get_hints_starting_with(nucleotide) {
            all_hints.push((hint.pattern.len(), hint.pattern, hint.frequency));
        }
    }

    // Remove duplicates, keeping highest frequency
    let mut unique_hints: HashMap<(usize, String), (usize, String, usize)> = HashMap::new();
    for (length, sequence, freq) in all_hints {
        let key = (length, sequence.clone());
        if let Some(existing) = unique_hints.get(&key) {
            if freq > existing.2 {
                unique_hints.insert(key, (length, sequence, freq));
            }
        } else {
            unique_hints.insert(key, (length, sequence, freq));
        }
    }

    let hints: Vec<(usize, String, usize)> = unique_hints.into_values().collect();

    // Try anchor graph decomposition first
    if let Some((monomers, cut_seq, period, cv)) =
        decompose_with_anchor_graph(array, &hints, verbose)
    {
        if verbose {
            eprintln!(
                "  Using anchor graph decomposition: {} monomers, period={}, CV={:.3}",
                monomers.len(),
                period,
                cv
            );
        }

        let monomers_with_cut = monomers.iter().filter(|m| m.starts_with(&cut_seq)).count();
        let score = monomers_with_cut as f64 / monomers.len() as f64;

        return Decomposition {
            monomers,
            cut_sequence: cut_seq,
            score,
            period,
            was_reversed: false,
            cv,
        };
    }

    // Fallback to FS-tree decomposition
    if verbose {
        eprintln!("  Anchor graph failed, falling back to FS-tree decomposition");
    }

    let (best_cut_seq, best_cut_score, best_period) =
        compute_cuts(array, &hints, 0.05, 0.5);

    let (monomers, _repeats2count) =
        decompose_array_iter1(array, &best_cut_seq, best_period, verbose);

    // Calculate CV
    let lengths: Vec<usize> = monomers
        .iter()
        .filter(|m| m.starts_with(&best_cut_seq))
        .map(|m| m.len())
        .collect();

    let cv = if lengths.len() > 1 {
        let mean_len: f64 = lengths.iter().sum::<usize>() as f64 / lengths.len() as f64;
        let variance: f64 = lengths
            .iter()
            .map(|&x| (x as f64 - mean_len).powi(2))
            .sum::<f64>()
            / lengths.len() as f64;
        if mean_len > 0.0 {
            variance.sqrt() / mean_len
        } else {
            f64::INFINITY
        }
    } else {
        0.0
    };

    Decomposition {
        monomers,
        cut_sequence: best_cut_seq,
        score: best_cut_score,
        period: best_period,
        was_reversed: false,
        cv,
    }
}

/// Decompose array using predefined cut sequences
pub fn decompose_array_with_cuts(
    array: &str,
    cut_sequences: &[String],
    verbose: bool,
) -> Decomposition {
    if cut_sequences.is_empty() {
        panic!("No cut sequences provided");
    }

    let mut best_result: Option<(String, f64, usize, usize)> = None;
    let mut best_score: f64 = -1.0;

    for cut_seq in cut_sequences {
        if !array.contains(cut_seq.as_str()) {
            continue;
        }

        let parts: Vec<&str> = array.split(cut_seq.as_str()).collect();
        let mut periods: Vec<usize> = Vec::new();

        for (i, part) in parts.iter().enumerate() {
            if i < parts.len() - 1 || !part.is_empty() {
                periods.push(part.len() + cut_seq.len());
            }
        }

        if periods.is_empty() {
            continue;
        }

        let period_counts = count_frequencies(&periods);
        let (mode_period, mode_count) = period_counts
            .iter()
            .max_by_key(|(_, count)| *count)
            .map(|(k, v)| (*k, *v))
            .unwrap_or((0, 0));

        let total_segments = periods.len();
        let base_score = mode_count as f64 / total_segments as f64;

        let short_threshold = (mode_period as f64 * 0.5) as usize;
        let short_fragments = periods.iter().filter(|&&p| p < short_threshold).count();
        let fragmentation = short_fragments as f64 / total_segments as f64;
        let adjusted_score = base_score * (1.0 - fragmentation * 0.5);

        if adjusted_score > best_score {
            best_score = adjusted_score;
            best_result = Some((cut_seq.clone(), base_score, mode_period, parts.len()));
        }
    }

    if best_result.is_none() {
        return Decomposition {
            monomers: vec![array.to_string()],
            cut_sequence: String::new(),
            score: 0.0,
            period: array.len(),
            was_reversed: false,
            cv: 0.0,
        };
    }

    let (cut_seq, score, period, _) = best_result.unwrap();

    let (monomers, _repeats2count) = decompose_array_iter1(array, &cut_seq, period, verbose);

    // Calculate CV
    let lengths: Vec<usize> = monomers
        .iter()
        .filter(|m| m.starts_with(&cut_seq))
        .map(|m| m.len())
        .collect();

    let cv = if lengths.len() > 1 {
        let mean_len: f64 = lengths.iter().sum::<usize>() as f64 / lengths.len() as f64;
        let variance: f64 = lengths
            .iter()
            .map(|&x| (x as f64 - mean_len).powi(2))
            .sum::<f64>()
            / lengths.len() as f64;
        if mean_len > 0.0 {
            variance.sqrt() / mean_len
        } else {
            f64::INFINITY
        }
    } else {
        0.0
    };

    Decomposition {
        monomers,
        cut_sequence: cut_seq,
        score,
        period,
        was_reversed: false,
        cv,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_canonical_orientation() {
        // A > T should be canonical
        assert!(is_canonical_orientation("AAAA"));
        // T > A should not be canonical
        assert!(!is_canonical_orientation("TTTT"));
        // A == T, C > G should be canonical
        assert!(is_canonical_orientation("ACGT"));
    }

    #[test]
    fn test_decompose_simple() {
        let array = "ACGTACGTACGTACGT";
        let result = decompose_array(array, 100, None, false);

        // Should decompose into 4bp monomers
        assert!(result.monomers.len() >= 2);

        // Verify reconstruction
        let reconstructed: String = result.monomers.join("");
        assert_eq!(reconstructed, array);
    }

    #[test]
    fn test_decompose_with_predefined_cuts() {
        let array = "ACGTACGTACGTACGT";
        let cuts = vec!["ACGT".to_string()];
        let result = decompose_array_with_cuts(array, &cuts, false);

        assert_eq!(result.cut_sequence, "ACGT");
        assert_eq!(result.period, 4);

        // Verify reconstruction
        let reconstructed: String = result.monomers.join("");
        assert_eq!(reconstructed, array);
    }
}
