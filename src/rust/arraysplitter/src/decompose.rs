//! Core decomposition algorithm for satellite DNA arrays
//!
//! Implements the main decomposition pipeline that:
//! 1. Builds frequency suffix tree to find candidate cut sequences (all nucleotides)
//! 2. Uses anchor graph to select optimal cuts (primary method)
//! 3. Falls back to FS-tree decomposition if anchor graph fails
//! 4. Applies post-processing heuristics to refine decomposition

use std::collections::HashMap;

use crate::sequence::find_gcd_of_list;
use crate::fs_tree::iterate_hints;
use crate::anchor_graph::AnchorGraphDecomposer;
use crate::autocorrelation;
use crate::anchor_by_period;
use crate::multiplet_split;

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
    /// Alternative anchors with same period (for redundant anchor refinement)
    pub alternative_anchors: Vec<String>,
    /// Period detected by autocorrelation (new pipeline)
    pub autocorr_period: Option<usize>,
    /// Autocorrelation value at detected period
    pub autocorr_value: Option<f64>,
    /// Source label per monomer: "anchor", "split_2x", "anchor_graph", "fs_tree", etc.
    pub monomer_sources: Vec<String>,
}

/// Count frequencies of items
fn count_frequencies<T: std::hash::Hash + Eq + Clone>(items: &[T]) -> HashMap<T, usize> {
    let mut counts = HashMap::new();
    for item in items {
        *counts.entry(item.clone()).or_insert(0) += 1;
    }
    counts
}

/// Get canonical orientation of a sequence (matching Python exactly)
///
/// Returns True if sequence is already canonical (A > T, or A == T and C > G)
/// When A == T and C == G, returns True (consider it canonical)
pub fn is_canonical_orientation(sequence: &str) -> bool {
    let sequence_upper = sequence.to_uppercase();
    let a_count = sequence_upper.chars().filter(|&c| c == 'A').count();
    let t_count = sequence_upper.chars().filter(|&c| c == 'T').count();
    let c_count = sequence_upper.chars().filter(|&c| c == 'C').count();
    let g_count = sequence_upper.chars().filter(|&c| c == 'G').count();

    // Primary criterion: A > T
    if a_count != t_count {
        return a_count > t_count;
    }

    // Secondary criterion: C > G (when A == T)
    // Python: return c_count > g_count
    c_count > g_count
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

/// Candidate for cut sequence selection (matching Python structure)
#[derive(Debug, Clone)]
pub struct CutCandidate {
    pub cut: String,
    pub mode_period: usize,
    pub base_score: f64,
    pub adjusted_score: f64,
    pub fragmentation: f64,
    pub period_gcd: usize,
    pub period_distribution: HashMap<usize, usize>,
    pub num_segments: usize,
    pub num_parts: usize,
    pub empty_ratio: f64,
}

/// Compute best cut sequence from hints (matching Python compute_cuts exactly)
///
/// Key features:
/// - Handles empty_ratio for perfect/near-perfect repeats
/// - Calculates fragmentation penalty
/// - Checks for fundamental period using GCD
/// - Sorts by (-adjusted_score, num_segments, mode_period)
pub fn compute_cuts(
    array: &str,
    hints: &[(usize, String, usize)],
    score_threshold: f64,
    fragmentation_threshold: f64,
) -> (String, f64, usize) {
    let mut candidates: Vec<CutCandidate> = Vec::new();

    // Calculate metrics for each hint (matching Python exactly)
    for (_, cut_sequence, _) in hints {
        let parts: Vec<&str> = array.split(cut_sequence.as_str()).collect();
        let mut periods: Vec<usize> = Vec::new();
        let mut non_empty_periods: Vec<usize> = Vec::new();

        for (i, part) in parts.iter().enumerate() {
            // Handle edge cases for first/last parts
            // Python: if i < len(parts) - 1 or part:
            if i < parts.len() - 1 || !part.is_empty() {
                let period = part.len() + cut_sequence.len();
                periods.push(period);
                // Track non-empty parts separately
                if !part.is_empty() {
                    non_empty_periods.push(period);
                }
            }
        }

        if periods.is_empty() {
            continue;
        }

        // Determine if this is a perfect/near-perfect repeat
        let empty_ratio = (periods.len() - non_empty_periods.len()) as f64 / periods.len() as f64;

        let (period_counts, mode_period, mode_count, total_segments) = if empty_ratio >= 0.8 {
            // 80% or more empty parts = perfect/near-perfect repeat
            // Use the cut sequence length as the period
            let mut counts = HashMap::new();
            counts.insert(cut_sequence.len(), periods.len());
            (counts, cut_sequence.len(), periods.len(), periods.len())
        } else if !non_empty_periods.is_empty() {
            // Use only non-empty parts for period calculation
            let counts = count_frequencies(&non_empty_periods);
            let (mode, count) = counts
                .iter()
                .max_by_key(|(_, count)| *count)
                .map(|(k, v)| (*k, *v))
                .unwrap_or((0, 0));
            (counts, mode, count, non_empty_periods.len())
        } else {
            // Fallback: use all periods
            let counts = count_frequencies(&periods);
            let (mode, count) = counts
                .iter()
                .max_by_key(|(_, count)| *count)
                .map(|(k, v)| (*k, *v))
                .unwrap_or((0, 0));
            (counts, mode, count, periods.len())
        };

        // Base score (uniformity)
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
            period_distribution: period_counts,
            num_segments: total_segments,
            num_parts: parts.len(),
            empty_ratio,
        });
    }

    if candidates.is_empty() {
        return (String::new(), 0.0, array.len());
    }

    // Group candidates by score
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
            // Check if most periods are multiples of GCD
            let multiples: usize = c.period_distribution
                .keys()
                .filter(|&&p| p % c.period_gcd == 0)
                .count();
            if multiples >= ((c.num_segments as f64 * 0.8) as usize).max(1) {
                fundamental_candidates.push(c);
            }
        }
    }

    // If we found fundamental periods, use those
    if !fundamental_candidates.is_empty() {
        let best = fundamental_candidates
            .iter()
            .min_by_key(|c| c.period_gcd)
            .unwrap();
        return (best.cut.clone(), best.base_score, best.period_gcd);
    }

    // Otherwise, penalize fragmentation and choose best
    // Sort by adjusted score, then by number of segments (fewer is better), then by period (smaller is better)
    let mut similar: Vec<CutCandidate> = similar_candidates.iter().map(|c| (*c).clone()).collect();
    similar.sort_by(|a, b| {
        // Python: similar_candidates.sort(key=lambda x: (-x['adjusted_score'], x['num_segments'], x['mode_period']))
        let score_cmp = b
            .adjusted_score
            .partial_cmp(&a.adjusted_score)
            .unwrap_or(std::cmp::Ordering::Equal);
        if score_cmp == std::cmp::Ordering::Equal {
            let seg_cmp = a.num_segments.cmp(&b.num_segments);
            if seg_cmp == std::cmp::Ordering::Equal {
                a.mode_period.cmp(&b.mode_period)
            } else {
                seg_cmp
            }
        } else {
            score_cmp
        }
    });

    let best = &similar[0];
    (best.cut.clone(), best.base_score, best.mode_period)
}

/// Decompose array using the cut sequence (matching Python decompose_array_iter1)
pub fn decompose_array_iter1(
    array: &str,
    best_cut_seq: &str,
    _best_period: usize,
    verbose: bool,
) -> (Vec<String>, HashMap<String, usize>) {
    let mut repeats2count: HashMap<String, usize> = HashMap::new();
    let mut decomposition: Vec<String> = Vec::new();

    if best_cut_seq.is_empty() || !array.contains(best_cut_seq) {
        // No cut sequence or not found, return whole array
        decomposition.push(array.to_string());
        *repeats2count.entry(array.to_string()).or_insert(0) += 1;
        return (decomposition, repeats2count);
    }

    // Split by cut sequence
    let parts: Vec<&str> = array.split(best_cut_seq).collect();

    // First part is a special case (flank)
    if !parts[0].is_empty() {
        // First fragment (before first cut) - this is a flank
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

/// Get candidates from hints for anchor graph (matching Python get_candidates_from_hints)
fn get_candidates_from_hints(array: &str, hints: &[(usize, String, usize)]) -> Vec<CutCandidate> {
    let mut candidates: Vec<CutCandidate> = Vec::new();

    for (_l, cut_sequence, _n) in hints {
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
            period_distribution: period_counts,
            num_segments: total_segments,
            num_parts: parts.len(),
            empty_ratio,
        });
    }

    candidates
}

/// Get all hints from all nucleotides (matching Python get_all_hints_for_graph)
fn get_all_hints_for_graph(array: &str, depth: usize, cutoff: usize) -> Vec<(usize, String, usize)> {
    let mut all_hints: Vec<(usize, String, usize)> = Vec::new();

    for nucleotide in "ACTG".chars() {
        // Get positions of this nucleotide
        let positions: Vec<usize> = array
            .char_indices()
            .filter(|(_, c)| c.to_ascii_uppercase() == nucleotide)
            .map(|(i, _)| i)
            .collect();

        if positions.len() <= cutoff {
            continue;
        }

        // Get hints from fs_tree
        let hints = iterate_hints(array, nucleotide, cutoff, depth);
        for hint in hints {
            all_hints.push((hint.length, hint.pattern, hint.frequency));
        }
    }

    // Remove duplicates, keeping highest frequency
    let mut unique: HashMap<(usize, String), (usize, String, usize)> = HashMap::new();
    for (length, anchor, freq) in all_hints {
        let key = (length, anchor.clone());
        if let Some(existing) = unique.get(&key) {
            if freq > existing.2 {
                unique.insert(key, (length, anchor, freq));
            }
        } else {
            unique.insert(key, (length, anchor, freq));
        }
    }

    unique.into_values().collect()
}

/// Decompose array using anchor graph (primary method, matching Python)
fn decompose_with_anchor_graph(
    array: &str,
    hints: &[(usize, String, usize)],
    verbose: bool,
) -> Option<(Vec<String>, String, usize, f64, Vec<String>)> {
    // Build candidates from hints
    let candidates = get_candidates_from_hints(array, hints);

    if candidates.is_empty() {
        return None;
    }

    // Get anchor strings from candidates (max length 11)
    // Sort by adjusted_score descending first
    let mut sorted_candidates = candidates;
    sorted_candidates.sort_by(|a, b| {
        b.adjusted_score
            .partial_cmp(&a.adjusted_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let anchors: Vec<String> = sorted_candidates
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

    // Verify reconstruction
    let reconstructed: String = result.monomers.join("");
    if reconstructed != array {
        if verbose {
            eprintln!("  Graph reconstruction failed");
        }
        return None;
    }

    Some((result.monomers, result.cut_sequence, result.period, result.cv, result.alternative_anchors))
}

/// Main decomposition function (matching Python decompose_array exactly)
///
/// Decomposes a satellite DNA array into monomers using:
/// 1. Gets hints from ALL nucleotides (ACTG)
/// 2. Tries anchor graph decomposition first (primary method)
/// 3. Falls back to FS-tree decomposition if anchor graph fails
/// 4. Post-processing heuristics are applied separately
pub fn decompose_array(
    array: &str,
    depth: usize,
    cutoff: Option<usize>,
    verbose: bool,
) -> Decomposition {
    // Set cutoff based on array size if not provided (matching Python)
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

    // Step 1-3. Get hints from all nucleotides instead of just the most frequent
    let mut all_hints: Vec<(usize, String, usize)> = Vec::new();
    let mut hint_sources: HashMap<(usize, String), (usize, char)> = HashMap::new();

    for nucleotide in "ACTG".chars() {
        // Get positions of this nucleotide
        let positions: Vec<usize> = array
            .char_indices()
            .filter(|(_, c)| c.to_ascii_uppercase() == nucleotide)
            .map(|(i, _)| i)
            .collect();

        if positions.len() <= cutoff {
            continue;
        }

        // Get hints from fs_tree for this nucleotide
        let hints = iterate_hints(array, nucleotide, cutoff, depth);
        for hint in hints {
            all_hints.push((hint.length, hint.pattern.clone(), hint.frequency));
            let hint_key = (hint.length, hint.pattern.clone());
            if let Some(existing) = hint_sources.get(&hint_key) {
                if hint.frequency > existing.0 {
                    hint_sources.insert(hint_key, (hint.frequency, nucleotide));
                }
            } else {
                hint_sources.insert(hint_key, (hint.frequency, nucleotide));
            }
        }
    }

    // Remove duplicates, keeping the one with highest frequency
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

    // Step 4. PRIMARY: Try anchor graph decomposition first
    // For large sequences, get more hints with smaller cutoff
    let graph_hints = if array.len() > 10_000 {
        let extra_hints = get_all_hints_for_graph(array, 100, 3);
        if extra_hints.len() > hints.len() {
            extra_hints
        } else {
            hints.clone()
        }
    } else {
        hints.clone()
    };

    if let Some((monomers, cut_seq, period, cv, alt_anchors)) =
        decompose_with_anchor_graph(array, &graph_hints, verbose)
    {
        if verbose {
            eprintln!(
                "  Using anchor graph decomposition: {} monomers, period={}, CV={:.3}",
                monomers.len(),
                period,
                cv
            );
        }

        // Calculate best_cut_score as fraction of monomers starting with cut
        let monomers_with_cut = monomers.iter().filter(|m| m.starts_with(&cut_seq)).count();
        let score = monomers_with_cut as f64 / monomers.len() as f64;

        return Decomposition {
            monomers,
            cut_sequence: cut_seq,
            score,
            period,
            was_reversed: false,
            cv,
            alternative_anchors: alt_anchors,
            autocorr_period: None,
            autocorr_value: None,
            monomer_sources: Vec::new(),
        };
    }

    // Step 5. FALLBACK: Use FS-tree decomposition if anchor graph fails
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
        alternative_anchors: Vec::new(),
        autocorr_period: None,
        autocorr_value: None,
        monomer_sources: Vec::new(),
    }
}

/// Decompose array using predefined cut sequences (matching Python decompose_array_with_cuts)
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

        // Calculate score (same logic as compute_cuts)
        let period_counts = count_frequencies(&periods);
        let (mode_period, mode_count) = period_counts
            .iter()
            .max_by_key(|(_, count)| *count)
            .map(|(k, v)| (*k, *v))
            .unwrap_or((0, 0));

        let total_segments = periods.len();
        let base_score = mode_count as f64 / total_segments as f64;

        // Fragmentation penalty
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
        // No cuts found, return whole array
        return Decomposition {
            monomers: vec![array.to_string()],
            cut_sequence: String::new(),
            score: 0.0,
            period: array.len(),
            was_reversed: false,
            cv: 0.0,
            alternative_anchors: Vec::new(),
            autocorr_period: None,
            autocorr_value: None,
            monomer_sources: Vec::new(),
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
        alternative_anchors: Vec::new(),
        autocorr_period: None,
        autocorr_value: None,
        monomer_sources: Vec::new(),
    }
}

/// Decompose array using the new autocorrelation-based pipeline.
///
/// Pipeline:
/// 1. Detect orientation (canonical form)
/// 2. Autocorrelation → find period
/// 3. Anchor by period → find anchor k-mer and positions
/// 4. Multiplet split → split segments spanning multiple periods
/// 5. Extract monomer sequences from boundaries
pub fn decompose_array_autocorr(
    array: &str,
    verbose: bool,
    max_period_cap: usize,
) -> Decomposition {
    let seq = array.as_bytes();
    let min_period = 5;
    let max_period = (array.len() / 3).min(max_period_cap); // need >=3 copies, cap at user-supplied --max-period

    // Step 1: Autocorrelation period detection
    // Use find_period (not find_period_refined) to preserve HOR structure.
    // find_period returns the period with maximum autocorrelation excess,
    // while find_period_refined would "refine" down to smaller sub-periods.
    let period_result = autocorrelation::find_period(seq, min_period, max_period);

    let (period, autocorr_value) = match period_result {
        Some((p, ac, _excess)) => (p, ac),
        None => {
            if verbose {
                eprintln!("  Autocorrelation: no periodicity detected");
            }
            // Return whole array as single "monomer"
            return Decomposition {
                monomers: vec![array.to_string()],
                cut_sequence: String::new(),
                score: 0.0,
                period: array.len(),
                was_reversed: false,
                cv: 0.0,
                alternative_anchors: Vec::new(),
                autocorr_period: None,
                autocorr_value: None,
                monomer_sources: vec!["no_period".to_string()],
            };
        }
    };

    if verbose {
        eprintln!("  Autocorrelation: period={}, autocorr={:.4}", period, autocorr_value);
    }

    // Step 2: Find anchor k-mer using the period constraint
    let anchor_result = anchor_by_period::find_anchor_by_period_with_fallback(seq, period);

    let (anchor_kmer, anchor_positions) = match anchor_result {
        Some(result) => {
            if verbose {
                eprintln!(
                    "  Anchor: k-mer={}, k={}, positions={}, score={:.3} (uniq={:.3}, reg={:.3})",
                    result.kmer, result.k, result.positions.len(),
                    result.score, result.uniqueness, result.regularity
                );
            }
            (result.kmer, result.positions)
        }
        None => {
            if verbose {
                eprintln!("  Anchor: no suitable anchor found, using even split");
            }
            // Fallback: even split based on period
            let n_monomers = array.len() / period;
            let mut monomers = Vec::new();
            let mut sources = Vec::new();
            for i in 0..n_monomers {
                let start = i * period;
                let end = if i == n_monomers - 1 { array.len() } else { (i + 1) * period };
                monomers.push(array[start..end].to_string());
                sources.push("even_split".to_string());
            }

            let cv = compute_cv_from_lengths(&monomers.iter().map(|m| m.len()).collect::<Vec<_>>());

            return Decomposition {
                monomers,
                cut_sequence: String::new(),
                score: 0.0,
                period,
                was_reversed: false,
                cv,
                alternative_anchors: Vec::new(),
                autocorr_period: Some(period),
                autocorr_value: Some(autocorr_value),
                monomer_sources: sources,
            };
        }
    };

    // Step 3: Multiplet splitting
    let boundaries = multiplet_split::split_multiplets(seq, &anchor_positions, period);

    if verbose {
        let n_anchor = boundaries.iter().filter(|b| b.source == "anchor").count();
        let n_split: usize = boundaries.iter().filter(|b| b.source.starts_with("split_")).count();
        let n_flank = boundaries.iter().filter(|b| b.source.contains("flank")).count();
        eprintln!(
            "  Multiplet split: {} monomers ({} anchor, {} split, {} flank)",
            boundaries.len(), n_anchor, n_split, n_flank
        );
    }

    // Step 4: Extract monomer sequences and build result
    let mut monomers: Vec<String> = Vec::new();
    let mut monomer_sources: Vec<String> = Vec::new();

    for boundary in &boundaries {
        let monomer_seq = &array[boundary.start..boundary.end];
        monomers.push(monomer_seq.to_string());
        monomer_sources.push(boundary.source.clone());
    }

    // Verify reconstruction
    let reconstructed: String = monomers.join("");
    if reconstructed != array {
        if verbose {
            eprintln!(
                "  WARNING: Reconstruction mismatch! array={} vs reconstructed={}",
                array.len(), reconstructed.len()
            );
        }
    }

    // Compute CV from monomer lengths (excluding flanks)
    let monomer_lengths: Vec<usize> = boundaries
        .iter()
        .zip(monomers.iter())
        .filter(|(b, _)| !b.source.contains("flank"))
        .map(|(_, m)| m.len())
        .collect();

    let cv = compute_cv_from_lengths(&monomer_lengths);

    // Score: fraction of boundaries from anchors (high = good)
    let anchor_count = boundaries.iter().filter(|b| b.source == "anchor").count();
    let score = if boundaries.is_empty() {
        0.0
    } else {
        anchor_count as f64 / boundaries.len() as f64
    };

    Decomposition {
        monomers,
        cut_sequence: anchor_kmer,
        score,
        period,
        was_reversed: false,
        cv,
        alternative_anchors: Vec::new(),
        autocorr_period: Some(period),
        autocorr_value: Some(autocorr_value),
        monomer_sources,
    }
}

/// Compute coefficient of variation from a list of lengths
fn compute_cv_from_lengths(lengths: &[usize]) -> f64 {
    if lengths.len() < 2 {
        return 0.0;
    }

    let mean: f64 = lengths.iter().sum::<usize>() as f64 / lengths.len() as f64;
    if mean == 0.0 {
        return f64::INFINITY;
    }

    let variance: f64 = lengths
        .iter()
        .map(|&x| (x as f64 - mean).powi(2))
        .sum::<f64>()
        / lengths.len() as f64;

    variance.sqrt() / mean
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
        assert!(is_canonical_orientation("ATCC"));
        // A == T, C < G should not be canonical
        assert!(!is_canonical_orientation("ATGG"));
        // A == T, C == G - Python returns c_count > g_count which is False
        assert!(!is_canonical_orientation("ATCG"));
    }

    #[test]
    fn test_decompose_simple() {
        let array = "ACGTACGTACGTACGT";
        let result = decompose_array(array, 100, None, false);

        // Should decompose into monomers
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

    #[test]
    fn test_compute_cuts_with_empty_ratio() {
        // Test perfect repeat (all parts empty except start/end)
        let array = "AAAAAAAAAA";
        let hints = vec![(1, "A".to_string(), 10)];
        let (cut, score, period) = compute_cuts(&array, &hints, 0.05, 0.5);

        assert_eq!(cut, "A");
        // For perfect repeat, period should be cut length
        assert!(period >= 1);
    }

    // === Edge case tests for autocorr pipeline ===

    #[test]
    fn test_autocorr_no_tandem_repeat() {
        // Random non-periodic sequence (LCG-generated)
        let mut seq = Vec::with_capacity(2000);
        let mut state: u64 = 42;
        for _ in 0..2000 {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let nucleotide = match (state >> 33) % 4 {
                0 => 'A', 1 => 'C', 2 => 'G', _ => 'T',
            };
            seq.push(nucleotide);
        }
        let array: String = seq.into_iter().collect();

        let result = decompose_array_autocorr(&array, true, 5000);

        // Should return single "monomer" = whole array, source = "no_period"
        assert_eq!(result.monomers.len(), 1, "Non-periodic should give 1 monomer");
        assert_eq!(result.monomer_sources[0], "no_period");
        assert_eq!(result.monomers[0], array);
    }

    #[test]
    fn test_autocorr_two_different_repeats() {
        // TTAGGG×100 followed by AATGG×120 (mixed telomere + satIII)
        let mut array = String::new();
        for _ in 0..100 {
            array.push_str("TTAGGG");
        }
        for _ in 0..120 {
            array.push_str("AATGG");
        }
        // Total: 600 + 600 = 1200bp

        let result = decompose_array_autocorr(&array, true, 5000);

        assert!(!result.monomers.is_empty(), "Should produce at least 1 monomer");
        eprintln!("Two repeats: period={}, n_monomers={}, sources={:?}",
            result.period, result.monomers.len(),
            result.monomer_sources.iter().take(5).collect::<Vec<_>>());

        // For chimeric array: period=30 = LCM(5,6) is the correct shared period.
        // Both halves have AC=1.0 at d=30, while d=5 and d=6 only work for half.
        // This is the expected limitation: chimeric arrays need pre-splitting.
        assert!(result.period == 5 || result.period == 6 || result.period == 30,
            "Period should be 5, 6, or 30 (LCM), got {}", result.period);
    }

    #[test]
    fn test_autocorr_alternating_ab_pattern() {
        // A-block B-block A-block B-block (A=ACGTTAGC ×10, B=TTGGCCAA ×10)
        let unit_a = "ACGTTAGC"; // period 8
        let unit_b = "TTGGCCAA"; // different 8bp unit
        let mut array = String::new();
        for _ in 0..3 {
            for _ in 0..10 { array.push_str(unit_a); }
            for _ in 0..10 { array.push_str(unit_b); }
        }
        // 3 × (80 + 80) = 480bp, super-period = 160bp

        let result = decompose_array_autocorr(&array, true, 5000);

        assert!(!result.monomers.is_empty());
        eprintln!("A-B-A-B: period={}, n_monomers={}, cv={:.3}",
            result.period, result.monomers.len(), result.cv);

        // Should detect either period=8 (fundamental) or period=160 (super-period)
        // Both are valid interpretations
        assert!(result.period == 8 || result.period == 160 || result.period == 80,
            "Period should be 8, 80, or 160, got {}", result.period);
    }

    #[test]
    fn test_autocorr_abc_a_pattern() {
        // A B C A: different units concatenated, A repeats
        // This tests what happens with non-uniform structure
        let unit_a = "ACGTACGTACGTACGT"; // 16bp
        let unit_b = "TTTTGGGGCCCCAAAA"; // 16bp different
        let unit_c = "AACCGGTTAACCGGTT"; // 16bp different
        let mut array = String::new();
        // A×20 + B×20 + C×20 + A×20 = 1280bp
        for _ in 0..20 { array.push_str(unit_a); }
        for _ in 0..20 { array.push_str(unit_b); }
        for _ in 0..20 { array.push_str(unit_c); }
        for _ in 0..20 { array.push_str(unit_a); }

        let result = decompose_array_autocorr(&array, true, 5000);

        assert!(!result.monomers.is_empty());
        eprintln!("A-B-C-A: period={}, n_monomers={}, cv={:.3}",
            result.period, result.monomers.len(), result.cv);

        // Each unit has internal period of 4 (ACGT, TTTT|GGGG|CCCC|AAAA, AACCGGTT)
        // or period 8 (AACCGGTT). The overall array might detect 4 or 16 or something else.
        // Main check: doesn't crash, produces reasonable output
        assert!(result.period <= 320, "Period should be reasonable, got {}", result.period);
    }

    #[test]
    fn test_autocorr_repeat_with_gap() {
        // A-repeat, then unique gap, then A-repeat again
        // ACGTTAGC×20 + random_50bp + ACGTTAGC×20
        let unit = "ACGTTAGC";
        let mut array = String::new();
        for _ in 0..20 { array.push_str(unit); }
        // 50bp "gap" (non-periodic)
        array.push_str("TGCAATGCAATTTGGGCCCAAAGGGCCCAAATTTGGGCCCAAAGGGCCCTT");
        for _ in 0..20 { array.push_str(unit); }
        // Total: 160 + 50 + 160 = 370bp

        let result = decompose_array_autocorr(&array, true, 5000);

        assert!(!result.monomers.is_empty());
        eprintln!("Repeat-gap-repeat: period={}, n_monomers={}, cv={:.3}",
            result.period, result.monomers.len(), result.cv);

        // Should detect period=8 (dominant repeat), gap will be a multiplet
        assert_eq!(result.period, 8,
            "Should detect period 8 through the gap, got {}", result.period);
    }

    #[test]
    fn test_autocorr_telomere_simple() {
        // Simple TTAGGG×200
        let mut array = String::new();
        for _ in 0..200 {
            array.push_str("TTAGGG");
        }

        let result = decompose_array_autocorr(&array, true, 5000);

        assert_eq!(result.period, 6, "Telomere period should be 6, got {}", result.period);
        // Should produce ~200 monomers (+ possible short flanks from anchor offset)
        let n = result.monomers.len();
        assert!(n >= 195 && n <= 210,
            "Expected ~200 monomers, got {}", n);
        // Most monomers should be 6bp; first/last may be flanks (shorter)
        let full_monomers: Vec<_> = result.monomers.iter()
            .filter(|m| m.len() >= 5 && m.len() <= 7)
            .collect();
        assert!(full_monomers.len() >= 195,
            "Expected >=195 full-size monomers, got {}", full_monomers.len());
        // Total sequence length should be preserved
        let total_len: usize = result.monomers.iter().map(|m| m.len()).sum();
        assert_eq!(total_len, 1200, "Total length should be 1200bp");
    }
}
