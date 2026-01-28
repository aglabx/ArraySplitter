//! Recursive HOR (Higher Order Repeat) decomposition
//!
//! Decomposes HOR monomers into their constituent base-level monomers.
//! Uses autocorrelation to detect periodicity within each HOR and recursively
//! splits until reaching base monomers (no detectable periodicity).

use crate::autocorrelation::{find_period_refined, random_expectation, autocorrelation};
use crate::anchor_by_period::find_anchor_by_period_with_fallback;
use crate::multiplet_split::split_multiplets;

/// A base-level monomer after recursive decomposition
#[derive(Debug, Clone)]
pub struct BaseMonomer {
    /// The DNA sequence of this monomer
    pub sequence: String,
    /// Index of the parent HOR from primary decomposition
    pub hor_idx: usize,
    /// Global index within the entire array (0, 1, 2, ...)
    pub global_idx: usize,
    /// Index within the parent HOR (0, 1, 2, ...)
    pub sub_idx: usize,
    /// Recursion depth (1 = direct child of HOR, 2 = grandchild, etc.)
    pub level: usize,
    /// Detected period at this level (0 if this is a base monomer)
    pub period: usize,
    /// Autocorrelation value at the detected period
    pub autocorr: f64,
    /// Source of this monomer: "recursive_anchor", "recursive_split", "base", "flank"
    pub source: String,
    /// Edit distance to template (consensus of submonomers within same parent)
    pub ed_tmpl: Option<usize>,
    /// Edit distance to previous submonomer
    pub ed_prev: Option<usize>,
    /// Edit distance to next submonomer
    pub ed_next: Option<usize>,
    /// Normalized edit distance (ed_tmpl / length)
    pub ed_per_bp: f64,
    /// Coefficient of variation for length within this parent group
    pub cv: f64,
}

/// Result of recursive HOR decomposition
#[derive(Debug, Clone)]
pub struct RecursiveResult {
    /// All base-level monomers from the recursive decomposition
    pub base_monomers: Vec<BaseMonomer>,
    /// Maximum recursion depth reached
    pub max_depth: usize,
    /// Total expected count of base monomers (sum across all HORs)
    pub n_expected: usize,
    /// Median period of base monomers
    pub median_period: usize,
    /// Mean autocorrelation of base monomers
    pub mean_autocorr: f64,
    /// Consensus sequence of base monomers
    pub consensus_seq: String,
    /// IUPAC ambiguity codes
    pub iupac_str: String,
    /// Quality string (digit 0-9 per position)
    pub quality_str: String,
}

/// Default minimum submonomer length (bp)
pub const DEFAULT_MIN_SUBMONOMER_LEN: usize = 5;

/// Default autocorrelation threshold for further decomposition
pub const DEFAULT_AUTOCORR_THRESHOLD: f64 = 0.5;

/// Maximum length for edit distance computation to avoid quadratic blowup
const MAX_ED_LEN: usize = 10000;

/// Calculate edit distance between two sequences
fn edit_distance(s1: &str, s2: &str) -> usize {
    let len1 = s1.len();
    let len2 = s2.len();

    // Skip ED for very long sequences
    if len1 > MAX_ED_LEN || len2 > MAX_ED_LEN {
        return usize::MAX / 2;
    }

    if len1 == 0 { return len2; }
    if len2 == 0 { return len1; }

    let s1_chars: Vec<char> = s1.chars().collect();
    let s2_chars: Vec<char> = s2.chars().collect();

    let mut prev_row: Vec<usize> = (0..=len2).collect();
    let mut curr_row: Vec<usize> = vec![0; len2 + 1];

    for i in 1..=len1 {
        curr_row[0] = i;
        for j in 1..=len2 {
            let cost = if s1_chars[i - 1] == s2_chars[j - 1] { 0 } else { 1 };
            curr_row[j] = (prev_row[j] + 1)
                .min(curr_row[j - 1] + 1)
                .min(prev_row[j - 1] + cost);
        }
        std::mem::swap(&mut prev_row, &mut curr_row);
    }

    prev_row[len2]
}

/// Compute consensus sequence from a list of sequences (simple majority vote)
/// Used internally for template building
fn compute_consensus(sequences: &[&str]) -> String {
    if sequences.is_empty() {
        return String::new();
    }

    // Use median length as consensus length
    let mut lengths: Vec<usize> = sequences.iter().map(|s| s.len()).collect();
    lengths.sort();
    let consensus_len = lengths[lengths.len() / 2];

    let mut consensus = String::with_capacity(consensus_len);

    for pos in 0..consensus_len {
        let mut counts = [0u32; 4]; // A, C, G, T

        for seq in sequences {
            let bytes = seq.as_bytes();
            if pos < bytes.len() {
                match bytes[pos] {
                    b'A' | b'a' => counts[0] += 1,
                    b'C' | b'c' => counts[1] += 1,
                    b'G' | b'g' => counts[2] += 1,
                    b'T' | b't' => counts[3] += 1,
                    _ => {},
                }
            }
        }

        let (best_idx, _) = counts.iter().enumerate()
            .max_by_key(|(_, &c)| c)
            .unwrap();

        let base = match best_idx {
            0 => 'A', 1 => 'C', 2 => 'G', 3 => 'T', _ => 'N',
        };
        consensus.push(base);
    }

    consensus
}

/// Compute full consensus with IUPAC codes and quality digits
/// Returns (consensus, iupac, quality)
pub fn compute_consensus_full(sequences: &[&str]) -> (String, String, String) {
    if sequences.is_empty() {
        return (String::new(), String::new(), String::new());
    }

    // Use median length as consensus length
    let mut lengths: Vec<usize> = sequences.iter().map(|s| s.len()).collect();
    lengths.sort();
    let consensus_len = lengths[lengths.len() / 2];

    let mut consensus = String::with_capacity(consensus_len);
    let mut iupac = String::with_capacity(consensus_len);
    let mut quality = String::with_capacity(consensus_len);

    for pos in 0..consensus_len {
        let mut counts = [0u32; 4]; // A, C, G, T
        let mut total = 0u32;

        for seq in sequences {
            let bytes = seq.as_bytes();
            if pos < bytes.len() {
                match bytes[pos] {
                    b'A' | b'a' => { counts[0] += 1; total += 1; },
                    b'C' | b'c' => { counts[1] += 1; total += 1; },
                    b'G' | b'g' => { counts[2] += 1; total += 1; },
                    b'T' | b't' => { counts[3] += 1; total += 1; },
                    _ => {},
                }
            }
        }

        if total == 0 {
            consensus.push('N');
            iupac.push('N');
            quality.push('0');
            continue;
        }

        // Find most common base
        let (best_idx, &best_count) = counts.iter().enumerate()
            .max_by_key(|(_, &c)| c)
            .unwrap();

        let support = best_count as f64 / total as f64;

        // Consensus: uppercase
        let base = match best_idx {
            0 => 'A', 1 => 'C', 2 => 'G', 3 => 'T', _ => 'N',
        };
        consensus.push(base);

        // Quality: digit 0-9
        let digit = (support * 10.0).min(9.99) as u8;
        quality.push((b'0' + digit) as char);

        // IUPAC: bases with frequency >= 20%
        let threshold = (total as f64 * 0.2) as u32;
        let present: Vec<usize> = (0..4)
            .filter(|&i| counts[i] >= threshold.max(1))
            .collect();

        let iupac_char = match present.as_slice() {
            [0] => 'A', [1] => 'C', [2] => 'G', [3] => 'T',
            [0, 2] => 'R', // A+G purine
            [1, 3] => 'Y', // C+T pyrimidine
            [1, 2] => 'S', // G+C strong
            [0, 3] => 'W', // A+T weak
            [2, 3] => 'K', // G+T keto
            [0, 1] => 'M', // A+C amino
            [1, 2, 3] => 'B', // not A
            [0, 2, 3] => 'D', // not C
            [0, 1, 3] => 'H', // not G
            [0, 1, 2] => 'V', // not T
            [0, 1, 2, 3] => 'N',
            _ => 'N',
        };
        iupac.push(iupac_char);
    }

    (consensus, iupac, quality)
}

/// Compute edit distance metrics for base monomers
/// Groups monomers by parent HOR, computes template (consensus), and calculates
/// ed_tmpl, ed_prev, ed_next, ed_per_bp, and cv for each monomer
pub fn compute_submonomer_metrics(base_monomers: &mut Vec<BaseMonomer>) {
    if base_monomers.is_empty() {
        return;
    }

    // Group monomers by hor_idx
    let mut groups: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::new();
    for (i, mono) in base_monomers.iter().enumerate() {
        groups.entry(mono.hor_idx).or_default().push(i);
    }

    // Process each group
    for (_hor_idx, indices) in &groups {
        if indices.is_empty() {
            continue;
        }

        // Build template (consensus) for this group
        let sequences: Vec<&str> = indices.iter()
            .map(|&i| base_monomers[i].sequence.as_str())
            .collect();
        let template = compute_consensus(&sequences);

        // Compute lengths for CV calculation
        let lengths: Vec<f64> = indices.iter()
            .map(|&i| base_monomers[i].sequence.len() as f64)
            .collect();
        let mean_len = if !lengths.is_empty() {
            lengths.iter().sum::<f64>() / lengths.len() as f64
        } else {
            0.0
        };
        let cv = if lengths.len() > 1 && mean_len > 0.0 {
            let variance: f64 = lengths.iter()
                .map(|&x| (x - mean_len).powi(2))
                .sum::<f64>() / lengths.len() as f64;
            variance.sqrt() / mean_len
        } else {
            0.0
        };

        // Compute ed_tmpl for each monomer in group
        for &idx in indices {
            let ed = edit_distance(&template, &base_monomers[idx].sequence);
            let len = base_monomers[idx].sequence.len();
            base_monomers[idx].ed_tmpl = Some(ed);
            base_monomers[idx].ed_per_bp = if len > 0 { ed as f64 / len as f64 } else { 0.0 };
            base_monomers[idx].cv = cv;
        }
    }

    // Compute ed_prev and ed_next between adjacent monomers (globally, not per group)
    let n = base_monomers.len();
    for i in 0..n {
        if i > 0 {
            let ed = edit_distance(&base_monomers[i - 1].sequence, &base_monomers[i].sequence);
            base_monomers[i].ed_prev = Some(ed);
        }
        if i + 1 < n {
            let ed = edit_distance(&base_monomers[i].sequence, &base_monomers[i + 1].sequence);
            base_monomers[i].ed_next = Some(ed);
        }
    }
}

/// Decompose HOR monomers recursively into base-level monomers
///
/// # Arguments
/// * `hor_monomers` - List of HOR monomer sequences from primary decomposition
/// * `min_submonomer_len` - Minimum length to attempt further decomposition (default: 5bp)
/// * `autocorr_threshold` - Autocorrelation threshold above which to decompose (default: 0.5)
///
/// # Returns
/// A `RecursiveResult` containing all base-level monomers and max recursion depth
pub fn decompose_hors_to_base(
    hor_monomers: &[String],
    min_submonomer_len: usize,
    autocorr_threshold: f64,
) -> RecursiveResult {
    let mut all_base_monomers: Vec<BaseMonomer> = Vec::new();
    let mut max_depth: usize = 0;

    for (hor_idx, hor_seq) in hor_monomers.iter().enumerate() {
        let base_monomers = decompose_single_hor(
            hor_seq,
            hor_idx,
            min_submonomer_len,
            autocorr_threshold,
            1, // Start at level 1
            &mut max_depth,
        );
        all_base_monomers.extend(base_monomers);
    }

    // Assign global indices
    for (i, mono) in all_base_monomers.iter_mut().enumerate() {
        mono.global_idx = i;
    }

    // Compute edit distance metrics
    compute_submonomer_metrics(&mut all_base_monomers);

    // Compute summary statistics
    let n_expected = all_base_monomers.len();

    // Median period
    let mut periods: Vec<usize> = all_base_monomers.iter()
        .filter(|m| m.period > 0)
        .map(|m| m.period)
        .collect();
    periods.sort();
    let median_period = if !periods.is_empty() {
        periods[periods.len() / 2]
    } else {
        // Use median length if no periods
        let mut lengths: Vec<usize> = all_base_monomers.iter().map(|m| m.sequence.len()).collect();
        lengths.sort();
        if !lengths.is_empty() { lengths[lengths.len() / 2] } else { 0 }
    };

    // Mean autocorrelation
    let autocorr_values: Vec<f64> = all_base_monomers.iter()
        .filter(|m| m.autocorr > 0.0)
        .map(|m| m.autocorr)
        .collect();
    let mean_autocorr = if !autocorr_values.is_empty() {
        autocorr_values.iter().sum::<f64>() / autocorr_values.len() as f64
    } else {
        0.0
    };

    // Compute consensus for base monomers
    let sequences: Vec<&str> = all_base_monomers.iter()
        .map(|m| m.sequence.as_str())
        .collect();
    let (consensus_seq, iupac_str, quality_str) = if sequences.len() >= 2 {
        compute_consensus_full(&sequences)
    } else {
        (String::new(), String::new(), String::new())
    };

    RecursiveResult {
        base_monomers: all_base_monomers,
        max_depth,
        n_expected,
        median_period,
        mean_autocorr,
        consensus_seq,
        iupac_str,
        quality_str,
    }
}

/// Decompose a single HOR monomer recursively
fn decompose_single_hor(
    hor_seq: &str,
    hor_idx: usize,
    min_len: usize,
    autocorr_threshold: f64,
    current_level: usize,
    max_depth: &mut usize,
) -> Vec<BaseMonomer> {
    // Update max depth
    if current_level > *max_depth {
        *max_depth = current_level;
    }

    // Too short to decompose further
    if hor_seq.len() < min_len * 2 {
        return vec![BaseMonomer {
            sequence: hor_seq.to_string(),
            hor_idx,
            global_idx: 0, // Will be assigned later
            sub_idx: 0,
            level: current_level,
            period: 0,
            autocorr: 0.0,
            source: "base".to_string(),
            ed_tmpl: None,
            ed_prev: None,
            ed_next: None,
            ed_per_bp: 0.0,
            cv: 0.0,
        }];
    }

    let seq_bytes = hor_seq.as_bytes();
    let min_period = min_len;
    let max_period = hor_seq.len() / 2; // Need at least 2 copies

    // Check for periodicity
    let period_result = find_period_refined(seq_bytes, min_period, max_period);

    match period_result {
        Some((period, autocorr, _excess)) => {
            // Check if autocorrelation exceeds threshold
            if autocorr <= autocorr_threshold {
                // No strong periodicity - this is a base monomer
                return vec![BaseMonomer {
                    sequence: hor_seq.to_string(),
                    hor_idx,
                    global_idx: 0, // Will be assigned later
                    sub_idx: 0,
                    level: current_level,
                    period: 0,
                    autocorr,
                    source: "base".to_string(),
                    ed_tmpl: None,
                    ed_prev: None,
                    ed_next: None,
                    ed_per_bp: 0.0,
                    cv: 0.0,
                }];
            }

            // Strong periodicity detected - decompose further
            let sub_monomers = decompose_hor_by_period(hor_seq, period, autocorr);

            // Recursively decompose each sub-monomer
            let mut result: Vec<BaseMonomer> = Vec::new();
            for (sub_idx, (sub_seq, sub_source)) in sub_monomers.iter().enumerate() {
                if sub_seq.len() < min_len {
                    // Too short - keep as is
                    result.push(BaseMonomer {
                        sequence: sub_seq.clone(),
                        hor_idx,
                        global_idx: 0, // Will be assigned later
                        sub_idx,
                        level: current_level,
                        period,
                        autocorr,
                        source: sub_source.clone(),
                        ed_tmpl: None,
                        ed_prev: None,
                        ed_next: None,
                        ed_per_bp: 0.0,
                        cv: 0.0,
                    });
                } else {
                    // Try to decompose further
                    let sub_bytes = sub_seq.as_bytes();
                    let sub_min_period = min_len;
                    let sub_max_period = sub_seq.len() / 2;

                    if sub_max_period >= sub_min_period {
                        if let Some((_, sub_autocorr, _)) = find_period_refined(sub_bytes, sub_min_period, sub_max_period) {
                            if sub_autocorr > autocorr_threshold {
                                // Recurse
                                let deeper = decompose_single_hor(
                                    sub_seq,
                                    hor_idx,
                                    min_len,
                                    autocorr_threshold,
                                    current_level + 1,
                                    max_depth,
                                );
                                // Adjust sub_idx for the deeper monomers
                                for mut m in deeper {
                                    m.sub_idx = sub_idx * 1000 + m.sub_idx; // Encode hierarchy
                                    result.push(m);
                                }
                                continue;
                            }
                        }
                    }

                    // No further periodicity - this is a base monomer
                    let final_autocorr = if sub_max_period >= sub_min_period {
                        find_period_refined(sub_bytes, sub_min_period, sub_max_period)
                            .map(|(_, ac, _)| ac)
                            .unwrap_or(0.0)
                    } else {
                        0.0
                    };

                    result.push(BaseMonomer {
                        sequence: sub_seq.clone(),
                        hor_idx,
                        global_idx: 0, // Will be assigned later
                        sub_idx,
                        level: current_level,
                        period,
                        autocorr: final_autocorr,
                        source: sub_source.clone(),
                        ed_tmpl: None,
                        ed_prev: None,
                        ed_next: None,
                        ed_per_bp: 0.0,
                        cv: 0.0,
                    });
                }
            }
            result
        }
        None => {
            // No periodicity detected - return as base monomer
            vec![BaseMonomer {
                sequence: hor_seq.to_string(),
                hor_idx,
                global_idx: 0, // Will be assigned later
                sub_idx: 0,
                level: current_level,
                period: 0,
                autocorr: compute_max_autocorr(seq_bytes, min_len),
                source: "base".to_string(),
                ed_tmpl: None,
                ed_prev: None,
                ed_next: None,
                ed_per_bp: 0.0,
                cv: 0.0,
            }]
        }
    }
}

/// Decompose a HOR sequence using already-detected period (no re-detection)
fn decompose_hor_by_period(hor_seq: &str, period: usize, _autocorr: f64) -> Vec<(String, String)> {
    let seq_len = hor_seq.len();
    let seq_bytes = hor_seq.as_bytes();

    // For very short sequences or period too large, use even-split
    if seq_len < 50 || period >= seq_len / 2 {
        return even_split(hor_seq, period);
    }

    // Find anchor using the known period (no autocorr re-detection)
    match find_anchor_by_period_with_fallback(seq_bytes, period) {
        Some(anchor_result) if anchor_result.positions.len() >= 2 => {
            let boundaries = split_multiplets(seq_bytes, &anchor_result.positions, period);

            if boundaries.len() > 1 {
                boundaries.iter()
                    .map(|b| {
                        let seq = hor_seq[b.start..b.end].to_string();
                        let source = if b.source.contains("flank") {
                            "recursive_flank".to_string()
                        } else {
                            "recursive_anchor".to_string()
                        };
                        (seq, source)
                    })
                    .collect()
            } else {
                even_split(hor_seq, period)
            }
        }
        _ => even_split(hor_seq, period),
    }
}

/// Split sequence evenly by period
fn even_split(seq: &str, period: usize) -> Vec<(String, String)> {
    if period == 0 || period >= seq.len() {
        return vec![(seq.to_string(), "recursive_split".to_string())];
    }

    let mut result = Vec::new();
    let mut pos = 0;

    while pos < seq.len() {
        let end = (pos + period).min(seq.len());
        result.push((seq[pos..end].to_string(), "recursive_split".to_string()));
        pos = end;
    }

    result
}

/// Compute the maximum autocorrelation in a reasonable period range
fn compute_max_autocorr(seq: &[u8], min_period: usize) -> f64 {
    let max_period = seq.len() / 2;
    if max_period < min_period {
        return 0.0;
    }

    let random_exp = random_expectation(seq);
    let mut max_ac = 0.0f64;

    for d in min_period..=max_period.min(500) {
        let ac = autocorrelation(seq, d);
        let excess = ac - random_exp;
        if excess > max_ac {
            max_ac = excess;
        }
    }

    // Return actual autocorrelation, not excess
    max_ac + random_exp
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_base_monomer_no_periodicity() {
        // Short sequence without enough length for decomposition
        // The key is that it should just verify total length preservation
        let hor = "ACGTTAGCAGTCGATCAGTCAGTCGATCGATCGATCAGTCAGTCAGTCAGT";
        let result = decompose_hors_to_base(&[hor.to_string()], 5, 0.5);

        // Total length should be preserved regardless of decomposition
        let total_len: usize = result.base_monomers.iter().map(|m| m.sequence.len()).sum();
        assert_eq!(total_len, hor.len(), "Total length should be preserved");

        // All monomers should have hor_idx = 0
        for mono in &result.base_monomers {
            assert_eq!(mono.hor_idx, 0);
        }
    }

    #[test]
    fn test_hor_with_periodicity() {
        // HOR with 10x copies of 8bp monomer
        let monomer = "ACGTTAGC";
        let hor: String = (0..10).map(|_| monomer).collect();

        let result = decompose_hors_to_base(&[hor.clone()], 5, 0.5);

        // Should detect periodicity and decompose
        assert!(result.base_monomers.len() >= 3, "Expected at least 3 monomers, got {}", result.base_monomers.len());

        // Verify total length is preserved
        let total_len: usize = result.base_monomers.iter().map(|m| m.sequence.len()).sum();
        assert_eq!(total_len, hor.len());
    }

    #[test]
    fn test_short_sequence() {
        // Sequence too short for meaningful decomposition
        let hor = "ACGT";
        let result = decompose_hors_to_base(&[hor.to_string()], 5, 0.5);

        // Should return as-is
        assert_eq!(result.base_monomers.len(), 1);
        assert_eq!(result.base_monomers[0].sequence, hor);
    }

    #[test]
    fn test_multiple_hors() {
        // Multiple HOR monomers
        let hor1: String = (0..5).map(|_| "ACGT").collect();
        let hor2: String = (0..5).map(|_| "TTAGGG").collect();

        let result = decompose_hors_to_base(&[hor1.clone(), hor2.clone()], 5, 0.5);

        // Should decompose both (total length preserved)
        let total_len: usize = result.base_monomers.iter().map(|m| m.sequence.len()).sum();
        assert_eq!(total_len, hor1.len() + hor2.len(), "Total length should be preserved");

        // Check hor_idx is preserved
        let hor0_monomers: Vec<_> = result.base_monomers.iter().filter(|m| m.hor_idx == 0).collect();
        let hor1_monomers: Vec<_> = result.base_monomers.iter().filter(|m| m.hor_idx == 1).collect();

        assert!(!hor0_monomers.is_empty());
        assert!(!hor1_monomers.is_empty());
    }

    #[test]
    fn test_alpha_satellite_like() {
        // Simulate a 512bp HOR (3x ~171bp monomers)
        // Use a realistic alpha satellite monomer sequence
        let base_monomer = "AATGGTTTCAAAGTTATTTTTAAAATTGTAAAAAGACTTTCGATTTTTTTTATCTTTTTGACTGAAAATATTTCTTTTGTAAGATTTGAGATCTCAGTGTATAATCCTTTCATAAAAAATTAAAATTGGGATATTGAGGGAATAACATTCTTATG";
        let hor: String = (0..3).map(|_| base_monomer).collect();

        let result = decompose_hors_to_base(&[hor.clone()], 5, 0.5);

        // Total length should be preserved
        let total_len: usize = result.base_monomers.iter().map(|m| m.sequence.len()).sum();
        assert_eq!(total_len, hor.len(), "Total length should be preserved");

        // Should produce at least one monomer
        assert!(!result.base_monomers.is_empty(), "Should produce at least one monomer");
    }

    #[test]
    fn test_threshold_behavior() {
        // Sequence with moderate periodicity
        let monomer = "ACGTTAGC";
        let hor: String = (0..20).map(|_| monomer).collect();

        // Low threshold - should decompose
        let result_low = decompose_hors_to_base(&[hor.clone()], 5, 0.3);
        // High threshold - might not decompose (depends on autocorr)
        let result_high = decompose_hors_to_base(&[hor.clone()], 5, 0.99);

        // With perfect repeats, autocorr should be high, so both should decompose
        assert!(result_low.base_monomers.len() >= result_high.base_monomers.len());

        // Both should preserve total length
        let total_low: usize = result_low.base_monomers.iter().map(|m| m.sequence.len()).sum();
        let total_high: usize = result_high.base_monomers.iter().map(|m| m.sequence.len()).sum();
        assert_eq!(total_low, hor.len());
        assert_eq!(total_high, hor.len());
    }

    #[test]
    fn test_length_preservation() {
        // Key test: regardless of decomposition, total length must be preserved
        let sequences = vec![
            "ACGT".to_string(),
            "ACGTACGTACGT".to_string(),
            "TTAGGGTTAGGGTTAGGG".to_string(),
            "AATGGAATGGAATGGAATGGAATGG".to_string(),
        ];

        for seq in &sequences {
            let result = decompose_hors_to_base(&[seq.clone()], 5, 0.5);
            let total_len: usize = result.base_monomers.iter().map(|m| m.sequence.len()).sum();
            assert_eq!(total_len, seq.len(), "Length not preserved for: {}", seq);
        }
    }
}
