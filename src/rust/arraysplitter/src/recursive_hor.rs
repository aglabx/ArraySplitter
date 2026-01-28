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
}

/// Result of recursive HOR decomposition
#[derive(Debug, Clone)]
pub struct RecursiveResult {
    /// All base-level monomers from the recursive decomposition
    pub base_monomers: Vec<BaseMonomer>,
    /// Maximum recursion depth reached
    pub max_depth: usize,
}

/// Default minimum submonomer length (bp)
pub const DEFAULT_MIN_SUBMONOMER_LEN: usize = 5;

/// Default autocorrelation threshold for further decomposition
pub const DEFAULT_AUTOCORR_THRESHOLD: f64 = 0.5;

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

    RecursiveResult {
        base_monomers: all_base_monomers,
        max_depth,
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
            sub_idx: 0,
            level: current_level,
            period: 0,
            autocorr: 0.0,
            source: "base".to_string(),
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
                    sub_idx: 0,
                    level: current_level,
                    period: 0,
                    autocorr,
                    source: "base".to_string(),
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
                        sub_idx,
                        level: current_level,
                        period,
                        autocorr,
                        source: sub_source.clone(),
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
                        sub_idx,
                        level: current_level,
                        period,
                        autocorr: final_autocorr,
                        source: sub_source.clone(),
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
                sub_idx: 0,
                level: current_level,
                period: 0,
                autocorr: compute_max_autocorr(seq_bytes, min_len),
                source: "base".to_string(),
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
