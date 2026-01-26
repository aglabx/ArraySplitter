//! Anchor discovery with known period constraint
//!
//! Given a known period P from autocorrelation, finds the best k-mer anchor
//! that occurs approximately once per P bases with regular spacing.

use std::collections::HashMap;

/// Result of anchor discovery
#[derive(Debug, Clone)]
pub struct AnchorResult {
    /// The anchor k-mer sequence
    pub kmer: String,
    /// Positions of anchor occurrences in the sequence
    pub positions: Vec<usize>,
    /// Length of the anchor k-mer
    pub k: usize,
    /// Score = uniqueness × regularity
    pub score: f64,
    /// Fraction of occurrences that are the sole one in their P-window
    pub uniqueness: f64,
    /// Fraction of gaps that are approximately P (±30%)
    pub regularity: f64,
}

/// Find the best anchor k-mer with a known period constraint.
///
/// Strategy:
/// - N_expected = seq.len() / period
/// - For each k in [min_k, max_k]:
///   - Count all k-mers
///   - Filter: freq in [0.5*N, 1.5*N]
///   - For each candidate: compute uniqueness and regularity scores
///   - Score = uniqueness × regularity
/// - Return best candidate (highest score × longest k, score > 0.2)
pub fn find_anchor_by_period(
    seq: &[u8],
    period: usize,
    min_k: usize,
    max_k: usize,
) -> Option<AnchorResult> {
    if seq.len() < period * 2 {
        return None; // Need at least 2 copies
    }

    let n_expected = seq.len() / period;
    if n_expected < 2 {
        return None;
    }

    let freq_min = (n_expected as f64 * 0.5) as usize;
    let freq_max = (n_expected as f64 * 1.5) as usize;

    // For short periods (e.g., TTAGGG=6bp), allow k up to period-1
    // For long periods, cap at max_k (typically 15)
    let effective_max_k = max_k.min(period.saturating_sub(1));
    let effective_min_k = min_k.min(effective_max_k);

    let mut best_result: Option<AnchorResult> = None;
    let mut best_combined_score: f64 = 0.0;

    for k in effective_min_k..=effective_max_k {
        if k > seq.len() {
            break;
        }

        // Count k-mers and record positions
        let mut kmer_positions: HashMap<&[u8], Vec<usize>> = HashMap::new();
        for i in 0..=(seq.len() - k) {
            let kmer = &seq[i..i + k];
            // Skip k-mers with non-ACGT characters
            if kmer.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')) {
                kmer_positions.entry(kmer).or_default().push(i);
            }
        }

        // Filter by frequency and sort for deterministic iteration
        let mut candidates: Vec<(&[u8], &Vec<usize>)> = kmer_positions
            .iter()
            .filter(|(_, positions)| {
                let freq = positions.len();
                freq >= freq_min && freq <= freq_max
            })
            .map(|(kmer, positions)| (*kmer, positions))
            .collect();

        // Sort by k-mer bytes for deterministic tie-breaking
        candidates.sort_by_key(|(kmer, _)| *kmer);

        for (kmer_bytes, positions) in candidates {
            let uniqueness = compute_uniqueness(positions, period);
            let regularity = compute_regularity(positions, period);
            let score = uniqueness * regularity;

            if score < 0.2 {
                continue;
            }

            // Combined score: favor longer k-mers with good scores
            // Use score * sqrt(k) to prefer longer anchors without overwhelming score
            let combined = score * (k as f64).sqrt();

            if combined > best_combined_score {
                best_combined_score = combined;
                let kmer_str = String::from_utf8_lossy(kmer_bytes).to_uppercase();
                best_result = Some(AnchorResult {
                    kmer: kmer_str,
                    positions: positions.clone(),
                    k,
                    score,
                    uniqueness,
                    regularity,
                });
            }
        }
    }

    best_result
}

/// Compute uniqueness: fraction of occurrences that are the sole one in their P-window.
/// A "P-window" is [pos - P/2, pos + P/2].
/// High uniqueness means the k-mer doesn't appear multiple times per monomer.
fn compute_uniqueness(positions: &[usize], period: usize) -> f64 {
    if positions.is_empty() {
        return 0.0;
    }

    let half_p = period / 2;
    let mut unique_count = 0;

    for (i, &pos) in positions.iter().enumerate() {
        let window_start = pos.saturating_sub(half_p);
        let window_end = pos + half_p;

        let mut is_unique = true;
        // Check neighbors
        if i > 0 && positions[i - 1] >= window_start {
            is_unique = false;
        }
        if i + 1 < positions.len() && positions[i + 1] <= window_end {
            is_unique = false;
        }

        if is_unique {
            unique_count += 1;
        }
    }

    unique_count as f64 / positions.len() as f64
}

/// Compute regularity: fraction of consecutive gaps that are approximately P (±30%).
fn compute_regularity(positions: &[usize], period: usize) -> f64 {
    if positions.len() < 2 {
        return 0.0;
    }

    let tolerance = (period as f64 * 0.3) as usize;
    let min_gap = period.saturating_sub(tolerance);
    let max_gap = period + tolerance;

    let mut regular_count = 0;
    let total_gaps = positions.len() - 1;

    for i in 0..total_gaps {
        let gap = positions[i + 1] - positions[i];
        // Gap should be approximately P, or a multiple of P (for missed anchors)
        if gap >= min_gap && gap <= max_gap {
            regular_count += 1;
        } else {
            // Check if it's a multiple of P (missed anchor in between)
            let n_periods = (gap as f64 / period as f64).round() as usize;
            if n_periods >= 2 {
                let expected = n_periods * period;
                let deviation = if gap > expected { gap - expected } else { expected - gap };
                if deviation <= tolerance {
                    regular_count += 1; // Count as regular (multiplet gap)
                }
            }
        }
    }

    regular_count as f64 / total_gaps as f64
}

/// Find anchor with fallback: if primary anchor search fails,
/// try with relaxed parameters.
pub fn find_anchor_by_period_with_fallback(
    seq: &[u8],
    period: usize,
) -> Option<AnchorResult> {
    // First try with default parameters
    if let Some(result) = find_anchor_by_period(seq, period, 5, 15) {
        return Some(result);
    }

    // Fallback: try shorter k-mers
    if let Some(result) = find_anchor_by_period(seq, period, 3, 8) {
        return Some(result);
    }

    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_anchor_perfect_repeat() {
        // "ACGTTAGC" repeated 20 times → period 8, anchor should be found
        let unit = b"ACGTTAGC";
        let mut seq = Vec::new();
        for _ in 0..20 {
            seq.extend_from_slice(unit);
        }

        let result = find_anchor_by_period(&seq, 8, 3, 10);
        assert!(result.is_some());
        let anchor = result.unwrap();
        assert!(anchor.score > 0.5);
        assert!(anchor.positions.len() >= 19); // ~One per copy (boundary k-mers may have N-1)
        // Verify positions are approximately P apart
        for i in 1..anchor.positions.len() {
            let gap = anchor.positions[i] - anchor.positions[i - 1];
            assert_eq!(gap, 8);
        }
    }

    #[test]
    fn test_find_anchor_with_mutations() {
        // Repeat with ~10% mutations - anchor should still be found
        let unit = b"ACGTTAGCAA";
        let mut seq = Vec::new();
        for _ in 0..20 {
            seq.extend_from_slice(unit);
        }
        // Mutate some positions (but not the anchor region)
        for i in (5..seq.len()).step_by(11) {
            if i < seq.len() {
                seq[i] = b'N'; // Will be skipped as non-ACGT
            }
        }
        // Remove N's and replace with valid bases
        for b in seq.iter_mut() {
            if *b == b'N' {
                *b = b'G';
            }
        }

        let result = find_anchor_by_period(&seq, 10, 3, 8);
        assert!(result.is_some());
        let anchor = result.unwrap();
        assert!(anchor.score > 0.2);
    }

    #[test]
    fn test_compute_uniqueness_perfect() {
        // Positions exactly P apart → all unique in their windows
        let positions: Vec<usize> = (0..10).map(|i| i * 100).collect();
        let u = compute_uniqueness(&positions, 100);
        assert!((u - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_compute_uniqueness_duplicates() {
        // Some positions are too close (within P/2)
        let positions = vec![0, 10, 100, 110, 200, 300];
        let u = compute_uniqueness(&positions, 100);
        assert!(u < 1.0); // Some are not unique
    }

    #[test]
    fn test_compute_regularity_perfect() {
        let positions: Vec<usize> = (0..10).map(|i| i * 100).collect();
        let r = compute_regularity(&positions, 100);
        assert!((r - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_compute_regularity_with_multiplets() {
        // Gaps: 100, 200 (missed one), 100, 100, 300 (missed two)
        let positions = vec![0, 100, 300, 400, 500, 800];
        let r = compute_regularity(&positions, 100);
        // Gaps: 100(ok), 200(2P ok), 100(ok), 100(ok), 300(3P ok)
        assert!((r - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_find_anchor_too_short_sequence() {
        let seq = b"ACGT";
        let result = find_anchor_by_period(seq, 100, 5, 15);
        assert!(result.is_none());
    }

    #[test]
    fn test_find_anchor_period_170() {
        // Simulate alpha satellite-like array: period ~170
        let mut unit = Vec::new();
        // Build a 170bp "monomer" with some structure
        for i in 0..170 {
            unit.push(match i % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            });
        }
        // Insert a distinctive anchor at position 0
        unit[0] = b'T';
        unit[1] = b'T';
        unit[2] = b'A';
        unit[3] = b'G';
        unit[4] = b'G';
        unit[5] = b'G';

        let mut seq = Vec::new();
        for _ in 0..50 {
            seq.extend_from_slice(&unit);
        }

        let result = find_anchor_by_period(&seq, 170, 5, 15);
        assert!(result.is_some());
        let anchor = result.unwrap();
        assert!(anchor.positions.len() >= 40); // Should find most copies
    }
}
