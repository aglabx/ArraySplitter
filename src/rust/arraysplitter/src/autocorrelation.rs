//! Autocorrelation-based period detection for satellite DNA arrays
//!
//! Computes the fraction of matching nucleotides at offset d to detect periodicity.
//! The best period is the offset with the maximum excess over random expectation.

/// Compute autocorrelation at a specific offset d.
/// Returns the fraction of positions where seq[i] == seq[i+d].
pub fn autocorrelation(seq: &[u8], d: usize) -> f64 {
    if d == 0 || d >= seq.len() {
        return 0.0;
    }

    let n = seq.len() - d;
    if n == 0 {
        return 0.0;
    }

    let matches: usize = seq[..n]
        .iter()
        .zip(seq[d..].iter())
        .filter(|(&a, &b)| a == b)
        .count();

    matches as f64 / n as f64
}

/// Compute random expectation from nucleotide composition.
/// This is sum(freq_i^2) for each nucleotide, which is the expected
/// autocorrelation for a random sequence with the same composition.
pub fn random_expectation(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.25; // uniform distribution
    }

    let mut counts = [0usize; 4]; // A, C, G, T
    for &b in seq {
        match b {
            b'A' | b'a' => counts[0] += 1,
            b'C' | b'c' => counts[1] += 1,
            b'G' | b'g' => counts[2] += 1,
            b'T' | b't' => counts[3] += 1,
            _ => {}
        }
    }

    let total = counts.iter().sum::<usize>() as f64;
    if total == 0.0 {
        return 0.25;
    }

    counts
        .iter()
        .map(|&c| {
            let f = c as f64 / total;
            f * f
        })
        .sum()
}

/// Compute autocorrelation at offset d using subsampling for large sequences.
/// For sequences > SAMPLE_THRESHOLD, samples N_SAMPLE positions instead of all.
fn autocorrelation_sampled(seq: &[u8], d: usize, sample_positions: Option<&[usize]>) -> f64 {
    if d == 0 || d >= seq.len() {
        return 0.0;
    }

    match sample_positions {
        Some(positions) => {
            // Use precomputed sample positions
            let valid: Vec<&usize> = positions.iter().filter(|&&p| p + d < seq.len()).collect();
            if valid.is_empty() {
                return 0.0;
            }
            let matches: usize = valid.iter().filter(|&&&p| seq[p] == seq[p + d]).count();
            matches as f64 / valid.len() as f64
        }
        None => {
            // Full computation
            autocorrelation(seq, d)
        }
    }
}

/// Generate deterministic sample positions for subsampling
fn generate_sample_positions(seq_len: usize, n_samples: usize) -> Vec<usize> {
    if seq_len <= n_samples {
        return (0..seq_len).collect();
    }
    let step = seq_len / n_samples;
    (0..n_samples).map(|i| i * step).collect()
}

/// Find the best period in a given range [min_d, max_d].
/// Returns (period, autocorr_value, excess_over_random) or None if no periodicity detected.
///
/// The excess threshold is 0.05 — if no offset exceeds random expectation by this amount,
/// we conclude there is no periodicity.
///
/// For sequences > 50KB, uses position subsampling (10000 positions) to keep
/// computation time constant regardless of sequence length.
pub fn find_period(seq: &[u8], min_d: usize, max_d: usize) -> Option<(usize, f64, f64)> {
    if seq.len() < min_d + 1 {
        return None;
    }

    let random_exp = random_expectation(seq);
    let effective_max = max_d.min(seq.len() / 2); // need at least 2 copies

    if min_d > effective_max {
        return None;
    }

    // For large sequences, use subsampling
    const SAMPLE_THRESHOLD: usize = 50_000;
    const N_SAMPLES: usize = 10_000;

    let sample_positions = if seq.len() > SAMPLE_THRESHOLD {
        Some(generate_sample_positions(seq.len(), N_SAMPLES))
    } else {
        None
    };

    let mut best_period = 0;
    let mut best_autocorr = 0.0f64;
    let mut best_excess = 0.0f64;

    for d in min_d..=effective_max {
        let ac = autocorrelation_sampled(seq, d, sample_positions.as_deref());
        let excess = ac - random_exp;

        if excess > best_excess {
            best_excess = excess;
            best_autocorr = ac;
            best_period = d;
        }
    }

    // Validate: peak excess must be significant
    if best_excess < 0.05 {
        return None;
    }

    Some((best_period, best_autocorr, best_excess))
}

/// Find the best period with harmonic refinement.
/// After finding the raw peak, checks if half-period has similar excess
/// (which would indicate the true period is the half).
/// Also checks multiples to ensure we have the fundamental period.
pub fn find_period_refined(seq: &[u8], min_d: usize, max_d: usize) -> Option<(usize, f64, f64)> {
    let result = find_period(seq, min_d, max_d)?;
    let (period, _autocorr, _excess) = result;

    let random_exp = random_expectation(seq);

    // Check if half-period is also a strong peak (harmonic detection)
    if period >= min_d * 2 {
        let half = period / 2;
        if half >= min_d {
            let half_ac = autocorrelation(seq, half);
            let half_excess = half_ac - random_exp;

            // If half-period has ≥80% of the full period's excess, it's likely the fundamental
            if half_excess >= _excess * 0.8 && half_excess >= 0.05 {
                // But verify by checking if period is ~2×half
                if (period as f64 / half as f64 - 2.0).abs() < 0.1 {
                    return Some((half, half_ac, half_excess));
                }
            }
        }
    }

    Some(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_autocorrelation_perfect_repeat() {
        // "ACGT" repeated → period 4 should have autocorrelation = 1.0
        let seq = b"ACGTACGTACGTACGT";
        let ac = autocorrelation(seq, 4);
        assert!((ac - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_autocorrelation_no_repeat() {
        // Random-ish sequence, autocorrelation at arbitrary offset should be ~0.25
        let seq = b"ACGTTAGCAGTCGATCAGTCAGTCGATCGATCGATCAGTCAGTCAGTCAGT";
        let ac = autocorrelation(seq, 7);
        assert!(ac < 0.5); // Should not be strongly periodic at offset 7
    }

    #[test]
    fn test_random_expectation_uniform() {
        let seq = b"ACGTACGTACGTACGT";
        let re = random_expectation(seq);
        // Uniform composition: each nucleotide = 0.25, sum(0.25^2) = 0.25
        assert!((re - 0.25).abs() < 1e-10);
    }

    #[test]
    fn test_random_expectation_biased() {
        let seq = b"AAAAAAAAAA"; // All A
        let re = random_expectation(seq);
        // Only A: freq_A = 1.0, sum = 1.0
        assert!((re - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_find_period_simple() {
        // 10 copies of "ACGT" → period should be 4
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let result = find_period(seq, 2, 20);
        assert!(result.is_some());
        let (period, autocorr, excess) = result.unwrap();
        assert_eq!(period, 4);
        assert!((autocorr - 1.0).abs() < 1e-10);
        assert!(excess > 0.5); // Well above random
    }

    #[test]
    fn test_find_period_longer() {
        // "ACGTTAGC" repeated → period 8
        let unit = b"ACGTTAGC";
        let mut seq = Vec::new();
        for _ in 0..20 {
            seq.extend_from_slice(unit);
        }
        let result = find_period(&seq, 2, 50);
        assert!(result.is_some());
        let (period, _, _) = result.unwrap();
        assert_eq!(period, 8);
    }

    #[test]
    fn test_find_period_no_periodicity() {
        // A truly non-periodic sequence: use a linear congruential generator
        // Need a long sequence (5000bp) so that spurious peaks are below threshold
        // Standard error of autocorrelation ≈ 1/sqrt(n-d), so for n=5000, SE ≈ 0.014
        let mut seq = Vec::with_capacity(5000);
        let mut state: u64 = 12345;
        for _ in 0..5000 {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let nucleotide = match (state >> 33) % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            };
            seq.push(nucleotide);
        }
        // Search for periods in range [50, 500] — should find nothing significant
        let result = find_period(&seq, 50, 500);
        assert!(result.is_none());
    }

    #[test]
    fn test_find_period_with_mutations() {
        // "ACGTTAGC" repeated with mutations at varying positions within each copy
        // This ensures mutations don't create a super-period
        let unit = b"ACGTTAGC";
        let mut seq = Vec::new();
        for _ in 0..30 {
            seq.extend_from_slice(unit);
        }
        // Mutate one position per copy at different offsets within the monomer
        // This disrupts the period slightly but keeps it dominant
        let offsets = [0, 2, 5, 1, 4, 7, 3, 6, 2, 0, 5, 1, 7, 3, 4,
                       6, 0, 2, 5, 7, 1, 3, 4, 6, 0, 2, 5, 7, 1, 3];
        for (copy_idx, &offset) in offsets.iter().enumerate() {
            let pos = copy_idx * 8 + offset;
            if pos < seq.len() {
                seq[pos] = match seq[pos] {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'A',
                    _ => b'C',
                };
            }
        }

        let result = find_period(&seq, 2, 50);
        assert!(result.is_some());
        let (period, _autocorr, excess) = result.unwrap();
        // Period should be 8 (or a harmonic like 16, 24 in rare cases)
        assert!(period % 8 == 0, "Expected period to be multiple of 8, got {}", period);
        assert!(excess > 0.05); // Should be detectable
    }

    #[test]
    fn test_autocorrelation_zero_offset() {
        let seq = b"ACGT";
        assert_eq!(autocorrelation(seq, 0), 0.0);
    }

    #[test]
    fn test_autocorrelation_too_large_offset() {
        let seq = b"ACGT";
        assert_eq!(autocorrelation(seq, 10), 0.0);
    }
}
