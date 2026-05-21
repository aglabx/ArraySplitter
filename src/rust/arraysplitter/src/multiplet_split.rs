//! Multiplet splitting using template edit-distance alignment
//!
//! After anchor-based decomposition, some segments span multiple periods.
//! This module splits those multiplet segments by sliding a template monomer
//! and finding optimal split points.

/// A monomer boundary with its source label
#[derive(Debug, Clone)]
pub struct MonomerBoundary {
    /// Start position in the original sequence
    pub start: usize,
    /// End position (exclusive) in the original sequence
    pub end: usize,
    /// Source label: "anchor", "split_2x", "split_3x", etc.
    pub source: String,
}

/// Split multiplet segments using template ED alignment.
///
/// Pipeline:
/// 1. Build segments from anchor positions (add 0 and seq_len as boundaries)
/// 2. Classify: 1×P = monomer (source="anchor"), N×P = multiplet
/// 3. From 1×P monomers, select template (closest length to period)
/// 4. For each multiplet, find optimal split points using ED minimization
pub fn split_multiplets(
    seq: &[u8],
    anchor_positions: &[usize],
    period: usize,
) -> Vec<MonomerBoundary> {
    if anchor_positions.is_empty() {
        // No anchors found - return entire sequence as single monomer
        return vec![MonomerBoundary {
            start: 0,
            end: seq.len(),
            source: "no_anchor".to_string(),
        }];
    }

    // Build boundaries: [0, anchor_pos_0, anchor_pos_1, ..., seq_len]
    let mut boundaries: Vec<usize> = Vec::new();
    boundaries.push(0);
    for &pos in anchor_positions {
        if pos > 0 && pos < seq.len() && (boundaries.is_empty() || *boundaries.last().unwrap() != pos) {
            boundaries.push(pos);
        }
    }
    if *boundaries.last().unwrap() != seq.len() {
        boundaries.push(seq.len());
    }

    // Build segments with classification
    let tolerance = (period as f64 * 0.3) as usize;
    let mut segments: Vec<(usize, usize, usize)> = Vec::new(); // (start, end, n_periods)

    for i in 0..boundaries.len() - 1 {
        let start = boundaries[i];
        let end = boundaries[i + 1];
        let len = end - start;

        let n_periods = if period > 0 {
            ((len as f64 / period as f64) + 0.5) as usize
        } else {
            1
        };

        // Validate: is this approximately n_periods × P?
        let n_periods = if n_periods == 0 {
            1
        } else {
            let expected_len = n_periods * period;
            let deviation = if len > expected_len {
                len - expected_len
            } else {
                expected_len - len
            };
            if deviation <= tolerance * n_periods {
                n_periods
            } else {
                // Recompute more carefully
                let best_n = (len as f64 / period as f64).round() as usize;
                if best_n == 0 { 1 } else { best_n }
            }
        };

        segments.push((start, end, n_periods));
    }

    // Collect 1×P segments as templates
    let single_segments: Vec<(usize, usize)> = segments
        .iter()
        .filter(|(_, _, n)| *n == 1)
        .map(|(s, e, _)| (*s, *e))
        .collect();

    // Select template: the 1×P monomer closest in length to P
    let template: Option<&[u8]> = if !single_segments.is_empty() {
        let best_idx = single_segments
            .iter()
            .enumerate()
            .min_by_key(|(_, (s, e))| {
                let len = e - s;
                if len > period { len - period } else { period - len }
            })
            .map(|(idx, _)| idx)
            .unwrap();

        let (ts, te) = single_segments[best_idx];
        Some(&seq[ts..te])
    } else {
        None
    };

    // Build final monomer boundaries
    let mut result: Vec<MonomerBoundary> = Vec::new();

    for (start, end, n_periods) in &segments {
        if *n_periods <= 1 {
            // Single monomer from anchor
            let source = if *start == 0 && !anchor_positions.contains(start) {
                "left_flank".to_string()
            } else if *end == seq.len() && !anchor_positions.contains(start) {
                "right_flank".to_string()
            } else {
                "anchor".to_string()
            };
            result.push(MonomerBoundary {
                start: *start,
                end: *end,
                source,
            });
        } else {
            // Multiplet: needs splitting
            let sub_boundaries = split_single_multiplet(
                seq,
                *start,
                *end,
                *n_periods,
                period,
                template,
            );
            result.extend(sub_boundaries);
        }
    }

    result
}

/// Split a single multiplet segment into n_expected monomers.
///
/// Strategy:
/// - If template is available, use ED-guided splitting
/// - Otherwise, split evenly based on period
fn split_single_multiplet(
    seq: &[u8],
    start: usize,
    end: usize,
    n_expected: usize,
    period: usize,
    template: Option<&[u8]>,
) -> Vec<MonomerBoundary> {
    let segment_len = end - start;

    if n_expected <= 1 || segment_len < period {
        return vec![MonomerBoundary {
            start,
            end,
            source: format!("split_{}x", n_expected),
        }];
    }

    let split_positions = if let Some(tmpl) = template {
        find_split_positions_by_ed(seq, start, end, n_expected, period, tmpl)
    } else {
        find_split_positions_even(start, end, n_expected)
    };

    // Build monomer boundaries from split positions
    let mut boundaries: Vec<usize> = Vec::new();
    boundaries.push(start);
    boundaries.extend(split_positions.iter());
    boundaries.push(end);

    let source_label = format!("split_{}x", n_expected);

    boundaries
        .windows(2)
        .map(|w| MonomerBoundary {
            start: w[0],
            end: w[1],
            source: source_label.clone(),
        })
        .collect()
}

/// Find split positions using edit-distance guided approach.
///
/// For each expected boundary position (segment_len/n * k for k=1..n-1):
/// - Search in a ±P/4 window around the expected position
/// - Slide the template and compute ED at each candidate position
/// - Choose position with minimum ED between the template and the subsequence
///
/// Optimization: for periods > 200bp, only compare the first COMPARE_PREFIX bytes
/// of template and candidate to keep ED computation O(prefix²) instead of O(period²).
fn find_split_positions_by_ed(
    seq: &[u8],
    start: usize,
    end: usize,
    n_expected: usize,
    period: usize,
    template: &[u8],
) -> Vec<usize> {
    let segment_len = end - start;
    let search_window = period / 4;

    // For large periods, truncate comparison to prefix for speed
    const COMPARE_PREFIX: usize = 200;
    let effective_tmpl = if template.len() > COMPARE_PREFIX {
        &template[..COMPARE_PREFIX]
    } else {
        template
    };
    let compare_len = effective_tmpl.len();

    let mut split_positions: Vec<usize> = Vec::new();

    for k in 1..n_expected {
        // Expected position of k-th boundary (relative to start)
        let expected_rel = (segment_len as f64 * k as f64 / n_expected as f64) as usize;
        let expected_abs = start + expected_rel;

        // Search window
        let window_start = expected_abs.saturating_sub(search_window).max(start + 1);
        let window_end = (expected_abs + search_window).min(end.saturating_sub(1));

        if window_start >= window_end {
            // Fallback to expected position
            split_positions.push(expected_abs.min(end - 1));
            continue;
        }

        // Find position with minimum ED between template prefix and the candidate prefix
        let mut best_pos = expected_abs;
        let mut best_ed = usize::MAX;

        for pos in window_start..=window_end {
            let monomer_end = (pos + compare_len).min(end);
            let candidate = &seq[pos..monomer_end];

            if candidate.is_empty() {
                continue;
            }

            let ed = fast_edit_distance_bounded(effective_tmpl, candidate, best_ed);
            if ed < best_ed {
                best_ed = ed;
                best_pos = pos;
            }
        }

        // Don't place split too close to previous split
        if let Some(&prev) = split_positions.last() {
            if best_pos <= prev + period / 4 {
                best_pos = prev + period / 2; // Minimum spacing
                if best_pos >= end {
                    continue;
                }
            }
        }

        split_positions.push(best_pos);
    }

    split_positions
}

/// Find split positions by even division (fallback when no template available).
fn find_split_positions_even(start: usize, end: usize, n_expected: usize) -> Vec<usize> {
    let segment_len = end - start;
    (1..n_expected)
        .map(|k| start + (segment_len * k / n_expected))
        .collect()
}

/// Bounded edit distance computation using triple_accel SIMD (k-bounded variant).
/// Returns `bound + 1` when ED exceeds the bound — preserves the existing sentinel semantics.
///
/// Pre/post-conditions match the previous naive bounded DP:
///   - len_diff > bound  → return bound + 1
///   - len > 5000        → fall back to approximate_edit_distance (unchanged)
///   - empty inputs      → return the non-empty length
///   - otherwise         → exact ED if ≤ bound, else bound + 1
fn fast_edit_distance_bounded(s1: &[u8], s2: &[u8], bound: usize) -> usize {
    let len1 = s1.len();
    let len2 = s2.len();

    let len_diff = if len1 > len2 { len1 - len2 } else { len2 - len1 };
    if len_diff > bound {
        return bound + 1;
    }

    if len1 > 5000 || len2 > 5000 {
        return approximate_edit_distance(s1, s2);
    }

    if len1 == 0 { return len2; }
    if len2 == 0 { return len1; }

    match triple_accel::levenshtein::levenshtein_simd_k(s1, s2, bound as u32) {
        Some(ed) => ed as usize,
        None => bound + 1,
    }
}

/// Approximate edit distance for very long sequences.
/// Uses block-based comparison.
fn approximate_edit_distance(s1: &[u8], s2: &[u8]) -> usize {
    let min_len = s1.len().min(s2.len());
    let max_len = s1.len().max(s2.len());
    let len_diff = max_len - min_len;

    // Count mismatches in the overlapping region
    let mismatches: usize = s1[..min_len]
        .iter()
        .zip(s2[..min_len].iter())
        .filter(|(&a, &b)| a != b)
        .count();

    mismatches + len_diff
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_split_multiplets_no_anchors() {
        let seq = b"ACGTACGTACGT";
        let result = split_multiplets(seq, &[], 4);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].start, 0);
        assert_eq!(result[0].end, 12);
        assert_eq!(result[0].source, "no_anchor");
    }

    #[test]
    fn test_split_multiplets_all_single() {
        // Anchors every 4bp in a 20bp sequence
        let seq = b"ACGTACGTACGTACGTACGT";
        let anchors = vec![0, 4, 8, 12, 16];
        let result = split_multiplets(seq, &anchors, 4);
        assert_eq!(result.len(), 5);
        for (i, m) in result.iter().enumerate() {
            assert_eq!(m.start, i * 4);
            assert_eq!(m.end, (i + 1) * 4);
            assert_eq!(m.source, "anchor");
        }
    }

    #[test]
    fn test_split_multiplets_with_doublet() {
        // Sequence with period 4, but one anchor missing → doublet
        let seq = b"ACGTACGTACGTACGTACGT"; // 20bp, period 4
        // Anchors at 0, 4, 12, 16 (missing 8 → segment 4..12 is 2×P)
        let anchors = vec![0, 4, 12, 16];
        let result = split_multiplets(seq, &anchors, 4);

        // Should have: [0,4), [4,8), [8,12), [12,16), [16,20)
        // The [4,12) segment should be split into two
        assert_eq!(result.len(), 5);

        // First segment: anchor
        assert_eq!(result[0].start, 0);
        assert_eq!(result[0].end, 4);
        assert_eq!(result[0].source, "anchor");

        // Middle segments: split from doublet
        assert_eq!(result[1].source, "split_2x");
        assert_eq!(result[2].source, "split_2x");
    }

    #[test]
    fn test_split_multiplets_with_flank() {
        // First segment is a flank (anchor doesn't start at 0)
        let seq = b"XXACGTACGTACGT"; // 14bp
        let anchors = vec![2, 6, 10];
        let result = split_multiplets(seq, &anchors, 4);

        // First segment [0,2) should be left_flank
        assert!(result.len() >= 4);
        assert_eq!(result[0].source, "left_flank");
    }

    #[test]
    fn test_fast_edit_distance_bounded() {
        let s1 = b"ACGT";
        let s2 = b"ACGT";
        assert_eq!(fast_edit_distance_bounded(s1, s2, 10), 0);

        let s1 = b"ACGT";
        let s2 = b"ACGA";
        assert_eq!(fast_edit_distance_bounded(s1, s2, 10), 1);

        let s1 = b"ACGT";
        let s2 = b"TGCA";
        let ed = fast_edit_distance_bounded(s1, s2, 10);
        assert!(ed <= 4);
    }

    #[test]
    fn test_fast_edit_distance_bounded_early_exit() {
        let s1 = b"AAAA";
        let s2 = b"TTTT";
        // ED = 4, but bound = 2 → should return > 2
        let ed = fast_edit_distance_bounded(s1, s2, 2);
        assert!(ed > 2);
    }

    #[test]
    fn test_find_split_positions_even() {
        let positions = find_split_positions_even(0, 100, 4);
        assert_eq!(positions, vec![25, 50, 75]);
    }

    #[test]
    fn test_split_multiplets_realistic() {
        // Simulate a 170bp period array with 10 copies
        // Missing 2 anchors → creates a triplet
        let period = 170;
        let mut seq = vec![b'A'; period * 10];
        // Create some variation at period boundaries
        for i in 0..10 {
            let pos = i * period;
            if pos + 5 < seq.len() {
                seq[pos] = b'T';
                seq[pos + 1] = b'T';
                seq[pos + 2] = b'A';
                seq[pos + 3] = b'G';
                seq[pos + 4] = b'G';
            }
        }

        // Anchors at all positions except 3 and 4 (creates a triplet at positions 2-5)
        let anchors: Vec<usize> = (0..10)
            .filter(|&i| i != 3 && i != 4)
            .map(|i| i * period)
            .collect();

        let result = split_multiplets(&seq, &anchors, period);

        // Should split the triplet segment into 3 monomers
        let total_monomers: usize = result.len();
        assert_eq!(total_monomers, 10); // All 10 monomers should be recovered
    }
}
