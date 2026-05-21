//! Recursive HOR (Higher Order Repeat) decomposition
//!
//! Decomposes HOR monomers into their constituent base-level monomers.
//! Uses autocorrelation to detect periodicity within each HOR and recursively
//! splits until reaching base monomers (no detectable periodicity).

use crate::autocorrelation::{find_period_refined, random_expectation, autocorrelation};
use crate::anchor_by_period::find_anchor_by_period_with_fallback;
use crate::multiplet_split::split_multiplets;
use triple_accel::levenshtein_exp;

/// A base-level monomer after recursive decomposition
#[derive(Debug, Clone)]
pub struct BaseMonomer {
    /// The DNA sequence of this monomer
    pub sequence: String,
    /// Index of the parent HOR from primary decomposition (top-level / level-1 HOR)
    pub hor_idx: usize,
    /// Global index within the entire array (0, 1, 2, ...)
    pub global_idx: usize,
    /// Index within the parent HOR (0, 1, 2, ...)
    pub sub_idx: usize,
    /// Recursion depth at which this leaf was emitted (1 = direct child of top HOR, 2 = grandchild, etc.)
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
    /// idx_within_level of the deepest HOR (in hors.tsv) enclosing this leaf.
    /// At level=1 this equals `hor_idx` (the top-level HOR row).
    /// At level>=2 this is the idx_within_level of the nearest `sub_hor` row.
    pub deepest_hor_idx: usize,
    /// Level of the deepest enclosing HOR (1 = top-level, 2+ = sub-HOR depth).
    pub deepest_hor_level: usize,
}

/// An intermediate HOR detected during recursive decomposition (a non-leaf
/// internal node of the decomposition tree). Top-level HORs (level=1) are
/// already emitted as `monomer` rows in hors.tsv by the primary writer, so
/// this struct only carries level>=2 sub-HORs.
#[derive(Debug, Clone)]
pub struct IntermediateHor {
    /// Depth in the decomposition tree (>= 2; level=1 lives in the primary writer).
    pub level: usize,
    /// Index within this level for the containing array (0-based, fresh per array+level).
    pub idx_within_level: usize,
    /// Index of the parent HOR at level (level - 1). For level=2 this is hor_idx (top-level).
    pub parent_idx: usize,
    /// Period detected for this HOR (its sub-monomers repeat with this length).
    pub period: usize,
    /// Autocorrelation value at the detected period.
    pub autocorr: f64,
    /// Number of sub-monomers this HOR was split into.
    pub n_children: usize,
    /// Length of the HOR sequence in bp.
    pub length: usize,
    /// Source of decomposition: "recursive_anchor" / "recursive_split" / etc.
    pub source: String,
    /// The HOR sequence itself.
    pub sequence: String,
}

/// Result of recursive HOR decomposition
#[derive(Debug, Clone)]
pub struct RecursiveResult {
    /// All base-level monomers from the recursive decomposition (leaves of the tree)
    pub base_monomers: Vec<BaseMonomer>,
    /// All intermediate HOR nodes from the recursive decomposition (level>=2 internal nodes).
    /// Top-level (level=1) HORs are not in this vec — they are written by the primary writer.
    pub intermediate_hors: Vec<IntermediateHor>,
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
    /// Period size classes: "period:count,period:count,..." sorted by count descending
    pub period_classes: String,
}

/// Default minimum submonomer length (bp)
pub const DEFAULT_MIN_SUBMONOMER_LEN: usize = 5;

/// Default autocorrelation threshold for further decomposition
pub const DEFAULT_AUTOCORR_THRESHOLD: f64 = 0.5;

/// Maximum length for edit distance computation to avoid quadratic blowup
const MAX_ED_LEN: usize = 10000;

/// Calculate edit distance between two sequences using SIMD (triple_accel::levenshtein_exp).
/// Returns usize::MAX/2 if either string exceeds MAX_ED_LEN to avoid quadratic blowup.
fn edit_distance(s1: &str, s2: &str) -> usize {
    if s1.len() > MAX_ED_LEN || s2.len() > MAX_ED_LEN {
        return usize::MAX / 2;
    }
    if s1.is_empty() { return s2.len(); }
    if s2.is_empty() { return s1.len(); }
    levenshtein_exp(s1.as_bytes(), s2.as_bytes()) as usize
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
    let mut intermediate_hors: Vec<IntermediateHor> = Vec::new();
    let mut next_idx_per_level: Vec<usize> = Vec::new();
    let mut max_depth: usize = 0;

    for (hor_idx, hor_seq) in hor_monomers.iter().enumerate() {
        // At level=1 the HOR is represented by an existing `monomer` row in
        // hors.tsv with idx=hor_idx, so leaves emitted here without further
        // decomposition still have a row to point to as their deepest HOR.
        let base_monomers = decompose_single_hor(
            hor_seq,
            hor_idx,
            (1, hor_idx),
            min_submonomer_len,
            autocorr_threshold,
            1,
            &mut max_depth,
            &mut intermediate_hors,
            &mut next_idx_per_level,
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

    // Compute period classes: cluster base monomer periods into size classes
    // Monomers within 10% of each other are grouped into the same class
    let period_classes = compute_period_classes(&all_base_monomers);

    RecursiveResult {
        base_monomers: all_base_monomers,
        intermediate_hors,
        max_depth,
        n_expected,
        median_period,
        mean_autocorr,
        consensus_seq,
        iupac_str,
        quality_str,
        period_classes,
    }
}

/// Compute period size classes from base monomers.
/// Groups periods within 10% of each other into clusters.
/// Returns string like "171:85,2420:3" (period:count sorted by count desc).
fn compute_period_classes(base_monomers: &[BaseMonomer]) -> String {
    if base_monomers.is_empty() {
        return "-".to_string();
    }

    // Collect all periods (use sequence length as period proxy)
    let mut lengths: Vec<usize> = base_monomers.iter()
        .map(|m| m.sequence.len())
        .collect();
    lengths.sort();

    // Cluster: greedy merge of sorted lengths within 10% of cluster median
    let mut clusters: Vec<(usize, usize)> = Vec::new(); // (representative_period, count)
    let mut cluster_start = 0;
    while cluster_start < lengths.len() {
        let repr = lengths[cluster_start];
        let threshold = (repr as f64 * 0.1).max(5.0) as usize; // 10% or at least 5bp
        let mut cluster_end = cluster_start + 1;
        while cluster_end < lengths.len() && lengths[cluster_end] <= repr + threshold {
            cluster_end += 1;
        }
        // Representative = median of cluster
        let mid = cluster_start + (cluster_end - cluster_start) / 2;
        let count = cluster_end - cluster_start;
        clusters.push((lengths[mid], count));
        cluster_start = cluster_end;
    }

    // Sort by count descending
    clusters.sort_by(|a, b| b.1.cmp(&a.1));

    // Format as "period:count,period:count,..."
    clusters.iter()
        .map(|(p, c)| format!("{}:{}", p, c))
        .collect::<Vec<_>>()
        .join(",")
}

/// Decompose a single HOR monomer recursively into base-level leaves.
///
/// `parent_hor_row` is `(level, idx_within_level)` of the deepest row in
/// `hors.tsv` that encloses `hor_seq`. At the top-level call (current_level=1)
/// it is `(1, top_hor_idx)` — the existing `monomer` row in hors.tsv that
/// represents `hor_seq` itself, so leaves emitted without further
/// decomposition still have a row to point to.
///
/// On every recursion entry where periodicity is detected and decomposition
/// proceeds, an `IntermediateHor` is pushed (only at current_level >= 2 —
/// the level=1 HOR is already in hors.tsv from the primary writer).
fn decompose_single_hor(
    hor_seq: &str,
    top_hor_idx: usize,
    parent_hor_row: (usize, usize),
    min_len: usize,
    autocorr_threshold: f64,
    current_level: usize,
    max_depth: &mut usize,
    intermediate_hors: &mut Vec<IntermediateHor>,
    next_idx_per_level: &mut Vec<usize>,
) -> Vec<BaseMonomer> {
    if current_level > *max_depth {
        *max_depth = current_level;
    }

    let make_leaf = |sequence: String,
                     sub_idx: usize,
                     period: usize,
                     autocorr: f64,
                     source: String,
                     deepest: (usize, usize)| -> BaseMonomer {
        BaseMonomer {
            sequence,
            hor_idx: top_hor_idx,
            global_idx: 0, // assigned later
            sub_idx,
            level: current_level,
            period,
            autocorr,
            source,
            ed_tmpl: None,
            ed_prev: None,
            ed_next: None,
            ed_per_bp: 0.0,
            cv: 0.0,
            deepest_hor_idx: deepest.1,
            deepest_hor_level: deepest.0,
        }
    };

    // Too short to decompose further — the monomer itself is a leaf.
    if hor_seq.len() < min_len * 2 {
        return vec![make_leaf(
            hor_seq.to_string(),
            0,
            hor_seq.len(),
            0.0,
            "base".to_string(),
            parent_hor_row,
        )];
    }

    let seq_bytes = hor_seq.as_bytes();
    let min_period = min_len;
    let max_period = hor_seq.len() / 2; // need at least 2 copies

    let period_result = find_period_refined(seq_bytes, min_period, max_period);

    match period_result {
        Some((period, autocorr, _excess)) => {
            if autocorr <= autocorr_threshold {
                // Weak periodicity — hor_seq is a leaf at this call.
                return vec![make_leaf(
                    hor_seq.to_string(),
                    0,
                    hor_seq.len(),
                    autocorr,
                    "base".to_string(),
                    parent_hor_row,
                )];
            }

            // Decompose further.
            let sub_monomers = decompose_hor_by_period(hor_seq, period, autocorr);
            let n_children = sub_monomers.len();

            // Reserve an `IntermediateHor` row for hor_seq at level >= 2.
            // Level=1 HORs are already in hors.tsv as `monomer` rows, so the
            // existing `parent_hor_row` (1, top_hor_idx) is reused.
            let this_hor_row: (usize, usize) = if current_level == 1 {
                parent_hor_row
            } else {
                while next_idx_per_level.len() <= current_level {
                    next_idx_per_level.push(0);
                }
                let idx = next_idx_per_level[current_level];
                next_idx_per_level[current_level] += 1;
                let source = sub_monomers
                    .first()
                    .map(|(_, s)| s.clone())
                    .unwrap_or_else(|| "recursive".to_string());
                intermediate_hors.push(IntermediateHor {
                    level: current_level,
                    idx_within_level: idx,
                    parent_idx: parent_hor_row.1,
                    period,
                    autocorr,
                    n_children,
                    length: hor_seq.len(),
                    source,
                    sequence: hor_seq.to_string(),
                });
                (current_level, idx)
            };

            let mut result: Vec<BaseMonomer> = Vec::new();
            for (sub_idx, (sub_seq, sub_source)) in sub_monomers.iter().enumerate() {
                if sub_seq.len() < min_len {
                    result.push(make_leaf(
                        sub_seq.clone(),
                        sub_idx,
                        period,
                        autocorr,
                        sub_source.clone(),
                        this_hor_row,
                    ));
                } else {
                    let sub_bytes = sub_seq.as_bytes();
                    let sub_min_period = min_len;
                    let sub_max_period = sub_seq.len() / 2;

                    if sub_max_period >= sub_min_period {
                        if let Some((_, sub_autocorr, _)) =
                            find_period_refined(sub_bytes, sub_min_period, sub_max_period)
                        {
                            if sub_autocorr > autocorr_threshold {
                                let deeper = decompose_single_hor(
                                    sub_seq,
                                    top_hor_idx,
                                    this_hor_row,
                                    min_len,
                                    autocorr_threshold,
                                    current_level + 1,
                                    max_depth,
                                    intermediate_hors,
                                    next_idx_per_level,
                                );
                                for mut m in deeper {
                                    m.sub_idx = sub_idx * 1000 + m.sub_idx;
                                    result.push(m);
                                }
                                continue;
                            }
                        }
                    }

                    let final_autocorr = if sub_max_period >= sub_min_period {
                        find_period_refined(sub_bytes, sub_min_period, sub_max_period)
                            .map(|(_, ac, _)| ac)
                            .unwrap_or(0.0)
                    } else {
                        0.0
                    };

                    result.push(make_leaf(
                        sub_seq.clone(),
                        sub_idx,
                        period,
                        final_autocorr,
                        sub_source.clone(),
                        this_hor_row,
                    ));
                }
            }
            result
        }
        None => {
            // No periodicity at all — hor_seq is a leaf at this call.
            vec![make_leaf(
                hor_seq.to_string(),
                0,
                hor_seq.len(),
                compute_max_autocorr(seq_bytes, min_len),
                "base".to_string(),
                parent_hor_row,
            )]
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
