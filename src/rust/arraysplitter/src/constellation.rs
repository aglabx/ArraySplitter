//! Per-level Homogeneity/Heterogeneity Regime Detector & Router, Indel-Tolerant Substring
//! ($k$-mer) Lag Autocorrelation, and Circularly Phased Substring Constellations.
//!
//! At each decomposition level (top-level array or recursive HOR unit), satellite DNA
//! sequences are classified into one of three regimes:
//! 1. [`LevelRegime::HomogeneousSingleAnchor`] / [`LevelRegime::HomogeneousLongHorAnchor`]:
//!    Fixed-lag nucleotide autocorrelation is strong and a single exact anchor $k$-mer
//!    ($k \ge 7$, or $k \in \{18, 21\}$ for multi-kilobase HORs) covers $\ge 50\%$ of copies
//!    with high regularity. Routed to the classic v1.8.2 single-anchor + `split_multiplets` path.
//! 2. [`LevelRegime::HeterogeneousConstellation`]:
//!    Monomers diverge by 20–45% with $\pm 1\text{..}4\text{ bp}$ indel jitter so fixed-lag
//!    nucleotide autocorrelation is weak (`0.28–0.34`) and no single 7-mer reaches 50%
//!    monomer coverage, **while** windowed $k$-mer ($k=6..7$, $\pm 4\%$ indel window,
//!    harmonic $D_w(d)+0.5 D_w(2d)$) lag density has $\text{SNR} \ge 3.5$ and a circularly
//!    phased constellation of periodic $k$-mers on $\mathbb{Z}/P\mathbb{Z}$ cuts $\ge 45\%$
//!    of expected copies into regular $1\times P$ intervals.
//! 3. [`LevelRegime::TerminalLeaf`]:
//!    Neither a homogeneous single anchor nor a verified phased constellation exists at any
//!    sub-period; recursion stops cleanly without blind `even_split`.

use std::collections::HashMap;
use crate::sequence::decode_kmer;

/// Maximum sequence bytes sampled per window when building lag histograms on mega-arrays.
const MAX_LAG_SAMPLE_BYTES: usize = 120_000;

/// Maximum occurrences per distinct k-mer in a single window to avoid quadratic blowup.
const MAX_OCC_PER_KMER: usize = 450;

/// Maximum number of periodic k-mers included in a phased constellation.
const MAX_CONSTELLATION_KMERS: usize = 80;

/// Classification of a sequence or HOR ensemble at a candidate decomposition level.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LevelRegime {
    /// Homogeneous repeat decomposable by a single exact k-mer anchor (k in 5..=15).
    HomogeneousSingleAnchor,
    /// Homogeneous multi-kilobase HOR (period >= 800 bp) decomposable by an exact 18..21-mer.
    HomogeneousLongHorAnchor,
    /// Heterogeneous repeat of diverged/indel-jittered units sharing a circularly phased
    /// constellation of substrings on Z / P Z.
    HeterogeneousConstellation,
    /// Terminal base monomer leaf (no supported sub-periodicity).
    TerminalLeaf,
}

impl LevelRegime {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::HomogeneousSingleAnchor => "homogeneous_anchor",
            Self::HomogeneousLongHorAnchor => "homogeneous_hor_anchor",
            Self::HeterogeneousConstellation => "heterogeneous_constellation",
            Self::TerminalLeaf => "terminal_leaf",
        }
    }
}

/// Encode ACGT base to 0..3.
#[inline]
fn encode_base(b: u8) -> Option<u32> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Check if a 2-bit encoded k-mer of length `k` is low-complexity:
/// - <= 2 distinct nucleotides (e.g. AAAAAAA, TTTTTTT, ACACACA, CCCCGGG), OR
/// - contains a homopolymer run of length >= 5.
fn is_low_complexity_code(code: u32, k: usize) -> bool {
    let mut mask = 0u8;
    let mut prev = 4u32;
    let mut run = 0usize;
    let mut max_run = 0usize;
    for i in 0..k {
        let b = (code >> (2 * (k - 1 - i))) & 3;
        mask |= 1u8 << b;
        if b == prev {
            run += 1;
        } else {
            prev = b;
            run = 1;
        }
        if run > max_run {
            max_run = run;
        }
    }
    mask.count_ones() <= 2 || max_run >= 5
}

/// Precompute low-complexity table for all $4^k$ codes ($k \le 8$).
fn build_low_complexity_table(k: usize) -> Vec<bool> {
    let n_codes = 1usize << (2 * k);
    (0..n_codes)
        .map(|c| is_low_complexity_code(c as u32, k))
        .collect()
}

/// Sample up to 3 slices totaling at most `MAX_LAG_SAMPLE_BYTES` from a long sequence.
fn sample_windows<'a>(seq: &'a [u8]) -> Vec<&'a [u8]> {
    if seq.len() <= MAX_LAG_SAMPLE_BYTES {
        vec![seq]
    } else {
        let w = MAX_LAG_SAMPLE_BYTES / 3;
        let mid = seq.len() / 2 - w / 2;
        vec![
            &seq[..w],
            &seq[mid..mid + w],
            &seq[seq.len() - w..],
        ]
    }
}

/// Candidate period detected by windowed substring lag autocorrelation.
#[derive(Debug, Clone, Copy)]
pub struct KmerPeriodCandidate {
    /// Fundamental period in bp
    pub period: usize,
    /// Windowed lag density at `period`
    pub primary_density: f64,
    /// Harmonic-reinforced score $D_w(P) + 0.5 D_w(2P)$
    pub harmonic_score: f64,
    /// Signal-to-noise ratio over background median density
    pub snr: f64,
}

/// Detect a heterogeneous monomer sub-period inside a single consensus HOR sequence
/// using 6-mer intra-unit lag density and distinct-kmer diversity.
///
/// Evaluating on the consensus sequence of a parent unit (rather than summing across
/// hundreds of identical copies of a homogeneous monomer) ensures that accidental single-copy
/// internal duplications in homogeneous monomers (e.g., 716-bp or 1249-bp zebra finch monomers)
/// are not multiplied by array copy number, while true heterogeneous HORs (composed of
/// $m \ge 5$ diverged monomers per HOR unit) exhibit strong multi-harmonic 6-mer lag density
/// (`harmonic_score >= 12.0`) supported by $\ge 50$ distinct repeating 6-mers.
pub fn find_consensus_subperiod_kmer(
    consensus_seq: &[u8],
    min_d: usize,
    max_d: usize,
) -> Option<KmerPeriodCandidate> {
    let len = consensus_seq.len();
    let min_d = min_d.max(50);
    let max_d = max_d.min(len / 3).min(600);
    if len < 300 || max_d < min_d {
        return None;
    }

    let k = 6usize;
    let n_codes = 1usize << (2 * k);
    let kmask = (n_codes as u32) - 1;
    // Exclude only pure 1-base homopolymers (AAAAAA, TTTTTT, etc.) so AT/GC-tract
    // satellite motifs with 2 distinct bases are preserved.
    let is_homopolymer = |c: u32| -> bool {
        let b0 = c & 3;
        let mut all_same = true;
        for i in 1..k {
            if ((c >> (2 * i)) & 3) != b0 {
                all_same = false;
                break;
            }
        }
        all_same
    };

    let mut buckets: Vec<Vec<u32>> = vec![Vec::new(); n_codes];
    let mut active_codes: Vec<usize> = Vec::new();
    let mut code = 0u32;
    let mut valid_len = 0usize;

    for (i, &b) in consensus_seq.iter().enumerate() {
        if let Some(v) = encode_base(b) {
            code = ((code << 2) | v) & kmask;
            valid_len += 1;
            if valid_len >= k {
                let idx = code as usize;
                if !is_homopolymer(code) {
                    if buckets[idx].is_empty() {
                        active_codes.push(idx);
                    }
                    if buckets[idx].len() < 200 {
                        buckets[idx].push((i + 1 - k) as u32);
                    }
                }
            }
        } else {
            code = 0;
            valid_len = 0;
        }
    }

    let mut lags = vec![0u32; len + 1];
    let mut distinct_at_lag: Vec< Vec<u16> > = vec![Vec::new(); len + 1];

    for &idx in &active_codes {
        let ps = &buckets[idx];
        if ps.len() < 2 || ps.len() > 200 {
            continue;
        }
        for a in 0..ps.len() {
            let pa = ps[a] as usize;
            for b in (a + 1)..ps.len() {
                let d = (ps[b] as usize) - pa;
                if d <= len {
                    lags[d] += 1;
                    distinct_at_lag[d].push(idx as u16);
                }
            }
        }
    }

    let mut dw = vec![0.0f64; len + 1];
    let mut raw = vec![0u32; len + 1];
    let mut u_cnt = vec![0usize; max_d + 1];
    let mut seen_stamp = vec![0u32; n_codes];
    let mut stamp = 1u32;

    for d in 20..len {
        let w_half = ((d as f64 * 0.04).round() as usize).max(2);
        let lo = d.saturating_sub(w_half).max(1);
        let hi = (d + w_half).min(len);
        let mut s = 0u32;
        for dd in lo..=hi {
            s += lags[dd];
        }
        raw[d] = s;
        dw[d] = s as f64 / ((hi - lo + 1) as f64);

        if d >= min_d && d <= max_d {
            stamp += 1;
            let mut u = 0usize;
            for dd in lo..=hi {
                for &c in &distinct_at_lag[dd] {
                    let ci = c as usize;
                    if seen_stamp[ci] != stamp {
                        seen_stamp[ci] = stamp;
                        u += 1;
                    }
                }
            }
            u_cnt[d] = u;
        }
    }

    let mut best_d = 0usize;
    let mut best_h = 0.0f64;
    for d in min_d..=max_d {
        let mut h = dw[d];
        if 2 * d < len {
            h += 0.50 * dw[2 * d];
        }
        if 3 * d < len {
            h += 0.33 * dw[3 * d];
        }
        if 4 * d < len {
            h += 0.25 * dw[4 * d];
        }
        if h > best_h {
            best_h = h;
            best_d = d;
        }
    }

    if best_d == 0 {
        return None;
    }

    // Refine to fundamental if best_d is an integer multiple of d >= min_d with dw[d] >= 0.70 * dw[best_d]
    for d in min_d..best_d {
        let ratio = (best_d as f64 / d as f64).round();
        if ratio >= 2.0 {
            let approx = (ratio as usize) * d;
            let tol = ((best_d as f64 * 0.05).round() as usize).max(3);
            if best_d.abs_diff(approx) <= tol && dw[d] >= 0.70 * dw[best_d] && u_cnt[d] >= 50 {
                best_d = d;
                best_h = dw[d]
                    + if 2 * d < len { 0.50 * dw[2 * d] } else { 0.0 }
                    + if 3 * d < len { 0.33 * dw[3 * d] } else { 0.0 }
                    + if 4 * d < len { 0.25 * dw[4 * d] } else { 0.0 };
                break;
            }
        }
    }

    // Center best_d on highest fixed ±2 bp lag sum within ±3% (avoids window-width step jump at d=113)
    let w_center = ((best_d as f64 * 0.03).round() as usize).max(2);
    let lo_c = best_d.saturating_sub(w_center).max(min_d).max(3);
    let hi_c = (best_d + w_center).min(max_d).min(len.saturating_sub(3));
    if let Some(centered) = (lo_c..=hi_c).max_by_key(|&d| {
        lags[d - 2] + lags[d - 1] + lags[d] + lags[d + 1] + lags[d + 2]
    }) {
        best_d = centered;
    }

    let bg_sum: f64 = dw[min_d..=max_d].iter().sum();
    let bg_mean = (bg_sum / ((max_d - min_d + 1) as f64)).max(1e-5);
    let snr = best_h / bg_mean;

    // Strict intra-unit heterogeneity gate:
    // Requires strong multi-harmonic lag density (best_h >= 9.5), high 6-mer diversity
    // (u_cnt >= 50 distinct 6-mers at lag best_d), and raw lag hits >= 55.
    if best_h < 9.5 || u_cnt[best_d] < 50 || raw[best_d] < 55 || snr < 2.8 {
        return None;
    }

    Some(KmerPeriodCandidate {
        period: best_d,
        primary_density: dw[best_d],
        harmonic_score: best_h,
        snr,
    })
}

/// Detect the fundamental monomer period on a single sequence using windowed $k$-mer
/// lag autocorrelation.
pub fn find_monomer_period_kmer(
    seq: &[u8],
    min_d: usize,
    max_d: usize,
) -> Option<KmerPeriodCandidate> {
    let windows = sample_windows(seq);
    find_subperiod_kmer_ensemble(&windows, min_d, max_d)
}

/// Detect the fundamental repeat period across an ensemble of sequences (e.g. HOR units
/// or windows of an array) by accumulating strictly intra-sequence $k$-mer pair lags.
pub fn find_subperiod_kmer_ensemble(
    seqs: &[&[u8]],
    min_d: usize,
    max_d: usize,
) -> Option<KmerPeriodCandidate> {
    let min_d = min_d.max(15);
    let total_len: usize = seqs.iter().map(|s| s.len()).sum();
    if max_d < min_d || total_len < min_d * 3 {
        return None;
    }

    let k: usize = if total_len < 3000 || max_d < 80 { 6 } else { 7 };
    let n_codes = 1usize << (2 * k);
    let kmask = (n_codes as u32) - 1;
    let low_comp = build_low_complexity_table(k);

    let max_seq_len = seqs.iter().map(|s| s.len()).max().unwrap_or(0);
    let max_lag = (max_d * 2 + 20).min(max_seq_len);
    if max_lag <= min_d {
        return None;
    }

    // If total_len across many HOR units exceeds 360 kb, subsample evenly
    let sampled_seqs: Vec<&[u8]> = if total_len > 360_000 && seqs.len() > 36 {
        let step = (seqs.len() / 36).max(1);
        seqs.iter().step_by(step).copied().collect()
    } else {
        seqs.to_vec()
    };

    let mut lag_counts = vec![0u64; max_lag + 1];
    let mut buckets: Vec<Vec<u32>> = vec![Vec::new(); n_codes];
    let mut active_codes: Vec<usize> = Vec::new();
    let mut sampled_window_lens: Vec<usize> = Vec::new();

    for &seq in &sampled_seqs {
        if seq.len() <= min_d + k {
            continue;
        }
        for win in sample_windows(seq) {
            if win.len() <= min_d + k {
                continue;
            }
            sampled_window_lens.push(win.len());
            for &c in &active_codes {
                buckets[c].clear();
            }
            active_codes.clear();

            let mut code = 0u32;
            let mut valid_len = 0usize;
            for (i, &b) in win.iter().enumerate() {
                if let Some(v) = encode_base(b) {
                    code = ((code << 2) | v) & kmask;
                    valid_len += 1;
                    if valid_len >= k {
                        let idx = code as usize;
                        if !low_comp[idx] {
                            let pos = (i + 1 - k) as u32;
                            if buckets[idx].is_empty() {
                                active_codes.push(idx);
                            }
                            if buckets[idx].len() < MAX_OCC_PER_KMER {
                                buckets[idx].push(pos);
                            }
                        }
                    }
                } else {
                    code = 0;
                    valid_len = 0;
                }
            }

            for &idx in &active_codes {
                let pos_list = &buckets[idx];
                let m = pos_list.len();
                if m < 2 {
                    continue;
                }
                for a in 0..m {
                    let pa = pos_list[a] as usize;
                    for b in (a + 1)..m {
                        let d = (pos_list[b] as usize) - pa;
                        if d > max_lag {
                            break;
                        }
                        lag_counts[d] += 1;
                    }
                }
            }
        }
    }

    if sampled_window_lens.is_empty() {
        return None;
    }

    // Prefix sum of lag_counts for O(1) window queries
    let mut prefix = vec![0u64; max_lag + 2];
    for d in 1..=max_lag {
        prefix[d + 1] = prefix[d] + lag_counts[d];
    }

    let pairs_at_lag = |d: usize| -> usize {
        sampled_window_lens
            .iter()
            .filter(|&&l| l > d)
            .map(|&l| l - d)
            .sum()
    };

    let windowed_density = |d: usize| -> (f64, u64) {
        if d == 0 || d > max_lag {
            return (0.0, 0);
        }
        let n_pairs = pairs_at_lag(d);
        if n_pairs < 30 {
            return (0.0, 0);
        }
        let w = ((d as f64 * 0.04).round() as usize).clamp(2, 6);
        let lo = d.saturating_sub(w).max(1);
        let hi = (d + w).min(max_lag);
        let sum = prefix[hi + 1] - prefix[lo];
        (sum as f64 / n_pairs as f64, sum)
    };

    let mut primary = vec![0.0f64; max_d + 1];
    let mut raw_hits = vec![0u64; max_d + 1];
    let mut scores = vec![0.0f64; max_d + 1];
    let mut non_zero_densities = Vec::with_capacity(max_d - min_d + 1);

    for d in min_d..=max_d {
        let (dens, hits) = windowed_density(d);
        primary[d] = dens;
        raw_hits[d] = hits;
        if dens > 0.0 {
            non_zero_densities.push(dens);
        }
        let (dens2, _) = windowed_density(2 * d);
        scores[d] = dens + 0.5 * dens2;
    }

    if non_zero_densities.is_empty() {
        return None;
    }
    non_zero_densities.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let bg_median = non_zero_densities[non_zero_densities.len() / 2].max(1e-5);

    // Collect local peaks of scores[d] that exceed Poisson random expectation by >= 3.5x and +12 hits
    let random_pair_prob = 1.0 / (n_codes as f64);
    let mut peaks: Vec<usize> = Vec::new();
    for d in (min_d + 1)..max_d {
        let w = ((d as f64 * 0.04).round() as usize).clamp(2, 6);
        let win_width = (2 * w + 1) as f64;
        let exp_random_hits = (pairs_at_lag(d) as f64) * win_width * random_pair_prob;
        let min_required_hits = ((exp_random_hits * 3.5).ceil() as u64)
            .max((exp_random_hits as u64) + 12)
            .max(15);
        if scores[d] >= scores[d - 1]
            && scores[d] >= scores[d + 1]
            && raw_hits[d] >= min_required_hits
        {
            peaks.push(d);
        }
    }
    if peaks.is_empty() {
        return None;
    }

    peaks.sort_by(|&a, &b| {
        scores[b]
            .partial_cmp(&scores[a])
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.cmp(&b))
    });

    let center_peak = |d0: usize| -> usize {
        let w = ((d0 as f64 * 0.03).round() as usize).max(2);
        let lo = d0.saturating_sub(w).max(min_d);
        let hi = (d0 + w).min(max_d);
        (lo..=hi)
            .max_by(|&a, &b| {
                primary[a]
                    .partial_cmp(&primary[b])
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or(d0)
    };

    let mut best_d = center_peak(peaks[0]);

    // Check if best_d is a 2x..8x multiple of a smaller fundamental period P = best_d / m
    let orig_best_d = best_d;
    for m in (2..=8usize).rev() {
        let sub = orig_best_d / m;
        if sub >= min_d.max(50) {
            let sub_centered = center_peak(sub);
            let sub_primary = primary[sub_centered];
            let best_primary = primary[orig_best_d];
            if sub_primary >= 0.50 * best_primary
                && scores[sub_centered] / bg_median >= 3.5
            {
                best_d = sub_centered;
                break;
            }
        }
    }

    // Conversely, check if best_d is a half-period P/2 inflated by 0.5 * D(2d)
    if 2 * best_d <= max_d {
        let d2 = center_peak(2 * best_d);
        let (d3, _) = windowed_density(3 * best_d);
        if primary[d2] > 1.75 * primary[best_d] && d3 < 0.35 * primary[d2] {
            best_d = d2;
        }
    }

    let snr = scores[best_d] / bg_median;
    if snr < 3.5 || primary[best_d] / bg_median < 2.2 {
        return None;
    }

    Some(KmerPeriodCandidate {
        period: best_d,
        primary_density: primary[best_d],
        harmonic_score: scores[best_d],
        snr,
    })
}

/// Detect a long-range HOR period $P_{\text{HOR}} \in [\max(300, 3 P_{\text{mono}}), \text{max\_d}]$
/// using exact 10-mer lag density ($\pm 2\text{ bp}$ window).
pub fn find_hor_period_kmer(
    seq: &[u8],
    mono_period: usize,
    max_d: usize,
) -> Option<(usize, f64)> {
    let min_hor = (mono_period * 3).max(300);
    let max_hor = max_d.min(seq.len() / 3);
    if max_hor <= min_hor || seq.len() < min_hor * 3 {
        return None;
    }

    let k = 10usize;
    let kmask = (1u32 << (2 * k)) - 1;
    let mut kmer_pos: HashMap<u32, Vec<u32>> = HashMap::new();

    let mut code = 0u32;
    let mut valid_len = 0usize;
    for (i, &b) in seq.iter().enumerate() {
        if let Some(v) = encode_base(b) {
            code = ((code << 2) | v) & kmask;
            valid_len += 1;
            if valid_len >= k && !is_low_complexity_code(code, k) {
                let entry = kmer_pos.entry(code).or_default();
                if entry.len() < 300 {
                    entry.push((i + 1 - k) as u32);
                }
            }
        } else {
            code = 0;
            valid_len = 0;
        }
    }

    let mut lag_counts = vec![0u32; max_hor + 5];
    for pos_list in kmer_pos.values() {
        let m = pos_list.len();
        if m < 2 || m > 250 {
            continue;
        }
        for a in 0..m {
            let pa = pos_list[a] as usize;
            for b in (a + 1)..m {
                let d = (pos_list[b] as usize) - pa;
                if d > max_hor + 2 {
                    break;
                }
                if d >= min_hor.saturating_sub(2) {
                    lag_counts[d] += 1;
                }
            }
        }
    }

    let mut best_d = 0usize;
    let mut best_density = 0.0f64;
    let mut best_hits = 0u32;
    let mut densities = Vec::new();

    for d in min_hor..=max_hor {
        let n_pairs = seq.len() - d;
        if n_pairs < 500 {
            continue;
        }
        let hits: u32 = (d.saturating_sub(2)..=(d + 2).min(max_hor + 2))
            .map(|x| lag_counts[x])
            .sum();
        let dens = hits as f64 / n_pairs as f64;
        if dens > 0.0 {
            densities.push(dens);
        }
        if dens > best_density {
            best_density = dens;
            best_d = d;
            best_hits = hits;
        }
    }

    if best_d == 0 || best_hits < 60 || densities.is_empty() {
        return None;
    }

    densities.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let bg = densities[densities.len() / 2].max(1e-5);
    let snr = best_density / bg;

    if snr < 12.0 || best_density < 0.008 {
        return None;
    }

    // Check if best_d is a 2x or 3x multiple of a smaller HOR period
    for div in [3usize, 2] {
        let sub = best_d / div;
        if sub >= min_hor {
            let sub_hits: u32 = ((sub - 2)..=(sub + 2))
                .map(|x| lag_counts[x])
                .sum();
            let sub_dens = sub_hits as f64 / (seq.len() - sub) as f64;
            if sub_dens >= 0.65 * best_density {
                best_d = sub;
                best_density = sub_dens;
            }
        }
    }

    Some((best_d, best_density))
}

/// Find an exact long anchor ($k \in [18, 21]$) for cutting multi-kilobase HORs
/// where 15-mers may repeat inside the HOR due to internal monomer structure.
pub fn find_long_hor_anchor(
    seq: &[u8],
    period: usize,
) -> Option<(String, Vec<usize>, f64)> {
    if period < 800 || seq.len() < period * 2 {
        return None;
    }
    let expected = (seq.len() / period).max(2);
    let min_freq = (expected / 2).max(2);
    let max_freq = expected * 2;
    let min_spacing = (period as f64 * 0.70) as usize;
    let tol = (period as f64 * 0.15).max(10.0) as usize;

    for &k in &[21usize, 18] {
        if seq.len() <= k {
            continue;
        }
        let mut table: HashMap<&[u8], Vec<usize>> = HashMap::new();
        for i in 0..=(seq.len() - k) {
            let slice = &seq[i..i + k];
            if slice.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')) {
                table.entry(slice).or_default().push(i);
            }
        }

        let mut best: Option<(&[u8], Vec<usize>, f64)> = None;
        for (kmer_bytes, positions) in table {
            let freq = positions.len();
            if freq < min_freq || freq > max_freq {
                continue;
            }
            let n_spacings = freq - 1;
            let mut n_close = 0usize;
            let mut n_good = 0usize;
            for w in positions.windows(2) {
                let sp = w[1] - w[0];
                if sp < min_spacing {
                    n_close += 1;
                }
                let dev = if sp > period { sp - period } else { period - sp };
                if dev <= tol {
                    n_good += 1;
                } else {
                    let dev2 = if sp > 2 * period { sp - 2 * period } else { 2 * period - sp };
                    if dev2 <= tol * 2 {
                        n_good += 1;
                    }
                }
            }
            let uniq = 1.0 - (n_close as f64 / n_spacings as f64);
            let reg = n_good as f64 / n_spacings as f64;
            if uniq < 0.85 || reg < 0.65 {
                continue;
            }
            let cov = (freq as f64 / expected as f64).min(1.0);
            let score = uniq * reg * cov;
            let is_better = match &best {
                None => true,
                Some((best_kmer, _, best_score)) => {
                    score > *best_score + 1e-9
                        || ((score - *best_score).abs() <= 1e-9 && kmer_bytes < *best_kmer)
                }
            };
            if is_better {
                best = Some((kmer_bytes, positions, score));
            }
        }

        if let Some((kmer_bytes, positions, score)) = best {
            let kmer_str = String::from_utf8_lossy(kmer_bytes).to_ascii_uppercase();
            return Some((kmer_str, positions, score));
        }
    }
    None
}

/// A circularly phased constellation of periodic $k$-mers on $\mathbb{Z}/P\mathbb{Z}$.
#[derive(Debug, Clone)]
pub struct PhasedConstellation {
    /// Target monomer period $P$
    pub period: usize,
    /// $k$-mer length ($k \in [6, 8]$)
    pub k: usize,
    /// Lookup table of size $4^k$: `Some((phase, weight))` for each phased $k$-mer
    phased_table: Vec<Option<(usize, f64)>>,
    /// Number of distinct $k$-mers phased into the constellation
    num_phased: usize,
    /// Representative seed $k$-mer string (at phase $\phi \approx 0$)
    pub seed_kmer: String,
    /// Top single-kmer coverage $C_1(P)$ observed during construction
    pub top1_coverage: f64,
}

impl PhasedConstellation {
    /// Number of phased $k$-mers in the constellation.
    pub fn num_phased(&self) -> usize {
        self.num_phased
    }

    /// Build a circularly phased constellation of periodic $k$-mers from one or more
    /// sequences sharing monomer period `period`.
    pub fn build(
        seqs: &[&[u8]],
        period: usize,
        align_to_sequence_starts: bool,
    ) -> Option<Self> {
        if period < 15 || period > 1200 {
            return None;
        }
        let total_len: usize = seqs.iter().map(|s| s.len()).sum();
        let n_exp = (total_len as f64 / period as f64).max(1.0);
        if n_exp < 4.0 {
            return None;
        }

        let k: usize = if n_exp < 25.0 || period < 50 {
            6
        } else if period > 250 {
            8
        } else {
            7
        };
        let n_codes = 1usize << (2 * k);
        let kmask = (n_codes as u32) - 1;
        let low_comp = build_low_complexity_table(k);

        let sampled_seqs: Vec<&[u8]> = if total_len > 400_000 && seqs.len() > 40 {
            let step = (seqs.len() / 40).max(1);
            seqs.iter().step_by(step).copied().collect()
        } else if total_len > 400_000 && seqs.len() == 1 {
            sample_windows(seqs[0])
        } else {
            seqs.to_vec()
        };

        let sampled_len: usize = sampled_seqs.iter().map(|s| s.len()).sum();
        let sampled_exp = (sampled_len as f64 / period as f64).max(1.0);

        let mut buckets: Vec<Vec<(u16, u32)>> = vec![Vec::new(); n_codes];
        for (s_idx, &seq) in sampled_seqs.iter().enumerate() {
            if seq.len() < k {
                continue;
            }
            let mut code = 0u32;
            let mut valid_len = 0usize;
            for (i, &b) in seq.iter().enumerate() {
                if let Some(v) = encode_base(b) {
                    code = ((code << 2) | v) & kmask;
                    valid_len += 1;
                    if valid_len >= k {
                        let idx = code as usize;
                        if !low_comp[idx] {
                            let pos = (i + 1 - k) as u32;
                            buckets[idx].push((s_idx as u16, pos));
                        }
                    }
                } else {
                    code = 0;
                    valid_len = 0;
                }
            }
        }

        let min_count = ((sampled_exp * 0.04).floor() as usize).max(2);
        let max_count = (sampled_exp * 1.40).ceil() as usize;
        let close_thresh = (period as f64 * 0.70) as u32;

        #[derive(Clone)]
        struct Cand {
            code: usize,
            count: usize,
            score: f64,
            per_ratio: f64,
            uniq: f64,
        }
        let mut candidates: Vec<Cand> = Vec::new();

        for (code, occs) in buckets.iter().enumerate() {
            let c = occs.len();
            if c < min_count || c > max_count {
                continue;
            }
            let mut n_close = 0usize;
            let mut good_lags = 0usize;
            let mut total_lags = 0usize;

            for a in 0..c {
                let (sa, pa) = occs[a];
                if a + 1 < c && occs[a + 1].0 == sa && (occs[a + 1].1 - pa) < close_thresh {
                    n_close += 1;
                }
                for b in (a + 1)..c {
                    let (sb, pb) = occs[b];
                    if sb != sa {
                        break;
                    }
                    let d = (pb - pa) as usize;
                    if d > (period * 11) / 2 {
                        break;
                    }
                    total_lags += 1;
                    for mult in 1..=5 {
                        let target = mult * period;
                        let tol = ((target as f64 * 0.065).round() as usize).max(3);
                        let dev = if d > target { d - target } else { target - d };
                        if dev <= tol {
                            good_lags += 1;
                            break;
                        }
                    }
                }
            }

            if total_lags < 2 || good_lags < 2 {
                continue;
            }
            let uniq = 1.0 - (n_close as f64 / c as f64);
            if uniq < 0.72 {
                continue;
            }
            let per_ratio = good_lags as f64 / total_lags as f64;
            let self_repeat_rate = good_lags as f64 / c as f64;
            if per_ratio < 0.25 || self_repeat_rate < 0.10 {
                continue;
            }
            let score = (c as f64) * per_ratio * self_repeat_rate * uniq * uniq;
            candidates.push(Cand {
                code,
                count: c,
                score,
                per_ratio,
                uniq,
            });
        }

        if candidates.is_empty() {
            return None;
        }

        candidates.sort_by(|a, b| {
            b.score
                .partial_cmp(&a.score)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| a.code.cmp(&b.code))
        });

        let top1_coverage = (candidates[0].count as f64 / sampled_exp).min(1.0);
        candidates.truncate(MAX_CONSTELLATION_KMERS);

        let seed_code = candidates[0].code;
        let seed_kmer = decode_kmer(seed_code as u64, k);

        let mut phased_table: Vec<Option<(usize, f64)>> = vec![None; n_codes];
        let seed_weight = (candidates[0].per_ratio * candidates[0].uniq).clamp(0.2, 1.0);
        phased_table[seed_code] = Some((0, seed_weight));
        let mut num_phased = 1usize;

        let mut ref_by_seq: Vec<Vec<(u32, usize)>> = vec![Vec::new(); sampled_seqs.len()];
        for &(s_idx, pos) in &buckets[seed_code] {
            ref_by_seq[s_idx as usize].push((pos, 0));
        }

        let min_votes = ((sampled_exp * 0.025).floor() as usize).max(2);
        let max_pair_dist = period - 4;

        for _pass in 0..2 {
            for cand in &candidates[1..] {
                if phased_table[cand.code].is_some() {
                    continue;
                }
                let mut votes = vec![0u32; period];
                let mut total_votes = 0u32;

                for &(s_idx, pos) in &buckets[cand.code] {
                    for &(ref_pos, ref_phase) in &ref_by_seq[s_idx as usize] {
                        let diff = (pos as isize) - (ref_pos as isize);
                        if diff.unsigned_abs() <= max_pair_dist {
                            let theta = (diff + (ref_phase as isize)).rem_euclid(period as isize) as usize;
                            votes[theta] += 1;
                            total_votes += 1;
                        }
                    }
                }

                if (total_votes as usize) < min_votes {
                    continue;
                }

                let mut best_center = 0usize;
                let mut best_win_votes = 0u32;
                for center in 0..period {
                    let mut wsum = 0u32;
                    for delta in -3isize..=3 {
                        let idx = (center as isize + delta).rem_euclid(period as isize) as usize;
                        wsum += votes[idx];
                    }
                    if wsum > best_win_votes {
                        best_win_votes = wsum;
                        best_center = center;
                    }
                }

                if (best_win_votes as usize) >= min_votes
                    && (best_win_votes as f64) / (total_votes as f64) >= 0.45
                {
                    let mut weighted_delta = 0isize;
                    for delta in -3isize..=3 {
                        let idx = (best_center as isize + delta).rem_euclid(period as isize) as usize;
                        weighted_delta += delta * (votes[idx] as isize);
                    }
                    let refined_phase = (best_center as isize
                        + (weighted_delta as f64 / best_win_votes as f64).round() as isize)
                        .rem_euclid(period as isize) as usize;

                    let weight = (cand.per_ratio * cand.uniq).clamp(0.2, 1.0);
                    phased_table[cand.code] = Some((refined_phase, weight));
                    num_phased += 1;
                    for &(s_idx, pos) in &buckets[cand.code] {
                        ref_by_seq[s_idx as usize].push((pos, refined_phase));
                    }
                }
            }
        }

        if align_to_sequence_starts && sampled_seqs.len() >= 2 {
            let mut start_votes = vec![0.0f64; period];
            for seq in &sampled_seqs {
                let scan_len = (period * 3).min(seq.len());
                if scan_len < k {
                    continue;
                }
                let mut code = 0u32;
                let mut valid_len = 0usize;
                for (i, &b) in seq[..scan_len].iter().enumerate() {
                    if let Some(v) = encode_base(b) {
                        code = ((code << 2) | v) & kmask;
                        valid_len += 1;
                        if valid_len >= k {
                            if let Some((phi, weight)) = phased_table[code as usize] {
                                let pos = i + 1 - k;
                                let offset = (pos as isize - phi as isize).rem_euclid(period as isize) as usize;
                                start_votes[offset] += weight;
                            }
                        }
                    } else {
                        code = 0;
                        valid_len = 0;
                    }
                }
            }
            let mut best_shift = 0usize;
            let mut best_sum = 0.0f64;
            for c in 0..period {
                let mut s = 0.0f64;
                for d in -2isize..=2 {
                    let idx = (c as isize + d).rem_euclid(period as isize) as usize;
                    s += start_votes[idx];
                }
                if s > best_sum {
                    best_sum = s;
                    best_shift = c;
                }
            }
            if best_sum > 0.0 {
                for slot in phased_table.iter_mut().flatten() {
                    slot.0 = (slot.0 + best_shift) % period;
                }
            }
        }

        Some(Self {
            period,
            k,
            phased_table,
            num_phased,
            seed_kmer,
            top1_coverage,
        })
    }

    /// Project phased $k$-mer occurrences in `seq` onto virtual cut positions and
    /// select the optimal cut lattice via dynamic programming.
    pub fn project_cuts(&self, seq: &[u8]) -> Vec<usize> {
        let n = seq.len();
        if n < self.period || n < self.k {
            return Vec::new();
        }

        let kmask = (1u32 << (2 * self.k)) - 1;
        let min_same_kmer_dist = (self.period * 2) / 3;

        let mut last_pos_for_code: HashMap<u32, usize> = HashMap::new();
        let mut raw_votes = vec![0.0f64; n];

        let mut code = 0u32;
        let mut valid_len = 0usize;
        let mut hits: Vec<(usize, usize, f64)> = Vec::new();

        for (i, &b) in seq.iter().enumerate() {
            if let Some(v) = encode_base(b) {
                code = ((code << 2) | v) & kmask;
                valid_len += 1;
                if valid_len >= self.k {
                    if let Some((phi, weight)) = self.phased_table[code as usize] {
                        let pos = i + 1 - self.k;
                        if let Some(&prev_pos) = last_pos_for_code.get(&code) {
                            if pos - prev_pos < min_same_kmer_dist {
                                last_pos_for_code.insert(code, pos);
                                continue;
                            }
                        }
                        last_pos_for_code.insert(code, pos);
                        hits.push((pos, phi, weight));
                    }
                }
            } else {
                code = 0;
                valid_len = 0;
            }
        }

        for &(pos, phi, weight) in &hits {
            if pos >= phi {
                let v = pos - phi;
                if v < n {
                    raw_votes[v] += weight;
                }
            } else {
                let v1 = pos + self.period - phi;
                if v1 < n {
                    raw_votes[v1] += weight * 0.75;
                }
            }
        }

        let mut smoothed = vec![0.0f64; n];
        for i in 0..n {
            if raw_votes[i] > 0.0 {
                let w0 = raw_votes[i];
                smoothed[i] += w0;
                if i >= 1 { smoothed[i - 1] += 0.60 * w0; }
                if i + 1 < n { smoothed[i + 1] += 0.60 * w0; }
                if i >= 2 { smoothed[i - 2] += 0.25 * w0; }
                if i + 2 < n { smoothed[i + 2] += 0.25 * w0; }
            }
        }

        let mut cand_peaks: Vec<(usize, f64)> = Vec::new();
        for i in 0..n {
            let s = smoothed[i];
            if s < 0.35 {
                continue;
            }
            let left = if i > 0 { smoothed[i - 1] } else { 0.0 };
            let right = if i + 1 < n { smoothed[i + 1] } else { 0.0 };
            if s >= left && s >= right && (s > left || s > right || raw_votes[i] > 0.0) {
                cand_peaks.push((i, s));
            }
        }

        if cand_peaks.is_empty() {
            return Vec::new();
        }

        let mut peaks: Vec<(usize, f64)> = Vec::with_capacity(cand_peaks.len());
        for (pos, s) in cand_peaks {
            if let Some(last) = peaks.last_mut() {
                if pos - last.0 <= 3 {
                    if s > last.1 {
                        *last = (pos, s);
                    }
                    continue;
                }
            }
            peaks.push((pos, s));
        }

        let m = peaks.len();
        let min_step = ((self.period as f64) * 0.72).round() as usize;
        let max_step = ((self.period as f64) * 3.40).round() as usize;
        let p_f = self.period as f64;

        let mut dp = vec![0.0f64; m];
        let mut prev_idx: Vec<Option<usize>> = vec![None; m];

        for j in 0..m {
            let (pj, sj) = peaks[j];
            let start_pen = ((pj / self.period) as f64) * 0.15;
            dp[j] = sj - start_pen;

            for i in (0..j).rev() {
                let (pi, _) = peaks[i];
                let step = pj - pi;
                if step < min_step {
                    continue;
                }
                if step > max_step {
                    break;
                }
                let mult = ((step as f64) / p_f).round().clamp(1.0, 3.0);
                let rel_err = ((step as f64) - mult * p_f).abs() / p_f;
                let mult_pen = (mult - 1.0) * 0.20;
                let step_bonus = (0.45 - 1.8 * rel_err - mult_pen).max(-0.80);

                let cand_score = dp[i] + sj + step_bonus;
                if cand_score > dp[j] {
                    dp[j] = cand_score;
                    prev_idx[j] = Some(i);
                }
            }
        }

        let mut best_end = 0usize;
        let mut best_total = f64::NEG_INFINITY;
        for j in 0..m {
            let rem_periods = ((n - peaks[j].0) / self.period) as f64;
            let total = dp[j] - rem_periods * 0.15;
            if total > best_total {
                best_total = total;
                best_end = j;
            }
        }

        let mut cuts = Vec::new();
        let mut cur = Some(best_end);
        while let Some(idx) = cur {
            cuts.push(peaks[idx].0);
            cur = prev_idx[idx];
        }
        cuts.reverse();
        cuts
    }

    /// Optimal all-position Viterbi dynamic-programming cutter on $\mathbb{Z}/P\mathbb{Z}$.
    /// Guarantees every monomer step lies within $[0.90\,P_{\text{eff}}, 1.10\,P_{\text{eff}}]$,
    /// eliminating phase-shifted pairs (e.g., 75 bp + 150 bp) and multiplet alignment artifacts.
    /// When `pin_endpoints` is true (e.g., cutting a HOR unit whose start `0` and end `n` are
    /// already at the same HOR anchor phase), the path is pinned from `0` to `n` so that the
    /// HOR copy is partitioned into an exact integer number of full monomers without edge fragments.
    pub fn cut_sequence_dp(&self, seq: &[u8], pin_endpoints: bool) -> Vec<usize> {
        let n = seq.len();
        if n < self.period || n < self.k {
            return Vec::new();
        }

        let kmask = (1u32 << (2 * self.k)) - 1;
        let min_same_kmer_dist = (self.period * 2) / 3;
        let mut last_pos_for_code: HashMap<u32, usize> = HashMap::new();
        let mut raw_votes = vec![0.0f64; n];

        let mut code = 0u32;
        let mut valid_len = 0usize;
        for (i, &b) in seq.iter().enumerate() {
            if let Some(v) = encode_base(b) {
                code = ((code << 2) | v) & kmask;
                valid_len += 1;
                if valid_len >= self.k {
                    if let Some((phi, weight)) = self.phased_table[code as usize] {
                        let pos = i + 1 - self.k;
                        if let Some(&prev_pos) = last_pos_for_code.get(&code) {
                            if pos - prev_pos < min_same_kmer_dist {
                                last_pos_for_code.insert(code, pos);
                                continue;
                            }
                        }
                        last_pos_for_code.insert(code, pos);
                        if pos >= phi {
                            let v_idx = pos - phi;
                            if v_idx < n {
                                raw_votes[v_idx] += weight;
                            }
                        } else {
                            let v1 = pos + self.period - phi;
                            if v1 < n {
                                raw_votes[v1] += weight * 0.75;
                            }
                        }
                    }
                }
            } else {
                code = 0;
                valid_len = 0;
            }
        }

        let mut smoothed = vec![0.0f64; n];
        for i in 0..n {
            if raw_votes[i] > 0.0 {
                let w0 = raw_votes[i];
                smoothed[i] += w0;
                if i >= 1 { smoothed[i - 1] += 0.60 * w0; }
                if i + 1 < n { smoothed[i + 1] += 0.60 * w0; }
                if i >= 2 { smoothed[i - 2] += 0.25 * w0; }
                if i + 2 < n { smoothed[i + 2] += 0.25 * w0; }
            }
        }

        // Normalize smoothed votes by the mean local window maximum so spring penalty
        // is scale-invariant across different constellation sizes and k-mer densities.
        let win = self.period.max(10);
        let mut win_max_sum = 0.0f64;
        let mut win_count = 0usize;
        let mut w_start = 0usize;
        while w_start + win <= n {
            let w_max = smoothed[w_start..w_start + win]
                .iter()
                .copied()
                .fold(0.0f64, f64::max);
            if w_max > 0.0 {
                win_max_sum += w_max;
                win_count += 1;
            }
            w_start += win;
        }
        let peak_scale = if win_count > 0 {
            (win_max_sum / win_count as f64).max(0.5)
        } else {
            1.0
        };

        let m_est = ((n as f64) / (self.period as f64)).round().max(2.0) as usize;
        let p_ratio = (n as f64) / (m_est as f64);
        let use_pinned = pin_endpoints
            && (p_ratio - (self.period as f64)).abs() <= 0.08 * (self.period as f64);
        let p_eff = if use_pinned { p_ratio } else { self.period as f64 };

        let d_min = ((p_eff * 0.94).round() as usize).max(5);
        let d_max = ((p_eff * 1.06).round() as usize).max(d_min + 2);

        if use_pinned {
            let mut dp = vec![f64::NEG_INFINITY; n + 1];
            let mut prev = vec![0usize; n + 1];
            dp[0] = 0.0;

            for i in d_min..=n {
                let v = if i < n {
                    (smoothed[i] / peak_scale).min(1.35) - 0.50
                } else {
                    0.0
                };
                let hi_d = d_max.min(i);
                for d in d_min..=hi_d {
                    let p_idx = i - d;
                    let prev_score = dp[p_idx];
                    if prev_score == f64::NEG_INFINITY {
                        continue;
                    }
                    let dev = (d as f64 - p_eff) / p_eff;
                    let cand = prev_score + v - 36.0 * dev * dev;
                    if cand > dp[i] {
                        dp[i] = cand;
                        prev[i] = p_idx;
                    }
                }
            }

            if dp[n] > f64::NEG_INFINITY {
                let mut cuts = Vec::with_capacity(m_est);
                let mut cur = n;
                while cur > 0 {
                    let p_idx = prev[cur];
                    cuts.push(p_idx);
                    cur = p_idx;
                }
                cuts.reverse();
                return cuts;
            }
        }

        // Unpinned endpoints (for whole arrays with flanks or when n is not a clean multiple)
        let mut dp = vec![f64::NEG_INFINITY; n];
        let mut prev: Vec<Option<usize>> = vec![None; n];
        let start_lim = d_max.min(n);
        for i in 0..start_lim {
            let v0 = (smoothed[i] / peak_scale).min(1.35);
            dp[i] = v0 - (i as f64 / p_eff) * 0.35;
        }
        for i in d_min..n {
            let v = (smoothed[i] / peak_scale).min(1.35) - 0.35;
            let hi_d = d_max.min(i);
            for d in d_min..=hi_d {
                let p_idx = i - d;
                let prev_score = dp[p_idx];
                if prev_score == f64::NEG_INFINITY {
                    continue;
                }
                let dev = (d as f64 - p_eff) / p_eff;
                let cand = prev_score + v - 36.0 * dev * dev;
                if cand > dp[i] {
                    dp[i] = cand;
                    prev[i] = Some(p_idx);
                }
            }
        }

        let end_start = n.saturating_sub(d_max);
        let mut best_end = end_start;
        let mut best_total = f64::NEG_INFINITY;
        for i in end_start..n {
            let total = dp[i] - ((n - i) as f64 / p_eff) * 0.35;
            if total > best_total {
                best_total = total;
                best_end = i;
            }
        }

        if best_total == f64::NEG_INFINITY {
            return self.project_cuts(seq);
        }

        let mut cuts = Vec::with_capacity(m_est);
        let mut cur = Some(best_end);
        while let Some(idx) = cur {
            cuts.push(idx);
            cur = prev[idx];
        }
        cuts.reverse();
        cuts
    }

    /// Evaluate whether the cuts projected on `seqs` pass the structural verification gate.
    /// Returns `(passed, coverage, frac_1x)`.
    pub fn verify_cuts_on_ensemble(&self, seqs: &[&[u8]]) -> (bool, f64, f64) {
        let total_len: usize = seqs.iter().map(|s| s.len()).sum();
        let n_exp = (total_len as f64 / self.period as f64).max(1.0);

        let mut total_cuts = 0usize;
        let mut n_1x = 0usize;
        let mut n_1to3x = 0usize;
        let mut n_intervals = 0usize;

        let lo_1x = ((self.period as f64) * 0.74).round() as usize;
        let hi_1x = ((self.period as f64) * 1.30).round() as usize;
        let hi_3x = ((self.period as f64) * 3.35).round() as usize;

        for &seq in seqs {
            let cuts = self.project_cuts(seq);
            total_cuts += cuts.len();
            for w in cuts.windows(2) {
                let d = w[1] - w[0];
                n_intervals += 1;
                if d >= lo_1x && d <= hi_1x {
                    n_1x += 1;
                }
                if d >= lo_1x && d <= hi_3x {
                    n_1to3x += 1;
                }
            }
        }

        let coverage = total_cuts as f64 / n_exp;
        let frac_1x = if n_intervals > 0 {
            n_1x as f64 / n_intervals as f64
        } else {
            0.0
        };
        let frac_1to3x = if n_intervals > 0 {
            n_1to3x as f64 / n_intervals as f64
        } else {
            0.0
        };

        let passed = total_cuts >= 4
            && self.num_phased >= 3
            && coverage >= 0.35
            && frac_1x >= 0.40
            && frac_1to3x >= 0.75;

        (passed, coverage, frac_1x)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_diverged_monomer_period_and_constellation() {
        let base_template = b"TTTTCAAACACTTTTATCTTAGAGAAAAAAGATGTTAGTGAAAAGGTTTTTCCCCAAGTTTTTTTTCCAAAGGGTTGGTTATTGGCCCTTTCCCGCTAGGCCCCCCCTGGGT";
        assert_eq!(base_template.len(), 112);

        let mut array = Vec::new();
        let bases = [b'A', b'C', b'G', b'T'];
        let mut rng: u64 = 123456789;
        let mut next_u64 = || -> u64 {
            rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
            rng
        };

        for copy in 0..25 {
            let mut mono = base_template.to_vec();
            let spare1 = (copy % 3) * 32;
            let spare2 = ((copy + 1) % 3) * 32;
            for pos in 0..mono.len() {
                if (pos >= spare1 && pos < spare1 + 10) || (pos >= spare2 && pos < spare2 + 10) {
                    continue;
                }
                if (next_u64() % 100) < 38 {
                    mono[pos] = bases[(next_u64() % 4) as usize];
                }
            }
            if copy % 4 == 0 {
                mono.push(b'A');
            } else if copy % 5 == 0 {
                mono.pop();
            }
            array.extend_from_slice(&mono);
        }

        let cand = find_monomer_period_kmer(&array, 20, 500).expect("Should detect monomer period");
        assert!(
            (cand.period as isize - 112).abs() <= 3,
            "Expected period ~112, got {}",
            cand.period
        );

        let constel = PhasedConstellation::build(&[&array], cand.period, false)
            .expect("Should build phased constellation");
        let (passed, cov, frac_1x) = constel.verify_cuts_on_ensemble(&[&array]);
        assert!(passed, "Constellation verification failed: cov={}, frac_1x={}", cov, frac_1x);
        assert!(cov >= 0.80, "Expected >=80% cut coverage, got {}", cov);
    }
}
