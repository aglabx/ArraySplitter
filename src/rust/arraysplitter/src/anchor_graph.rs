//! Anchor Graph for monomer decomposition.
//!
//! Key idea: Build a graph where:
//! - Nodes = conserved parts (anchors) from FS-tree hints
//! - Edges = transitions between anchors with support (copy number)
//!
//! This solves two problems:
//! 1. Multiple conserved parts per monomer (overcutting)
//! 2. Mutations in anchor (missing monomers) - handled by using shorter anchors
//!    where longer ones are absent
//!
//! The graph structure reveals the monomer architecture.

use std::collections::{HashMap, HashSet};
use triple_accel::levenshtein_exp;

/// Single occurrence of an anchor in sequence.
#[derive(Debug, Clone)]
pub struct AnchorHit {
    pub anchor: String,
    pub position: usize,
    pub length: usize,
}

impl AnchorHit {
    pub fn new(anchor: String, position: usize) -> Self {
        let length = anchor.len();
        Self { anchor, position, length }
    }
}

impl PartialEq for AnchorHit {
    fn eq(&self, other: &Self) -> bool {
        self.position == other.position
    }
}

impl Eq for AnchorHit {}

impl PartialOrd for AnchorHit {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.position.cmp(&other.position))
    }
}

impl Ord for AnchorHit {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.position.cmp(&other.position)
    }
}

/// Result of anchor graph decomposition
#[derive(Debug, Clone)]
pub struct GraphDecomposition {
    pub monomers: Vec<String>,
    pub cut_sequence: String,
    pub period: usize,
    pub cv: f64,
    pub alternative_anchors: Vec<String>,
}

// =============================================================================
// Distance Clustering Functions (matching Python exactly)
// =============================================================================

/// Simple k-means clustering without external dependencies.
fn simple_kmeans(data: &[f64], k: usize, max_iter: usize) -> (Vec<Vec<f64>>, Vec<f64>) {
    if data.is_empty() || k == 0 {
        return (vec![Vec::new()], vec![0.0]);
    }

    if k == 1 {
        let mean = data.iter().sum::<f64>() / data.len() as f64;
        return (vec![data.to_vec()], vec![mean]);
    }

    let mut data_sorted: Vec<f64> = data.to_vec();
    data_sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = data_sorted.len();

    // Initialize centers evenly across the data range
    let mut centers: Vec<f64> = (0..k)
        .map(|i| data_sorted[(i * n / k).min(n - 1)])
        .collect();

    for _ in 0..max_iter {
        // Assign points to nearest center
        let mut clusters: Vec<Vec<f64>> = vec![Vec::new(); k];
        for &x in data {
            let nearest = centers
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| {
                    (x - *a).abs().partial_cmp(&(x - *b).abs()).unwrap()
                })
                .map(|(i, _)| i)
                .unwrap_or(0);
            clusters[nearest].push(x);
        }

        // Update centers
        let mut new_centers: Vec<f64> = Vec::with_capacity(k);
        for (i, c) in clusters.iter().enumerate() {
            if c.is_empty() {
                new_centers.push(centers[i]);
            } else {
                new_centers.push(c.iter().sum::<f64>() / c.len() as f64);
            }
        }

        if new_centers == centers {
            break;
        }
        centers = new_centers;
    }

    // Final assignment
    let mut result: Vec<Vec<f64>> = vec![Vec::new(); k];
    for &val in data {
        let closest = centers
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                (val - *a).abs().partial_cmp(&(val - *b).abs()).unwrap()
            })
            .map(|(i, _)| i)
            .unwrap_or(0);
        result[closest].push(val);
    }

    (result, centers)
}

/// Find optimal number of clusters (1 to max_k).
///
/// Criteria:
/// - All clusters must be significant (>15% of data)
/// - Cluster centers must be well separated (gap > 20% of max center)
fn find_optimal_clusters(data: &[f64], max_k: usize) -> (usize, Vec<Vec<f64>>, Vec<f64>) {
    if data.len() < 4 {
        let median = {
            let mut sorted = data.to_vec();
            sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
            sorted[sorted.len() / 2]
        };
        return (1, vec![data.to_vec()], vec![median]);
    }

    let mut best_k = 1;
    let mut best_score = 0.0_f64;
    let mut best_clusters = vec![data.to_vec()];
    let mut best_centers = vec![{
        let mut sorted = data.to_vec();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        sorted[sorted.len() / 2]
    }];

    for k in 2..=(max_k.min(data.len() / 2)) {
        let (clusters, centers) = simple_kmeans(data, k, 100);

        // Filter out empty clusters
        let non_empty: Vec<(Vec<f64>, f64)> = clusters
            .into_iter()
            .zip(centers.into_iter())
            .filter(|(c, _)| !c.is_empty())
            .collect();

        if non_empty.len() < k {
            continue;
        }

        let clusters: Vec<Vec<f64>> = non_empty.iter().map(|(c, _)| c.clone()).collect();
        let centers: Vec<f64> = non_empty.iter().map(|(_, ctr)| *ctr).collect();

        // Check that all clusters are significant (>15% of data)
        let min_size = clusters.iter().map(|c| c.len()).min().unwrap_or(0);
        if (min_size as f64) < (data.len() as f64 * 0.15) {
            continue;
        }

        // Check that centers are well separated
        let mut centers_sorted = centers.clone();
        centers_sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        if centers_sorted.len() < 2 {
            continue;
        }

        let min_gap = centers_sorted
            .windows(2)
            .map(|w| w[1] - w[0])
            .fold(f64::MAX, |a, b| a.min(b));
        let max_center = centers_sorted.iter().cloned().fold(f64::MIN, f64::max);

        if min_gap < 0.2 * max_center {
            continue;
        }

        // Score: separation quality
        let score = min_gap / max_center;
        if score > best_score {
            best_score = score;
            best_k = k;
            best_clusters = clusters;
            best_centers = centers;
        }
    }

    (best_k, best_clusters, best_centers)
}

/// Analyze distance distribution to detect multi-region patterns.
///
/// If anchor appears N times per monomer, distances will have N clusters.
/// True period = sum of cluster centers.
fn analyze_distance_distribution(distances: &[f64], _verbose: bool) -> (usize, Vec<f64>, f64) {
    if distances.len() < 4 {
        let med = if distances.is_empty() {
            0.0
        } else {
            let mut sorted = distances.to_vec();
            sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
            sorted[sorted.len() / 2]
        };
        return (1, vec![med], med);
    }

    let (n_clusters, _, centers) = find_optimal_clusters(distances, 4);

    // Sort centers by value for consistent ordering
    let mut centers_sorted = centers;
    centers_sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let true_period: f64 = centers_sorted.iter().sum();

    (n_clusters, centers_sorted, true_period)
}

/// Classify each position by the cluster of its distance to the next position.
fn classify_positions_by_distance(positions: &[usize], centers: &[f64]) -> Vec<i32> {
    if positions.len() < 2 {
        return vec![-1; positions.len()];
    }

    let mut labels: Vec<i32> = Vec::with_capacity(positions.len());

    for i in 0..positions.len() - 1 {
        let dist = (positions[i + 1] - positions[i]) as f64;
        // Find nearest cluster center
        let nearest = centers
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                (dist - *a).abs().partial_cmp(&(dist - *b).abs()).unwrap()
            })
            .map(|(i, _)| i as i32)
            .unwrap_or(0);
        labels.push(nearest);
    }

    // Last position has no "next", mark with -1
    labels.push(-1);

    labels
}

/// Select cut positions based on cluster analysis.
///
/// Strategy: Cut at positions with the same label as the first position.
/// This preserves the monomer structure starting from the left flank.
fn select_cut_positions_by_cluster(positions: &[usize], labels: &[i32], _verbose: bool) -> Vec<usize> {
    if positions.is_empty() || labels.is_empty() {
        return positions.to_vec();
    }

    // Use the label of the first position as the "start" label
    let start_label = labels[0];

    if start_label < 0 {
        return positions.to_vec();
    }

    // Select all positions with the start label
    let mut cut_positions: Vec<usize> = positions
        .iter()
        .zip(labels.iter())
        .filter(|(_, &label)| label == start_label)
        .map(|(&pos, _)| pos)
        .collect();

    // Also include the last position if it wasn't included
    if let Some(&last) = positions.last() {
        if !cut_positions.contains(&last) {
            cut_positions.push(last);
        }
    }

    cut_positions.sort();
    cut_positions
}

/// Find all anchor occurrences, using longer anchors first.
///
/// Logic:
/// 1. Sort anchors by length (longest first)
/// 2. Find all positions of longest anchor
/// 3. For shorter anchors - only find in regions NOT covered by longer ones
/// 4. This handles mutations: if long anchor has mutation, short variant will be found
fn find_anchor_hits_hierarchical(sequence: &str, anchors: &[String], verbose: bool) -> Vec<AnchorHit> {
    if anchors.is_empty() {
        return Vec::new();
    }

    // Sort by length descending
    let mut sorted_anchors = anchors.to_vec();
    sorted_anchors.sort_by(|a, b| b.len().cmp(&a.len()));

    // Track covered positions
    let mut covered: HashSet<usize> = HashSet::new();
    let mut all_hits: Vec<AnchorHit> = Vec::new();

    for anchor in &sorted_anchors {
        let anchor_len = anchor.len();

        // Find all occurrences
        let mut start = 0;
        let mut hits_for_anchor = 0;

        while let Some(pos) = sequence[start..].find(anchor.as_str()) {
            let absolute_pos = start + pos;

            // Check if this position is already covered by a longer anchor
            let is_covered = (absolute_pos..absolute_pos + anchor_len).any(|p| covered.contains(&p));

            if !is_covered {
                all_hits.push(AnchorHit::new(anchor.clone(), absolute_pos));
                // Mark positions as covered
                for p in absolute_pos..absolute_pos + anchor_len {
                    covered.insert(p);
                }
                hits_for_anchor += 1;
            }

            start = absolute_pos + 1;
            if start >= sequence.len() {
                break;
            }
        }

        if verbose && hits_for_anchor > 0 {
            eprintln!("  '{}' (len={}): {} hits", anchor, anchor_len, hits_for_anchor);
        }
    }

    // Sort by position
    all_hits.sort();

    all_hits
}

/// Build graph of transitions between consecutive anchor hits.
fn build_transition_graph(hits: &[AnchorHit]) -> (HashMap<(String, String), usize>, HashMap<(String, String), Vec<usize>>) {
    let mut edges: HashMap<(String, String), usize> = HashMap::new();
    let mut distances: HashMap<(String, String), Vec<usize>> = HashMap::new();

    for i in 0..hits.len().saturating_sub(1) {
        let curr = &hits[i];
        let next_hit = &hits[i + 1];

        let edge = (curr.anchor.clone(), next_hit.anchor.clone());
        let distance = next_hit.position - curr.position;

        *edges.entry(edge.clone()).or_insert(0) += 1;
        distances.entry(edge).or_insert_with(Vec::new).push(distance);
    }

    (edges, distances)
}

/// Calculate similarity ratio between two sequences using edit distance.
/// Returns a value between 0.0 (completely different) and 1.0 (identical).
fn sequence_similarity(s1: &str, s2: &str) -> f64 {
    if s1.is_empty() && s2.is_empty() {
        return 1.0;
    }
    let max_len = s1.len().max(s2.len());
    if max_len == 0 {
        return 1.0;
    }
    let edit_dist = levenshtein_exp(s1.as_bytes(), s2.as_bytes()) as f64;
    1.0 - (edit_dist / max_len as f64)
}

/// Calculate raw edit distance between two sequences.
fn raw_edit_distance(s1: &str, s2: &str) -> usize {
    levenshtein_exp(s1.as_bytes(), s2.as_bytes()) as usize
}

/// Find the best anchor for decomposition using distance clustering.
///
/// Strategy:
/// 1. For each anchor with self-loop, analyze distance distribution
/// 2. Cluster distances to detect multi-region patterns (A-B, A-B-C, etc.)
/// 3. True period = sum of cluster centers
/// 4. Check that resulting monomers are SIMILAR (low edit distance)
/// 5. Select anchor with lowest CV after cluster-based correction
fn find_monomer_cycle(
    sequence: &str,
    edges: &HashMap<(String, String), usize>,
    _distances: &HashMap<(String, String), Vec<usize>>,
    hits: &[AnchorHit],
    verbose: bool,
) -> (Vec<String>, f64, f64, usize, Vec<f64>, Vec<String>) {
    if edges.is_empty() {
        return (Vec::new(), 0.0, f64::INFINITY, 1, Vec::new(), Vec::new());
    }

    // Two-pass approach:
    // Pass 1: Quick filtering by CV and cut count (no expensive ED calculation)
    // Pass 2: Full similarity calculation only for top candidates

    const MIN_CUT_COUNT: usize = 10;
    const MAX_CV_FOR_ED: f64 = 2.0;  // Calculate ED for all candidates (was 1.0)
    const TOP_CANDIDATES_FOR_ED: usize = 10;  // Calculate ED for top N candidates by CV
    const MIN_SIMILARITY: f64 = 0.3;  // Alpha satellite monomers are highly variable
    const MAX_MONOMER_FOR_ED: usize = 2000;  // Skip ED for very long monomers
    const MIN_PERIOD: f64 = 30.0;  // Skip sub-word noise (real sat monomers are >30bp)

    // Tuple for pass 1: (anchor, cv, period, cut_count, n_clusters, centers, cut_positions)
    let mut pass1_candidates: Vec<(String, f64, f64, usize, usize, Vec<f64>, Vec<usize>)> = Vec::new();

    // Get unique anchors from graph with self-loops
    let mut anchors_to_try: HashSet<String> = HashSet::new();
    for ((from_a, to_a), &count) in edges.iter() {
        if count >= 3 {
            anchors_to_try.insert(from_a.clone());
            if from_a != to_a {
                anchors_to_try.insert(to_a.clone());
            }
        }
    }

    // Pass 1: Calculate CV and filter by cut count (fast)
    for anchor in &anchors_to_try {
        let mut positions: Vec<usize> = hits
            .iter()
            .filter(|h| &h.anchor == anchor)
            .map(|h| h.position)
            .collect();

        if positions.len() < 3 {
            continue;
        }
        positions.sort();

        let dist_list: Vec<f64> = positions
            .windows(2)
            .map(|w| (w[1] - w[0]) as f64)
            .collect();

        if dist_list.len() < 3 {
            continue;
        }

        let (n_clusters, centers, true_period) = analyze_distance_distribution(&dist_list, false);
        let labels = classify_positions_by_distance(&positions, &centers);
        let cut_positions = select_cut_positions_by_cluster(&positions, &labels, false);

        if cut_positions.len() < MIN_CUT_COUNT {
            if verbose {
                eprintln!("  Skipping '{}': only {} cuts (minimum {})",
                    &anchor[..anchor.len().min(20)], cut_positions.len(), MIN_CUT_COUNT);
            }
            continue;
        }

        // Skip candidates with very short period (sub-word noise, not real monomers)
        let monomer_len_est = true_period / (n_clusters as f64).max(1.0);
        if monomer_len_est < MIN_PERIOD {
            if verbose {
                eprintln!("  Skipping '{}': period too short ({:.0}bp < {:.0}bp)",
                    &anchor[..anchor.len().min(20)], monomer_len_est, MIN_PERIOD);
            }
            continue;
        }

        let lengths: Vec<f64> = cut_positions
            .windows(2)
            .map(|w| (w[1] - w[0]) as f64)
            .collect();

        if lengths.len() < 2 {
            continue;
        }

        let mean_len: f64 = lengths.iter().sum::<f64>() / lengths.len() as f64;
        let variance: f64 = lengths.iter().map(|&x| (x - mean_len).powi(2)).sum::<f64>() / lengths.len() as f64;
        let cv = if mean_len > 0.0 { variance.sqrt() / mean_len } else { f64::INFINITY };

        pass1_candidates.push((anchor.clone(), cv, true_period, cut_positions.len(), n_clusters, centers, cut_positions));
    }

    // Sort pass1 by CV to select top candidates for ED calculation
    pass1_candidates.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    if verbose && !pass1_candidates.is_empty() {
        eprintln!("  Pass 1: {} candidates after CV/cuts filter", pass1_candidates.len());
        eprintln!("  All candidates (sorted by CV):");
        eprintln!("  {:20} {:>8} {:>8} {:>6} {:>8}", "Anchor", "Period", "CV", "Cuts", "Clusters");
        for (anchor, cv, period, cuts, n_clusters, _, _) in &pass1_candidates {
            let anchor_display = if anchor.len() > 18 {
                format!("{}...", &anchor[..15])
            } else {
                anchor.clone()
            };
            eprintln!("  {:20} {:>8.0} {:>8.3} {:>6} {:>8}",
                anchor_display, period, cv, cuts, n_clusters);
        }
        eprintln!();
    }

    // Pass 2: Calculate full similarity only for top candidates
    // Tuple: (cycle, cv, period, cut_count, n_clusters, centers, mean_ed, std_ed, combined_score)
    //
    // Score formula modes (SCORE_MODE):
    //   0: mean_ed                    (baseline)
    //   1: mean_ed / sqrt(cuts)       (statistical confidence weighting)
    //   2: mean_ed / log2(cuts)       (gentler confidence weighting)
    //   3: mean_ed / cuts             (linear cut weighting)
    //   4: mean_ed * (1 + 10/sqrt(cuts))  (penalty for few cuts)
    const SCORE_MODE: usize = 3;  // mean_ed / cuts (best overall)
    let mut candidates: Vec<(Vec<String>, f64, f64, usize, usize, Vec<f64>, f64, f64, f64)> = Vec::new();

    let candidates_for_ed = pass1_candidates.iter()
        .filter(|(_, cv, _, _, _, _, _)| *cv <= MAX_CV_FOR_ED)
        .take(TOP_CANDIDATES_FOR_ED)
        .collect::<Vec<_>>();

    // Also include candidates with CV > MAX_CV_FOR_ED but skip ED for them
    let remaining_candidates = pass1_candidates.iter()
        .filter(|(_, cv, _, _, _, _, _)| *cv > MAX_CV_FOR_ED)
        .collect::<Vec<_>>();

    // Calculate full ED for top candidates
    // Now we use mean_ed + lambda*std_ed (lower is better) instead of similarity
    for (anchor, cv, true_period, cut_count, n_clusters, centers, cut_positions) in &candidates_for_ed {
        let mut ed_values: Vec<f64> = Vec::new();
        let mut similarity_sum = 0.0;

        // Calculate ED for ALL consecutive pairs (not just 10)
        for i in 0..(cut_positions.len().saturating_sub(2)) {
            let start1 = cut_positions[i];
            let end1 = cut_positions[i + 1];
            let monomer1 = &sequence[start1..end1];

            let start2 = cut_positions[i + 1];
            let end2 = cut_positions[i + 2];
            let monomer2 = &sequence[start2..end2];

            // Skip very long monomers (ED is O(n*m))
            if monomer1.len() <= MAX_MONOMER_FOR_ED && monomer2.len() <= MAX_MONOMER_FOR_ED {
                let ed = raw_edit_distance(monomer1, monomer2);
                let sim = sequence_similarity(monomer1, monomer2);
                ed_values.push(ed as f64);
                similarity_sum += sim;
            }
        }

        let ed_count = ed_values.len();

        let (mean_ed, std_ed) = if ed_count > 0 {
            let mean = ed_values.iter().sum::<f64>() / ed_count as f64;
            let variance = if ed_count > 1 {
                ed_values.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / ed_count as f64
            } else {
                0.0
            };
            (mean, variance.sqrt())
        } else {
            (f64::INFINITY, f64::INFINITY)
        };

        // Combined score: lower is better
        // Uses cut count as confidence weight (more cuts = more reliable mean_ed)
        // Normalizes ED by actual monomer length (period/n_clusters) to compare across lengths
        let cuts_f64 = *cut_count as f64;
        let monomer_len = *true_period / (*n_clusters as f64).max(1.0);
        let combined_score = if mean_ed.is_finite() && cuts_f64 > 0.0 && monomer_len > 0.0 {
            match SCORE_MODE {
                0 => mean_ed,                                    // baseline
                1 => mean_ed / cuts_f64.sqrt(),                  // sqrt weighting
                2 => mean_ed / cuts_f64.log2().max(1.0),         // log2 weighting
                3 => mean_ed / cuts_f64,                         // linear weighting (old best)
                4 => mean_ed * (1.0 + 10.0 / cuts_f64.sqrt()),  // penalty for few cuts
                5 => (mean_ed / monomer_len) / cuts_f64,         // ed_per_bp / cuts
                _ => mean_ed,
            }
        } else {
            f64::INFINITY
        };

        let mean_similarity = if ed_count > 0 {
            similarity_sum / ed_count as f64
        } else {
            0.0
        };

        // Filter by minimum similarity
        if mean_similarity < MIN_SIMILARITY {
            if verbose {
                eprintln!("  Skipping '{}': monomers too dissimilar (sim={:.2} < {:.2})",
                    &anchor[..anchor.len().min(20)], mean_similarity, MIN_SIMILARITY);
            }
            continue;
        }

        // Filter by ed_per_bp: if ED > monomer length, these are not real repeats
        let ed_per_bp = mean_ed / monomer_len;
        if ed_per_bp > 1.0 {
            if verbose {
                eprintln!("  Skipping '{}': ed_per_bp={:.2} > 1.0 (ED exceeds monomer length, not a repeat)",
                    &anchor[..anchor.len().min(20)], ed_per_bp);
            }
            continue;
        }

        if verbose {
            let cluster_info = if *n_clusters > 1 {
                format!(" ({} clusters)", n_clusters)
            } else {
                String::new()
            };
            eprintln!("  Candidate: {}... period={:.0}, CV={:.3}, meanED={:.1}, stdED={:.1}, score={:.1}, sim={:.2}, cuts={}{} [full ED]",
                &anchor[..anchor.len().min(20)], true_period, cv, mean_ed, std_ed, combined_score, mean_similarity, cut_count, cluster_info);
        }

        // Store mean_ed, std_ed, combined_score (lower is better)
        candidates.push((vec![anchor.clone()], *cv, *true_period, *cut_count, *n_clusters, centers.clone(), mean_ed, std_ed, combined_score));
    }

    // Add remaining candidates without full ED (combined_score=INFINITY so they lose to candidates with ED)
    for (anchor, cv, true_period, cut_count, n_clusters, centers, _cut_positions) in &remaining_candidates {
        if verbose {
            let cluster_info = if *n_clusters > 1 {
                format!(" ({} clusters)", n_clusters)
            } else {
                String::new()
            };
            eprintln!("  Candidate: {}... period={:.0}, CV={:.3}, meanED=N/A, stdED=N/A, score=N/A, cuts={}{} [no ED - high CV]",
                &anchor[..anchor.len().min(20)], true_period, cv, cut_count, cluster_info);
        }
        candidates.push((vec![anchor.clone()], *cv, *true_period, *cut_count, *n_clusters, centers.clone(), f64::INFINITY, f64::INFINITY, f64::INFINITY));
    }

    // Select best candidate
    if candidates.is_empty() {
        // Fallback: most frequent edge
        if !edges.is_empty() {
            let start_edge = edges
                .iter()
                .max_by_key(|(_, &count)| count)
                .map(|((from, _), _)| from.clone())
                .unwrap();
            if verbose {
                eprintln!("  Fallback to most frequent: {}...", &start_edge[..start_edge.len().min(20)]);
            }
            return (vec![start_edge], 0.0, f64::INFINITY, 1, Vec::new(), Vec::new());
        }
        return (Vec::new(), 0.0, f64::INFINITY, 1, Vec::new(), Vec::new());
    }

    // Show all candidates with ED before sorting
    if verbose && !candidates.is_empty() {
        let mode_name = match SCORE_MODE {
            0 => "mean_ed",
            1 => "mean_ed/sqrt(cuts)",
            2 => "mean_ed/log2(cuts)",
            3 => "mean_ed/cuts",
            4 => "mean_ed*(1+10/sqrt(cuts))",
            5 => "(mean_ed/period)/cuts",
            _ => "unknown",
        };
        eprintln!("  Final candidates (before sorting, score={}):", mode_name);
        eprintln!("  {:20} {:>8} {:>8} {:>8} {:>8} {:>8} {:>6}", "Anchor", "Period", "CV", "MeanED", "StdED", "Score", "Cuts");
        for (anchors, cv, period, cuts, _, _, mean_ed, std_ed, score) in &candidates {
            let anchor = &anchors[0];
            let anchor_display = if anchor.len() > 18 {
                format!("{}...", &anchor[..15])
            } else {
                anchor.clone()
            };
            let ed_str = if mean_ed.is_finite() { format!("{:.1}", mean_ed) } else { "N/A".to_string() };
            let std_str = if std_ed.is_finite() { format!("{:.1}", std_ed) } else { "N/A".to_string() };
            let score_str = if score.is_finite() { format!("{:.1}", score) } else { "N/A".to_string() };
            eprintln!("  {:20} {:>8.0} {:>8.3} {:>8} {:>8} {:>8} {:>6}",
                anchor_display, period, cv, ed_str, std_str, score_str, cuts);
        }
        eprintln!();
    }

    // Sort by: combined_score (ascending), then CV (ascending), then cuts (descending), then clusters (ascending)
    // combined_score includes cut count weighting (more cuts = more confidence in mean_ed)
    // Rationale: FASTAN-like approach - anchors with more cuts (higher coverage)
    // are preferred, analogous to k-mer diagonal density in FASTAN
    candidates.sort_by(|a, b| {
        // First compare combined_score (lower is better)
        let score_cmp = a.8.partial_cmp(&b.8).unwrap_or(std::cmp::Ordering::Equal);
        if score_cmp != std::cmp::Ordering::Equal {
            return score_cmp;
        }
        // If combined_score is similar, prefer lower CV
        let cv_cmp = a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal);
        if cv_cmp != std::cmp::Ordering::Equal {
            return cv_cmp;
        }
        // Then prefer more cuts (more reliable)
        let cuts_cmp = b.3.cmp(&a.3);
        if cuts_cmp != std::cmp::Ordering::Equal {
            return cuts_cmp;
        }
        // Finally, prefer simpler patterns (fewer clusters)
        a.4.cmp(&b.4)
    });

    let mut candidates_iter = candidates.into_iter();
    let (best_cycle, best_cv, best_period, _best_count, best_n_clusters, best_centers, best_mean_ed, best_std_ed, best_score) =
        candidates_iter.next().unwrap();

    // Collect alternative anchors with similar period (±10%) and finite score
    let period_min = best_period * 0.9;
    let period_max = best_period * 1.1;
    let alt_anchors: Vec<String> = candidates_iter
        .filter(|(_, _, period, _, _, _, _, _, score)| {
            score.is_finite() && *period >= period_min && *period <= period_max
        })
        .take(5)
        .map(|(cycle, _, _, _, _, _, _, _, _)| cycle[0].clone())
        .collect();

    if verbose {
        let cluster_info = if best_n_clusters > 1 {
            format!(" ({} clusters)", best_n_clusters)
        } else {
            String::new()
        };
        let ed_info = if best_mean_ed.is_finite() && best_std_ed.is_finite() {
            format!(", meanED={:.1}, stdED={:.1}, score={:.1}", best_mean_ed, best_std_ed, best_score)
        } else {
            String::new()
        };
        eprintln!("  Selected: {}... period={:.0}, CV={:.3}{}{}",
            &best_cycle[0][..best_cycle[0].len().min(20)], best_period, best_cv, ed_info, cluster_info);
        if !alt_anchors.is_empty() {
            eprintln!("  Alternative anchors (same period): {}", alt_anchors.len());
        }
    }

    (best_cycle, best_period, best_cv, best_n_clusters, best_centers, alt_anchors)
}

/// Decompose sequence into monomers based on the anchor cycle.
///
/// Strategy:
/// 1. Find all positions of the anchor
/// 2. If cluster_centers provided (multi-region pattern), use cluster analysis
///    to select only the "start" positions
/// 3. Cut at selected positions
fn decompose_by_cycle(
    sequence: &str,
    hits: &[AnchorHit],
    cycle: &[String],
    cluster_centers: &[f64],
    verbose: bool,
) -> Vec<String> {
    if cycle.is_empty() || hits.is_empty() {
        return if sequence.is_empty() { Vec::new() } else { vec![sequence.to_string()] };
    }

    let start_anchor = &cycle[0];

    // Find all positions where anchor appears
    let mut all_positions: Vec<usize> = hits
        .iter()
        .filter(|h| &h.anchor == start_anchor)
        .map(|h| h.position)
        .collect();

    if all_positions.is_empty() {
        return vec![sequence.to_string()];
    }

    all_positions.sort();

    // Determine cut positions based on cluster analysis
    let cut_positions = if cluster_centers.len() > 1 {
        // Multi-region pattern: use cluster-based selection
        let labels = classify_positions_by_distance(&all_positions, cluster_centers);
        let cuts = select_cut_positions_by_cluster(&all_positions, &labels, verbose);

        if verbose {
            eprintln!("  Multi-region pattern: {} clusters", cluster_centers.len());
            eprintln!("  All positions: {}, cut positions: {}", all_positions.len(), cuts.len());
        }
        cuts
    } else {
        // Single region: cut at every position
        all_positions
    };

    if verbose && !cut_positions.is_empty() {
        let display: Vec<String> = cut_positions.iter().take(10).map(|p| p.to_string()).collect();
        let suffix = if cut_positions.len() > 10 { "..." } else { "" };
        eprintln!("  Cut positions: {}{}", display.join(", "), suffix);
    }

    if cut_positions.is_empty() {
        return vec![sequence.to_string()];
    }

    // Cut sequence
    let mut monomers: Vec<String> = Vec::new();

    // Left flank
    if cut_positions[0] > 0 {
        monomers.push(sequence[..cut_positions[0]].to_string());
    }

    // Monomers
    for i in 0..cut_positions.len() - 1 {
        monomers.push(sequence[cut_positions[i]..cut_positions[i + 1]].to_string());
    }

    // Last part (may be partial)
    if let Some(&last_pos) = cut_positions.last() {
        if last_pos < sequence.len() {
            monomers.push(sequence[last_pos..].to_string());
        }
    }

    monomers
}

/// Main class for anchor graph-based decomposition.
pub struct AnchorGraphDecomposer {
    sequence: String,
    anchors: Vec<String>,
    hits: Vec<AnchorHit>,
    edges: HashMap<(String, String), usize>,
    distances: HashMap<(String, String), Vec<usize>>,
    cycle: Vec<String>,
    estimated_period: f64,
    estimated_cv: f64,
    n_clusters: usize,
    cluster_centers: Vec<f64>,
    alternative_anchors: Vec<String>,
}

impl AnchorGraphDecomposer {
    pub fn new() -> Self {
        Self {
            sequence: String::new(),
            anchors: Vec::new(),
            hits: Vec::new(),
            edges: HashMap::new(),
            distances: HashMap::new(),
            cycle: Vec::new(),
            estimated_period: 0.0,
            estimated_cv: f64::INFINITY,
            n_clusters: 1,
            cluster_centers: Vec::new(),
            alternative_anchors: Vec::new(),
        }
    }

    /// Build anchor graph from compute_cuts candidates.
    ///
    /// # Arguments
    /// * `sequence` - The array sequence
    /// * `candidates` - List of candidate cut sequences (max 11bp length)
    /// * `verbose` - Print debug info
    pub fn build_from_candidates(&mut self, sequence: &str, candidates: &[String], verbose: bool) {
        self.sequence = sequence.to_string();

        // Filter by anchor length (3-11bp)
        // Min 3bp: single nucleotides like "A" produce nonsense decomposition
        // Max 11bp: longer anchors are too specific and may miss mutated copies
        const MIN_ANCHOR_LENGTH: usize = 3;
        const MAX_ANCHOR_LENGTH: usize = 11;
        self.anchors = candidates
            .iter()
            .filter(|c| c.len() >= MIN_ANCHOR_LENGTH && c.len() <= MAX_ANCHOR_LENGTH)
            .take(15)
            .cloned()
            .collect();

        if verbose {
            eprintln!("Top {} anchors (length {}-{}bp):", self.anchors.len(), MIN_ANCHOR_LENGTH, MAX_ANCHOR_LENGTH);
            for anchor in &self.anchors {
                eprintln!("  '{}'", anchor);
            }
        }

        // Find hits hierarchically (longer first)
        if verbose {
            eprintln!("\nFinding anchor hits (longer anchors first):");
        }
        self.hits = find_anchor_hits_hierarchical(sequence, &self.anchors, verbose);

        if verbose {
            eprintln!("\nTotal hits: {}", self.hits.len());
        }

        // Build transition graph
        let (edges, distances) = build_transition_graph(&self.hits);
        self.edges = edges;
        self.distances = distances;

        if verbose {
            eprintln!("\nTransition graph:");
            let mut edge_vec: Vec<_> = self.edges.iter().collect();
            edge_vec.sort_by_key(|(_, &count)| std::cmp::Reverse(count));
            for ((from, to), count) in edge_vec.iter().take(10) {
                if let Some(dists) = self.distances.get(&(from.clone(), to.clone())) {
                    let mean_dist: f64 = dists.iter().sum::<usize>() as f64 / dists.len() as f64;
                    eprintln!("  {}... -> {}...: count={}, mean_dist={:.0}",
                        &from[..from.len().min(15)], &to[..to.len().min(15)], count, mean_dist);
                }
            }
        }

        // Find monomer cycle with cluster analysis
        if verbose {
            eprintln!("\nFinding monomer cycle:");
        }
        let (cycle, period, cv, n_clusters, centers, alt_anchors) =
            find_monomer_cycle(sequence, &self.edges, &self.distances, &self.hits, verbose);

        self.cycle = cycle;
        self.estimated_period = period;
        self.estimated_cv = cv;
        self.n_clusters = n_clusters;
        self.cluster_centers = centers;
        self.alternative_anchors = alt_anchors;

        if verbose && !self.cycle.is_empty() {
            let cluster_info = if self.n_clusters > 1 {
                format!(", {} clusters", self.n_clusters)
            } else {
                String::new()
            };
            eprintln!("\nEstimated monomer length: {:.0} bp (CV={:.3}{})",
                self.estimated_period, self.estimated_cv, cluster_info);
        }
    }

    /// Decompose sequence into monomers.
    pub fn decompose(&self, verbose: bool) -> Option<GraphDecomposition> {
        if self.cycle.is_empty() {
            if verbose {
                eprintln!("No cycle found, returning whole sequence");
            }
            return None;
        }

        if verbose {
            let cluster_info = if self.n_clusters > 1 {
                format!(" ({} clusters)", self.n_clusters)
            } else {
                String::new()
            };
            let cycle_str: Vec<String> = self.cycle.iter()
                .map(|a| format!("{}...", &a[..a.len().min(10)]))
                .collect();
            eprintln!("\nDecomposing by cycle: {}{}", cycle_str.join(" -> "), cluster_info);
        }

        let monomers = decompose_by_cycle(
            &self.sequence,
            &self.hits,
            &self.cycle,
            &self.cluster_centers,
            verbose,
        );

        // Verify reconstruction
        let reconstructed: String = monomers.join("");
        if reconstructed != self.sequence {
            if verbose {
                eprintln!("WARNING: Reconstruction mismatch! {} vs {}",
                    self.sequence.len(), reconstructed.len());
            }
            return None;
        } else if verbose {
            eprintln!("Reconstruction: PERFECT");
        }

        // Calculate and report ED statistics in verbose mode
        if verbose && monomers.len() > 1 {
            let cut_seq = &self.cycle[0];
            // Calculate flank threshold
            let monomer_lengths: Vec<usize> = monomers.iter()
                .filter(|m| m.starts_with(cut_seq))
                .map(|m| m.len())
                .collect();
            let flank_threshold = if !monomer_lengths.is_empty() {
                let avg: f64 = monomer_lengths.iter().sum::<usize>() as f64 / monomer_lengths.len() as f64;
                (avg * 0.7) as usize
            } else {
                (self.estimated_period * 0.7) as usize
            };

            // Identify true monomers (not flanks) and calculate ED between consecutive pairs
            let mut eds: Vec<usize> = Vec::new();
            let mut prev_monomer: Option<&String> = None;
            for (i, monomer) in monomers.iter().enumerate() {
                let is_left_flank = i == 0 && !monomer.starts_with(cut_seq);
                let is_right_flank = i == monomers.len() - 1 && monomer.len() < flank_threshold;

                if !is_left_flank && !is_right_flank {
                    if let Some(prev) = prev_monomer {
                        let ed = levenshtein_exp(prev.as_bytes(), monomer.as_bytes());
                        eds.push(ed as usize);
                    }
                    prev_monomer = Some(monomer);
                }
            }

            if !eds.is_empty() {
                let sum_ed: usize = eds.iter().sum();
                let mean_ed: f64 = sum_ed as f64 / eds.len() as f64;
                let max_ed = *eds.iter().max().unwrap();
                let min_ed = *eds.iter().min().unwrap();
                eprintln!("ED between consecutive monomers: sum={}, mean={:.1}, min={}, max={}, n_pairs={}",
                    sum_ed, mean_ed, min_ed, max_ed, eds.len());
            }
        }

        Some(GraphDecomposition {
            monomers,
            cut_sequence: self.cycle[0].clone(),
            period: self.estimated_period as usize,
            cv: self.estimated_cv,
            alternative_anchors: self.alternative_anchors.clone(),
        })
    }

    /// Get decomposer statistics.
    pub fn get_stats(&self) -> HashMap<String, serde_json::Value> {
        use serde_json::json;

        let mut stats = HashMap::new();
        stats.insert("num_anchors".to_string(), json!(self.anchors.len()));
        stats.insert("anchors".to_string(), json!(self.anchors));
        stats.insert("num_hits".to_string(), json!(self.hits.len()));
        stats.insert("cycle".to_string(), json!(self.cycle));
        stats.insert("estimated_monomer_length".to_string(), json!(self.estimated_period));
        stats.insert("estimated_cv".to_string(), json!(self.estimated_cv));
        stats.insert("n_clusters".to_string(), json!(self.n_clusters));
        stats.insert("cluster_centers".to_string(), json!(self.cluster_centers));

        stats
    }
}

impl Default for AnchorGraphDecomposer {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_kmeans() {
        let data = vec![100.0, 100.0, 100.0, 200.0, 200.0, 200.0];
        let (clusters, centers) = simple_kmeans(&data, 2, 100);

        assert_eq!(clusters.len(), 2);
        // One cluster should have 100s, another 200s
        let has_100s = clusters.iter().any(|c| c.iter().all(|&x| (x - 100.0).abs() < 1.0));
        let has_200s = clusters.iter().any(|c| c.iter().all(|&x| (x - 200.0).abs() < 1.0));
        assert!(has_100s);
        assert!(has_200s);
    }

    #[test]
    fn test_find_anchor_hits_hierarchical() {
        let sequence = "ACGTACGTACGTACGT";
        let anchors = vec!["ACGT".to_string(), "ACG".to_string()];
        let hits = find_anchor_hits_hierarchical(sequence, &anchors, false);

        // ACGT should be found 4 times, ACG should not be found (covered by ACGT)
        assert!(!hits.is_empty());
        let acgt_count = hits.iter().filter(|h| h.anchor == "ACGT").count();
        assert!(acgt_count > 0);
    }

    #[test]
    fn test_anchor_graph_decomposer() {
        let sequence = "ACGTACGTACGTACGT";
        let candidates = vec!["ACGT".to_string()];

        let mut decomposer = AnchorGraphDecomposer::new();
        decomposer.build_from_candidates(sequence, &candidates, false);

        let result = decomposer.decompose(false);
        assert!(result.is_some());

        let decomp = result.unwrap();
        assert_eq!(decomp.cut_sequence, "ACGT");

        // Verify reconstruction
        let reconstructed: String = decomp.monomers.join("");
        assert_eq!(reconstructed, sequence);
    }
}
