//! Anchor Graph based decomposition
//!
//! Uses anchor sequences to build a transition graph and find monomer cycles.

use rustc_hash::FxHashMap;
use std::collections::HashMap;

/// An anchor hit in the sequence
#[derive(Debug, Clone)]
pub struct AnchorHit {
    pub anchor: String,
    pub position: usize,
    pub distance_to_next: Option<usize>,
}

impl AnchorHit {
    pub fn new(anchor: String, position: usize) -> Self {
        Self {
            anchor,
            position,
            distance_to_next: None,
        }
    }
}

/// Result of anchor graph decomposition
#[derive(Debug, Clone)]
pub struct GraphDecomposition {
    pub monomers: Vec<String>,
    pub cut_sequence: String,
    pub period: usize,
    pub variance: f64,
    pub distances: Vec<usize>,
}

/// Anchor Graph Decomposer
///
/// Main class for anchor graph-based decomposition of satellite arrays.
pub struct AnchorGraphDecomposer {
    sequence: String,
    anchors: Vec<String>,
    hits: Vec<AnchorHit>,
    graph: FxHashMap<String, Vec<(String, usize)>>,  // from -> [(to, distance)]
}

impl AnchorGraphDecomposer {
    pub fn new() -> Self {
        Self {
            sequence: String::new(),
            anchors: Vec::new(),
            hits: Vec::new(),
            graph: FxHashMap::default(),
        }
    }

    /// Build decomposer from sequence and candidate anchors
    pub fn build_from_candidates(
        &mut self,
        sequence: &str,
        candidates: &[String],
        verbose: bool,
    ) {
        self.sequence = sequence.to_string();
        self.anchors = candidates.to_vec();

        // Find all anchor hits
        self.hits = self.find_anchor_hits(verbose);

        // Build transition graph
        self.build_transition_graph();
    }

    /// Find all anchor hits in the sequence
    fn find_anchor_hits(&self, verbose: bool) -> Vec<AnchorHit> {
        let mut all_hits: Vec<AnchorHit> = Vec::new();

        for anchor in &self.anchors {
            let mut pos = 0;
            while let Some(found_pos) = self.sequence[pos..].find(anchor) {
                let absolute_pos = pos + found_pos;
                all_hits.push(AnchorHit::new(anchor.clone(), absolute_pos));
                pos = absolute_pos + 1;
            }
        }

        // Sort by position
        all_hits.sort_by_key(|h| h.position);

        // Calculate distances to next hit
        for i in 0..all_hits.len().saturating_sub(1) {
            all_hits[i].distance_to_next = Some(all_hits[i + 1].position - all_hits[i].position);
        }

        if verbose {
            eprintln!("  Found {} anchor hits", all_hits.len());
        }

        all_hits
    }

    /// Build transition graph from anchor hits
    fn build_transition_graph(&mut self) {
        self.graph.clear();

        for i in 0..self.hits.len().saturating_sub(1) {
            let from = &self.hits[i].anchor;
            let to = &self.hits[i + 1].anchor;
            let distance = self.hits[i + 1].position - self.hits[i].position;

            self.graph
                .entry(from.clone())
                .or_insert_with(Vec::new)
                .push((to.clone(), distance));
        }
    }

    /// Find the most common transition cycle (monomer pattern)
    pub fn find_monomer_cycle(&self, verbose: bool) -> Option<(Vec<String>, usize, f64)> {
        if self.hits.len() < 3 {
            return None;
        }

        // Analyze distances to find period
        let distances: Vec<usize> = self.hits
            .iter()
            .filter_map(|h| h.distance_to_next)
            .collect();

        if distances.is_empty() {
            return None;
        }

        // Find most common distance (approximate period)
        let (mean, std, clusters) = self.analyze_distance_distribution(&distances);

        if verbose {
            eprintln!("  Distance analysis: mean={:.1}, std={:.1}", mean, std);
        }

        // Find cycle starting from most common anchor
        let mut anchor_counts: HashMap<&String, usize> = HashMap::new();
        for hit in &self.hits {
            *anchor_counts.entry(&hit.anchor).or_insert(0) += 1;
        }

        let start_anchor = anchor_counts
            .into_iter()
            .max_by_key(|(_, count)| *count)
            .map(|(anchor, _)| anchor.clone())?;

        // Simple cycle: just the starting anchor
        let cycle = vec![start_anchor.clone()];
        let period = mean as usize;
        let variance = std / mean;

        Some((cycle, period, variance))
    }

    /// Analyze distance distribution using simple clustering
    fn analyze_distance_distribution(&self, distances: &[usize]) -> (f64, f64, Vec<Vec<usize>>) {
        let mean: f64 = distances.iter().sum::<usize>() as f64 / distances.len() as f64;
        let variance: f64 = distances.iter()
            .map(|&d| (d as f64 - mean).powi(2))
            .sum::<f64>() / distances.len() as f64;
        let std = variance.sqrt();

        // Simple k-means with k=2 (period and 2*period typically)
        let clusters = self.simple_kmeans(distances, 2);

        (mean, std, clusters)
    }

    /// Simple k-means clustering for distance analysis
    fn simple_kmeans(&self, data: &[usize], k: usize) -> Vec<Vec<usize>> {
        if data.is_empty() || k == 0 {
            return vec![];
        }

        let mut data_f64: Vec<f64> = data.iter().map(|&x| x as f64).collect();
        data_f64.sort_by(|a, b| a.partial_cmp(b).unwrap());

        // Initialize centroids
        let mut centroids: Vec<f64> = (0..k)
            .map(|i| data_f64[i * data_f64.len() / k.max(1)])
            .collect();

        for _ in 0..100 {
            // Assign points to clusters
            let mut clusters: Vec<Vec<f64>> = vec![Vec::new(); k];

            for &val in &data_f64 {
                let closest = centroids
                    .iter()
                    .enumerate()
                    .min_by(|(_, a), (_, b)| {
                        (val - *a).abs().partial_cmp(&(val - *b).abs()).unwrap()
                    })
                    .map(|(i, _)| i)
                    .unwrap_or(0);

                clusters[closest].push(val);
            }

            // Update centroids
            let new_centroids: Vec<f64> = clusters
                .iter()
                .map(|c| {
                    if c.is_empty() {
                        0.0
                    } else {
                        c.iter().sum::<f64>() / c.len() as f64
                    }
                })
                .collect();

            // Check convergence
            let converged = centroids
                .iter()
                .zip(new_centroids.iter())
                .all(|(a, b)| (a - b).abs() < 1.0);

            centroids = new_centroids;

            if converged {
                break;
            }
        }

        // Final assignment
        let mut result: Vec<Vec<usize>> = vec![Vec::new(); k];
        for &val in data {
            let closest = centroids
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| {
                    (val as f64 - *a).abs().partial_cmp(&(val as f64 - *b).abs()).unwrap()
                })
                .map(|(i, _)| i)
                .unwrap_or(0);

            result[closest].push(val);
        }

        result
    }

    /// Decompose sequence using anchor graph
    pub fn decompose(&self, verbose: bool) -> Option<GraphDecomposition> {
        let (cycle, period, variance) = self.find_monomer_cycle(verbose)?;

        if cycle.is_empty() {
            return None;
        }

        let cut_sequence = cycle[0].clone();

        // Split sequence by cut positions
        let cut_positions: Vec<usize> = self.hits
            .iter()
            .filter(|h| h.anchor == cut_sequence)
            .map(|h| h.position)
            .collect();

        if cut_positions.is_empty() {
            return None;
        }

        // Extract monomers
        let mut monomers: Vec<String> = Vec::new();
        let mut distances: Vec<usize> = Vec::new();

        for i in 0..cut_positions.len() {
            let start = cut_positions[i];
            let end = if i + 1 < cut_positions.len() {
                cut_positions[i + 1]
            } else {
                self.sequence.len()
            };

            monomers.push(self.sequence[start..end].to_string());

            if i + 1 < cut_positions.len() {
                distances.push(cut_positions[i + 1] - cut_positions[i]);
            }
        }

        // Handle left flank
        if cut_positions[0] > 0 {
            let left_flank = self.sequence[..cut_positions[0]].to_string();
            monomers.insert(0, left_flank);
        }

        Some(GraphDecomposition {
            monomers,
            cut_sequence,
            period,
            variance,
            distances,
        })
    }

    /// Get statistics about the decomposition
    pub fn get_stats(&self) -> HashMap<String, f64> {
        let mut stats = HashMap::new();

        stats.insert("num_anchors".to_string(), self.anchors.len() as f64);
        stats.insert("num_hits".to_string(), self.hits.len() as f64);

        if let Some((_cycle, period, variance)) = self.find_monomer_cycle(false) {
            stats.insert("period".to_string(), period as f64);
            stats.insert("variance".to_string(), variance);
        }

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
    fn test_anchor_graph_simple() {
        let sequence = "ACGTACGTACGTACGT";
        let candidates = vec!["ACGT".to_string()];

        let mut decomposer = AnchorGraphDecomposer::new();
        decomposer.build_from_candidates(sequence, &candidates, false);

        let result = decomposer.decompose(false);
        assert!(result.is_some());

        let decomp = result.unwrap();
        assert_eq!(decomp.cut_sequence, "ACGT");
        assert_eq!(decomp.period, 4);
    }

    #[test]
    fn test_kmeans() {
        let decomposer = AnchorGraphDecomposer::new();
        let data = vec![100, 100, 100, 200, 200, 200];
        let clusters = decomposer.simple_kmeans(&data, 2);

        assert_eq!(clusters.len(), 2);
        // One cluster should have 100s, another 200s
        assert!(clusters.iter().any(|c| c.iter().all(|&x| x == 100)));
        assert!(clusters.iter().any(|c| c.iter().all(|&x| x == 200)));
    }
}
