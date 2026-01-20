//! Frequency Suffix Tree for finding frequent patterns in sequences
//!
//! Used to identify potential cut sequences for array decomposition.

use rustc_hash::FxHashMap;
use std::collections::HashMap;

/// A hint represents a frequent pattern found in the sequence
#[derive(Debug, Clone)]
pub struct Hint {
    pub pattern: String,
    pub positions: Vec<usize>,
    pub frequency: usize,
    pub score: f64,
}

impl Hint {
    pub fn new(pattern: String, positions: Vec<usize>) -> Self {
        let frequency = positions.len();
        Self {
            pattern,
            positions,
            frequency,
            score: 0.0,
        }
    }

    /// Calculate score based on frequency and pattern properties
    pub fn calculate_score(&mut self, seq_len: usize) {
        if self.positions.len() < 2 {
            self.score = 0.0;
            return;
        }

        // Calculate distances between consecutive occurrences
        let distances: Vec<usize> = self.positions
            .windows(2)
            .map(|w| w[1] - w[0])
            .collect();

        if distances.is_empty() {
            self.score = 0.0;
            return;
        }

        // Score based on consistency of distances (low variance = good)
        let mean_dist: f64 = distances.iter().sum::<usize>() as f64 / distances.len() as f64;
        let variance: f64 = distances.iter()
            .map(|&d| (d as f64 - mean_dist).powi(2))
            .sum::<f64>() / distances.len() as f64;

        let cv = if mean_dist > 0.0 {
            variance.sqrt() / mean_dist
        } else {
            f64::MAX
        };

        // Lower CV = more regular pattern = higher score
        // Also factor in frequency
        let coverage = self.positions.len() as f64 * mean_dist / seq_len as f64;
        self.score = coverage / (1.0 + cv);
    }
}

/// Frequency Suffix Tree implementation using HashMap
///
/// For each suffix length, stores the frequency of each pattern.
pub struct FsTree {
    /// Pattern -> positions mapping
    patterns: FxHashMap<String, Vec<usize>>,
    /// Sequence length
    seq_len: usize,
    /// Minimum pattern length
    min_len: usize,
    /// Maximum pattern length
    max_len: usize,
}

impl FsTree {
    /// Build frequency tree from sequence
    ///
    /// # Arguments
    /// * `sequence` - DNA sequence
    /// * `min_len` - Minimum pattern length to consider
    /// * `max_len` - Maximum pattern length (depth)
    /// * `cutoff` - Minimum frequency to keep a pattern
    pub fn new(sequence: &str, min_len: usize, max_len: usize, cutoff: usize) -> Self {
        let mut patterns: FxHashMap<String, Vec<usize>> = FxHashMap::default();
        let seq_bytes = sequence.as_bytes();
        let seq_len = sequence.len();

        // Scan all positions
        for pos in 0..seq_len {
            // For each pattern length
            for len in min_len..=max_len.min(seq_len - pos) {
                let pattern = &sequence[pos..pos + len];

                // Skip patterns with N
                if pattern.contains('N') || pattern.contains('n') {
                    continue;
                }

                patterns
                    .entry(pattern.to_string())
                    .or_insert_with(Vec::new)
                    .push(pos);
            }
        }

        // Filter by cutoff
        patterns.retain(|_, positions| positions.len() >= cutoff);

        Self {
            patterns,
            seq_len,
            min_len,
            max_len,
        }
    }

    /// Get all patterns with their positions
    pub fn get_patterns(&self) -> &FxHashMap<String, Vec<usize>> {
        &self.patterns
    }

    /// Get hints sorted by score
    pub fn get_hints(&self) -> Vec<Hint> {
        let mut hints: Vec<Hint> = self.patterns
            .iter()
            .map(|(pattern, positions)| {
                let mut hint = Hint::new(pattern.clone(), positions.clone());
                hint.calculate_score(self.seq_len);
                hint
            })
            .collect();

        // Sort by score descending
        hints.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap_or(std::cmp::Ordering::Equal));

        hints
    }

    /// Get hints filtered by top1 nucleotide start
    pub fn get_hints_starting_with(&self, nucleotide: char) -> Vec<Hint> {
        self.get_hints()
            .into_iter()
            .filter(|h| h.pattern.starts_with(nucleotide))
            .collect()
    }
}

/// Build frequency suffix tree with automatic cutoff calculation
pub fn get_fs_tree(array: &str, top1_nucleotide: char, depth: usize) -> FsTree {
    // Calculate cutoff based on sequence length
    // Expect at least 3 occurrences per kb
    let cutoff = (array.len() / 1000).max(3);

    FsTree::new(array, 3, depth, cutoff)
}

/// Iterate through hints and find best candidates for cut sequences
pub fn iterate_hints(
    array: &str,
    fs_tree: &FsTree,
    depth: usize,
    top1_nucleotide: char,
) -> Vec<Hint> {
    use crate::sequence::is_self_repeating;

    let mut hints = fs_tree.get_hints_starting_with(top1_nucleotide);

    // Filter out self-repeating patterns
    hints.retain(|h| !is_self_repeating(&h.pattern));

    // Filter by minimum length
    hints.retain(|h| h.pattern.len() >= 3);

    // Keep top candidates
    hints.truncate(100);

    hints
}

/// Analyze distances between pattern occurrences to estimate period
pub fn analyze_distances(positions: &[usize]) -> Option<(f64, f64, usize)> {
    if positions.len() < 2 {
        return None;
    }

    let distances: Vec<usize> = positions
        .windows(2)
        .map(|w| w[1] - w[0])
        .collect();

    if distances.is_empty() {
        return None;
    }

    let mean: f64 = distances.iter().sum::<usize>() as f64 / distances.len() as f64;
    let variance: f64 = distances.iter()
        .map(|&d| (d as f64 - mean).powi(2))
        .sum::<f64>() / distances.len() as f64;
    let std_dev = variance.sqrt();

    // Most common distance
    let mut dist_counts: HashMap<usize, usize> = HashMap::new();
    for &d in &distances {
        *dist_counts.entry(d).or_insert(0) += 1;
    }

    let mode = dist_counts
        .into_iter()
        .max_by_key(|(_, count)| *count)
        .map(|(dist, _)| dist)
        .unwrap_or(mean as usize);

    Some((mean, std_dev, mode))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fs_tree_simple() {
        let seq = "ACGTACGTACGTACGT";
        let tree = FsTree::new(seq, 3, 10, 2);

        let patterns = tree.get_patterns();
        assert!(patterns.contains_key("ACGT"));
        assert!(patterns.contains_key("CGT"));
    }

    #[test]
    fn test_hints() {
        let seq = "ACGTACGTACGTACGT";
        let tree = FsTree::new(seq, 3, 10, 2);
        let hints = tree.get_hints();

        assert!(!hints.is_empty());
        // ACGT should be a good pattern
        assert!(hints.iter().any(|h| h.pattern == "ACGT"));
    }

    #[test]
    fn test_analyze_distances() {
        let positions = vec![0, 100, 200, 300, 400];
        let (mean, std, mode) = analyze_distances(&positions).unwrap();

        assert!((mean - 100.0).abs() < 0.1);
        assert!(std < 1.0);
        assert_eq!(mode, 100);
    }
}
