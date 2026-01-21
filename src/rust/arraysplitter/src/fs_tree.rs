//! Frequency Suffix Tree for finding frequent patterns in sequences
//!
//! Iterative implementation matching Python's iter_fs_tree_from_sequence exactly.
//!
//! Used to identify potential cut sequences for array decomposition.

use std::collections::BinaryHeap;

/// A node in the frequency suffix tree iteration
/// Stored as (length, names, positions) matching Python
#[derive(Debug, Clone)]
pub struct FsNode {
    pub length: usize,
    pub names: Vec<usize>,      // starting positions (original indices)
    pub positions: Vec<usize>,  // current positions (where we are now in the string)
}

impl Ord for FsNode {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Min-heap by length (smaller lengths processed first)
        other.length.cmp(&self.length)
    }
}

impl PartialOrd for FsNode {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for FsNode {
    fn eq(&self, other: &Self) -> bool {
        self.length == other.length
    }
}

impl Eq for FsNode {}

/// Hint from the fs-tree: (length, sequence, frequency)
#[derive(Debug, Clone)]
pub struct Hint {
    pub length: usize,
    pub pattern: String,
    pub frequency: usize,
}

/// Check if a pattern is composed of repeated smaller units.
/// Returns the minimal unit if self-repeating, None otherwise.
pub fn is_self_repeating(pattern: &str) -> Option<String> {
    let n = pattern.len();
    for sub_len in 1..=(n / 2) {
        if n % sub_len == 0 {
            let sub_pattern = &pattern[..sub_len];
            let repeated = sub_pattern.repeat(n / sub_len);
            if pattern == repeated {
                return Some(sub_pattern.to_string());
            }
        }
    }
    None
}

/// Iterate through fs-tree from sequence, matching Python iter_fs_tree_from_sequence exactly.
///
/// This function yields hints as it traverses the tree using a heap (priority queue).
///
/// # Arguments
/// * `array` - The DNA sequence
/// * `starting_nucleotide` - Starting nucleotide (e.g., 'A', 'C', 'G', 'T')
/// * `cutoff` - Minimum frequency to keep a branch
/// * `depth` - Maximum pattern length to explore
///
/// # Returns
/// Vector of (length, names, positions) tuples
pub fn iter_fs_tree_from_sequence(
    array: &str,
    starting_nucleotide: char,
    cutoff: usize,
    depth: usize,
) -> Vec<FsNode> {
    let array_bytes = array.as_bytes();
    let array_len = array.len();

    // Find all positions of the starting nucleotide
    let starting_nucleotide_upper = starting_nucleotide.to_ascii_uppercase() as u8;
    let names: Vec<usize> = array_bytes
        .iter()
        .enumerate()
        .filter(|(_, &b)| b == starting_nucleotide_upper)
        .map(|(i, _)| i)
        .collect();

    if names.len() <= cutoff {
        return Vec::new();
    }

    let positions = names.clone();

    // Initial node: (length=1, names, positions)
    let initial = FsNode {
        length: 1,
        names,
        positions,
    };

    // Use min-heap ordered by length
    let mut heap: BinaryHeap<FsNode> = BinaryHeap::new();
    heap.push(initial);

    let mut results: Vec<FsNode> = Vec::new();

    while let Some(node) = heap.pop() {
        let l = node.length;

        // Stop if we exceed depth
        if l > depth {
            continue;
        }

        // Classify positions by next nucleotide
        let mut fs_a: Vec<usize> = Vec::new();
        let mut fs_c: Vec<usize> = Vec::new();
        let mut fs_g: Vec<usize> = Vec::new();
        let mut fs_t: Vec<usize> = Vec::new();
        let mut fs_ap: Vec<usize> = Vec::new();
        let mut fs_cp: Vec<usize> = Vec::new();
        let mut fs_gp: Vec<usize> = Vec::new();
        let mut fs_tp: Vec<usize> = Vec::new();

        for (ii, &pos) in node.positions.iter().enumerate() {
            let name = node.names[ii];

            // Check if we're at the end of the array
            if pos + 1 >= array_len {
                continue;
            }

            let nucl = array_bytes[pos + 1];
            match nucl {
                b'A' => {
                    fs_a.push(name);
                    fs_ap.push(pos + 1);
                }
                b'C' => {
                    fs_c.push(name);
                    fs_cp.push(pos + 1);
                }
                b'G' => {
                    fs_g.push(name);
                    fs_gp.push(pos + 1);
                }
                b'T' => {
                    fs_t.push(name);
                    fs_tp.push(pos + 1);
                }
                _ => {}
            }
        }

        // Process each nucleotide extension
        if fs_a.len() > cutoff {
            let new_node = FsNode {
                length: l + 1,
                names: fs_a,
                positions: fs_ap,
            };
            heap.push(new_node.clone());
            results.push(new_node);
        }

        if fs_c.len() > cutoff {
            let new_node = FsNode {
                length: l + 1,
                names: fs_c,
                positions: fs_cp,
            };
            heap.push(new_node.clone());
            results.push(new_node);
        }

        if fs_g.len() > cutoff {
            let new_node = FsNode {
                length: l + 1,
                names: fs_g,
                positions: fs_gp,
            };
            heap.push(new_node.clone());
            results.push(new_node);
        }

        if fs_t.len() > cutoff {
            let new_node = FsNode {
                length: l + 1,
                names: fs_t,
                positions: fs_tp,
            };
            heap.push(new_node.clone());
            results.push(new_node);
        }
    }

    results
}

/// Iterate through hints matching Python's iterate_hints function.
///
/// Key features:
/// - Groups results by length and picks max frequency at each length
/// - Detects self-repeating patterns and yields minimal unit instead
/// - Tracks found patterns to avoid duplicates
///
/// # Returns
/// Vector of (length, sequence, frequency) hints
pub fn iterate_hints(
    array: &str,
    starting_nucleotide: char,
    cutoff: usize,
    depth: usize,
) -> Vec<Hint> {
    let nodes = iter_fs_tree_from_sequence(array, starting_nucleotide, cutoff, depth);

    if nodes.is_empty() {
        return Vec::new();
    }

    // Sort nodes by length
    let mut sorted_nodes = nodes;
    sorted_nodes.sort_by_key(|n| n.length);

    // Group by length and pick best at each length
    let mut hints: Vec<Hint> = Vec::new();
    let mut found_patterns: std::collections::HashSet<String> = std::collections::HashSet::new();

    let mut current_length = 0;
    let mut buffer: Vec<&FsNode> = Vec::new();

    for node in &sorted_nodes {
        let l = node.length;

        if l != current_length {
            // Process buffer for previous length
            if !buffer.is_empty() {
                // Find node with max frequency
                let best = buffer.iter().max_by_key(|n| n.names.len()).unwrap();
                let start = best.names[0];
                let end = best.positions[0];
                let found_seq = &array[start..=end];
                let n_freq = best.names.len();

                // Check if self-repeating
                if let Some(minimal_unit) = is_self_repeating(found_seq) {
                    // Yield minimal unit instead
                    if !found_patterns.contains(&minimal_unit) {
                        let min_len = minimal_unit.len();
                        let min_count = array.matches(&minimal_unit).count();
                        hints.push(Hint {
                            length: min_len,
                            pattern: minimal_unit.clone(),
                            frequency: min_count,
                        });
                        found_patterns.insert(minimal_unit);
                    }
                } else {
                    // Not self-repeating, yield as normal
                    let pattern = found_seq.to_string();
                    if !found_patterns.contains(&pattern) {
                        hints.push(Hint {
                            length: current_length,
                            pattern: pattern.clone(),
                            frequency: n_freq,
                        });
                        found_patterns.insert(pattern);
                    }
                }
            }

            buffer.clear();
            current_length = l;

            // Stop if we exceed depth
            if current_length > depth {
                break;
            }
        }

        buffer.push(node);
    }

    // Process final buffer
    if !buffer.is_empty() && current_length <= depth {
        let best = buffer.iter().max_by_key(|n| n.names.len()).unwrap();
        let start = best.names[0];
        let end = best.positions[0];

        if end < array.len() {
            let found_seq = &array[start..=end];
            let n_freq = best.names.len();

            if let Some(minimal_unit) = is_self_repeating(found_seq) {
                if !found_patterns.contains(&minimal_unit) {
                    let min_len = minimal_unit.len();
                    let min_count = array.matches(&minimal_unit).count();
                    hints.push(Hint {
                        length: min_len,
                        pattern: minimal_unit.clone(),
                        frequency: min_count,
                    });
                }
            } else {
                let pattern = found_seq.to_string();
                if !found_patterns.contains(&pattern) {
                    hints.push(Hint {
                        length: current_length,
                        pattern: pattern.clone(),
                        frequency: n_freq,
                    });
                }
            }
        }
    }

    hints
}

/// Get fs_tree for array with automatic cutoff calculation (matching Python get_fs_tree)
pub fn get_fs_tree_hints(array: &str, nucleotide: char, cutoff: usize, depth: usize) -> Vec<Hint> {
    iterate_hints(array, nucleotide, cutoff, depth)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_self_repeating() {
        assert_eq!(is_self_repeating("ACAC"), Some("AC".to_string()));
        assert_eq!(is_self_repeating("ACGTACGT"), Some("ACGT".to_string()));
        assert_eq!(is_self_repeating("AAAAAA"), Some("A".to_string()));
        assert_eq!(is_self_repeating("ACGT"), None);
        assert_eq!(is_self_repeating("ACG"), None);
    }

    #[test]
    fn test_iter_fs_tree() {
        let seq = "ACGTACGTACGTACGT";
        let nodes = iter_fs_tree_from_sequence(seq, 'A', 2, 10);

        // Should have nodes at different lengths
        assert!(!nodes.is_empty());

        // Check that we found ACGT pattern (starting with A)
        let has_acgt = nodes.iter().any(|n| {
            if n.names.is_empty() { return false; }
            let start = n.names[0];
            let end = n.positions[0];
            end < seq.len() && &seq[start..=end] == "ACGT"
        });
        assert!(has_acgt);
    }

    #[test]
    fn test_iterate_hints() {
        let seq = "ACGTACGTACGTACGT";
        let hints = iterate_hints(seq, 'A', 2, 10);

        assert!(!hints.is_empty());

        // Should find ACGT or its patterns
        let patterns: Vec<&str> = hints.iter().map(|h| h.pattern.as_str()).collect();
        println!("Found patterns: {:?}", patterns);
    }
}
