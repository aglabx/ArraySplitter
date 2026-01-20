//! Sequence utilities for DNA manipulation
//!
//! Provides functions for reverse complement, canonical orientation,
//! k-mer encoding, and other sequence operations.

use std::collections::HashMap;

/// Get reverse complement of a DNA sequence
pub fn get_revcomp(sequence: &str) -> String {
    sequence
        .chars()
        .rev()
        .map(|c| match c {
            'A' | 'a' => 'T',
            'T' | 't' => 'A',
            'C' | 'c' => 'G',
            'G' | 'g' => 'C',
            'N' | 'n' => 'N',
            _ => c,
        })
        .collect()
}

/// Get canonical orientation of a sequence
///
/// Returns the sequence in canonical form where:
/// - First differing position between forward and revcomp determines orientation
/// - A < T and C < G (so we prefer A over T, C over G)
///
/// Returns (canonical_sequence, is_forward)
pub fn get_canonical_orientation(sequence: &str) -> (String, bool) {
    let seq_upper: String = sequence.to_uppercase();
    let revcomp = get_revcomp(&seq_upper);

    // Compare character by character
    for (fwd, rev) in seq_upper.chars().zip(revcomp.chars()) {
        if fwd != rev {
            // Use lexicographic ordering: A < C < G < T
            // But we want A < T and C < G, so:
            // A(0) < C(1) < G(2) < T(3) in standard order
            // We prefer: A over T, C over G
            let fwd_val = match fwd {
                'A' => 0,
                'C' => 1,
                'G' => 2,
                'T' => 3,
                _ => 4,
            };
            let rev_val = match rev {
                'A' => 0,
                'C' => 1,
                'G' => 2,
                'T' => 3,
                _ => 4,
            };

            if fwd_val <= rev_val {
                return (seq_upper, true);
            } else {
                return (revcomp, false);
            }
        }
    }

    // Sequences are identical (palindrome)
    (seq_upper, true)
}

/// Encode a k-mer as a u64 (2 bits per nucleotide)
///
/// Maximum k = 32 for u64
pub fn encode_kmer(kmer: &str) -> u64 {
    let mut val: u64 = 0;
    for c in kmer.chars() {
        val <<= 2;
        val |= match c {
            'A' | 'a' => 0,
            'C' | 'c' => 1,
            'G' | 'g' => 2,
            'T' | 't' => 3,
            _ => 0, // N treated as A
        };
    }
    val
}

/// Decode a k-mer from u64 back to string
pub fn decode_kmer(val: u64, k: usize) -> String {
    let mut v = val;
    let mut chars: Vec<char> = Vec::with_capacity(k);

    for _ in 0..k {
        let base = match v & 3 {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => 'N',
        };
        chars.push(base);
        v >>= 2;
    }

    chars.reverse();
    chars.into_iter().collect()
}

/// Get reverse complement of encoded k-mer
pub fn revcomp_kmer(kmer: u64, k: usize) -> u64 {
    let mut rc: u64 = 0;
    let mut val = kmer;

    for _ in 0..k {
        rc <<= 2;
        rc |= (val & 3) ^ 3; // Complement: A(0)↔T(3), C(1)↔G(2)
        val >>= 2;
    }

    rc
}

/// Get canonical k-mer (min of forward and revcomp)
pub fn canonical_kmer(kmer: u64, k: usize) -> u64 {
    let rc = revcomp_kmer(kmer, k);
    kmer.min(rc)
}

/// Count nucleotide frequencies in a sequence
pub fn count_nucleotides(sequence: &str) -> HashMap<char, usize> {
    let mut counts: HashMap<char, usize> = HashMap::new();

    for c in sequence.chars() {
        let c_upper = c.to_ascii_uppercase();
        *counts.entry(c_upper).or_insert(0) += 1;
    }

    counts
}

/// Get the most frequent nucleotide
pub fn get_top1_nucleotide(sequence: &str) -> char {
    let counts = count_nucleotides(sequence);

    counts
        .into_iter()
        .filter(|(c, _)| *c == 'A' || *c == 'C' || *c == 'G' || *c == 'T')
        .max_by_key(|(_, count)| *count)
        .map(|(c, _)| c)
        .unwrap_or('A')
}

/// Check if a pattern is self-repeating (e.g., "ATAT" -> "AT" repeated)
pub fn is_self_repeating(pattern: &str) -> bool {
    if pattern.len() < 2 {
        return false;
    }

    // Check if pattern can be expressed as repetition of a smaller unit
    for unit_len in 1..=pattern.len() / 2 {
        if pattern.len() % unit_len == 0 {
            let unit = &pattern[..unit_len];
            let repeated: String = unit.repeat(pattern.len() / unit_len);
            if repeated == pattern {
                return true;
            }
        }
    }

    false
}

/// Compute GCD of two numbers
pub fn gcd(a: usize, b: usize) -> usize {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

/// Find GCD of a list of numbers
pub fn find_gcd_of_list(numbers: &[usize]) -> usize {
    if numbers.is_empty() {
        return 0;
    }

    numbers.iter().fold(numbers[0], |acc, &x| gcd(acc, x))
}

/// Rotate a sequence to start from a given position
pub fn rotate_sequence(sequence: &str, start_pos: usize) -> String {
    if sequence.is_empty() || start_pos >= sequence.len() {
        return sequence.to_string();
    }

    format!("{}{}", &sequence[start_pos..], &sequence[..start_pos])
}

/// Find position of a substring in a sequence
pub fn find_substring(sequence: &str, pattern: &str) -> Option<usize> {
    sequence.find(pattern)
}

/// Find all positions of a substring in a sequence
pub fn find_all_positions(sequence: &str, pattern: &str) -> Vec<usize> {
    let mut positions = Vec::new();
    let mut start = 0;

    while let Some(pos) = sequence[start..].find(pattern) {
        positions.push(start + pos);
        start += pos + 1;
    }

    positions
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_revcomp() {
        assert_eq!(get_revcomp("ACGT"), "ACGT");
        assert_eq!(get_revcomp("AAAAAA"), "TTTTTT");
        assert_eq!(get_revcomp("ATCGATCG"), "CGATCGAT");
    }

    #[test]
    fn test_canonical_orientation() {
        let (canon, is_fwd) = get_canonical_orientation("TTTT");
        assert_eq!(canon, "AAAA");
        assert!(!is_fwd);

        let (canon, is_fwd) = get_canonical_orientation("ACGT");
        assert_eq!(canon, "ACGT");
        assert!(is_fwd);
    }

    #[test]
    fn test_encode_decode_kmer() {
        let kmer = "ACGT";
        let encoded = encode_kmer(kmer);
        let decoded = decode_kmer(encoded, 4);
        assert_eq!(decoded, kmer);
    }

    #[test]
    fn test_is_self_repeating() {
        assert!(is_self_repeating("ATAT"));
        assert!(is_self_repeating("ABCABC"));
        assert!(!is_self_repeating("ACGT"));
        assert!(!is_self_repeating("A"));
    }

    #[test]
    fn test_gcd() {
        assert_eq!(gcd(12, 8), 4);
        assert_eq!(gcd(100, 25), 25);
        assert_eq!(find_gcd_of_list(&[12, 18, 24]), 6);
    }
}
