//! ArraySplitter - De novo decomposition of satellite DNA arrays into monomers
//!
//! This is a Rust port of the Python ArraySplitter tool, designed for
//! high-performance parallel processing of satellite DNA arrays.
//!
//! ## Main Features
//!
//! - **De novo decomposition**: Automatically discovers cut sequences and decomposes arrays
//! - **Anchor graph algorithm**: Uses graph-based approach for optimal cut selection
//! - **Post-processing heuristics**: Refines decomposition for edge cases
//! - **Parallel processing**: Uses rayon for multi-threaded processing
//!
//! ## Example
//!
//! ```rust
//! use arraysplitter_rs::{decompose_array, apply_all_heuristics, Decomposition};
//!
//! let array = "ACGTACGTACGTACGT";
//! let result = decompose_array(array, 100, None, false);
//! println!("Found {} monomers with period {}bp", result.monomers.len(), result.period);
//! ```

pub mod sequence;
pub mod io;
pub mod fs_tree;
pub mod anchor_graph;
pub mod decompose;
pub mod heuristics;

// Re-exports for convenient access
pub use sequence::{
    get_revcomp,
    get_canonical_orientation,
    encode_kmer,
    decode_kmer,
    revcomp_kmer,
    canonical_kmer,
    is_self_repeating,
    gcd,
    find_gcd_of_list,
    count_nucleotides,
    get_top1_nucleotide,
    rotate_sequence,
    find_substring,
    find_all_positions,
};

pub use io::{
    FastaRecord,
    FastaReader,
    read_fasta,
    write_fasta,
    MonomerRecord,
    write_monomers_tsv,
    write_lengths,
    parse_cuts,
};

pub use fs_tree::{
    FsNode,
    Hint,
    is_self_repeating as fs_is_self_repeating,
    iter_fs_tree_from_sequence,
    iterate_hints,
    get_fs_tree_hints,
};

pub use anchor_graph::{
    AnchorHit,
    GraphDecomposition,
    AnchorGraphDecomposer,
};

pub use decompose::{
    Decomposition,
    decompose_array,
    decompose_array_with_cuts,
    is_canonical_orientation,
    rotate_monomers_to_cut,
};

pub use heuristics::{
    detect_and_fix_ab_pattern,
    optimize_monomer_lengths,
    find_mutant_anchor,
    split_long_monomers,
    split_duplicate_halves,
    split_by_popular_prefix,
    apply_all_heuristics,
};
