//! ArraySplitter CLI - De novo decomposition of satellite DNA arrays into monomers
//!
//! Usage:
//!   arraysplitter_rs -i input.fa -o output_prefix
//!   arraysplitter_rs -i input.fa -o output_prefix -c ATG,ATGATG,CGCG

use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use arraysplitter_rs::{
    decompose_array, decompose_array_with_cuts, apply_all_heuristics,
    read_fasta, parse_cuts, get_revcomp,
    decompose::is_canonical_orientation,
};

/// De novo decomposition of satellite DNA arrays into monomers
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input FASTA file
    #[arg(short, long)]
    input: String,

    /// Output prefix (will create .decomposed.fasta, .monomers.tsv, .lengths files)
    #[arg(short, long)]
    output: String,

    /// Input format: fasta, trf [default: fasta]
    #[arg(long, default_value = "fasta")]
    format: String,

    /// Number of threads [default: number of CPUs]
    #[arg(short, long)]
    threads: Option<usize>,

    /// Comma-separated list of predefined cut sequences (e.g., ATG,ATGATG)
    /// If provided, skips cut discovery
    #[arg(short, long)]
    cuts: Option<String>,

    /// Depth for hint discovery [default: 100]
    #[arg(short, long, default_value = "100")]
    depth: usize,

    /// Verbose output
    #[arg(short, long)]
    verbose: bool,
}

/// Result of processing a single array
struct ProcessedArray {
    header: String,
    decomposition: Vec<String>,
    cut_sequence: String,
    was_reversed: bool,
    period: usize,
}

/// Process a single array
fn process_array(
    header: &str,
    array: &str,
    predefined_cuts: &Option<Vec<String>>,
    depth: usize,
    verbose: bool,
) -> ProcessedArray {
    // Check canonical orientation
    let is_canonical = is_canonical_orientation(array);
    let (working_array, was_reversed) = if is_canonical {
        (array.to_string(), false)
    } else {
        (get_revcomp(array), true)
    };

    // Decompose
    let result = if let Some(cuts) = predefined_cuts {
        decompose_array_with_cuts(&working_array, cuts, verbose)
    } else {
        decompose_array(&working_array, depth, None, verbose)
    };

    // Apply heuristics
    let decomposition = apply_all_heuristics(&result.monomers, &result.cut_sequence, verbose);

    // Verify reconstruction
    let reconstructed: String = decomposition.join("");
    if reconstructed != working_array {
        eprintln!(
            "WARNING: Reconstruction mismatch for {}: {} != {}",
            header,
            working_array.len(),
            reconstructed.len()
        );
    }

    ProcessedArray {
        header: header.to_string(),
        decomposition,
        cut_sequence: result.cut_sequence,
        was_reversed,
        period: result.period,
    }
}

fn main() {
    let args = Args::parse();

    // Configure thread pool
    if let Some(threads) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .expect("Failed to configure thread pool");
    }

    // Parse predefined cuts
    let predefined_cuts: Option<Vec<String>> = args.cuts.as_ref().map(|s| parse_cuts(s));

    // Check input file exists
    if !Path::new(&args.input).exists() {
        eprintln!("Error: Input file '{}' not found", args.input);
        std::process::exit(1);
    }

    // Read all sequences first
    eprintln!("Reading input file: {}", args.input);
    let records = match read_fasta(&args.input) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("Error reading input file: {}", e);
            std::process::exit(1);
        }
    };

    let total = records.len();
    eprintln!("Found {} sequences", total);

    // Setup output prefix
    let mut output_prefix = args.output.clone();
    if output_prefix.ends_with(".fasta") {
        output_prefix = output_prefix[..output_prefix.len() - 6].to_string();
    } else if output_prefix.ends_with(".fa") {
        output_prefix = output_prefix[..output_prefix.len() - 3].to_string();
    }

    let output_file = format!("{}.decomposed.fasta", output_prefix);
    let detail_file = format!("{}.monomers.tsv", output_prefix);
    let lengths_file = format!("{}.lengths", output_prefix);

    eprintln!("Output files:");
    eprintln!("  {}", output_file);
    eprintln!("  {}", detail_file);
    eprintln!("  {}", lengths_file);

    if predefined_cuts.is_some() {
        eprintln!("Using predefined cuts: {:?}", predefined_cuts.as_ref().unwrap());
    } else {
        eprintln!("Will discover cuts automatically (depth={})", args.depth);
    }

    // Setup progress bar
    let pb = ProgressBar::new(total as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("#>-"),
    );

    // Process arrays in parallel
    let depth = args.depth;
    let verbose = args.verbose;
    let cuts_clone = predefined_cuts.clone();

    let results: Vec<ProcessedArray> = records
        .par_iter()
        .map(|record| {
            let result = process_array(
                &record.id,
                &record.sequence,
                &cuts_clone,
                depth,
                verbose,
            );
            pb.inc(1);
            result
        })
        .collect();

    pb.finish_with_message("Done!");

    // Write output files
    eprintln!("\nWriting output files...");

    // Open output files
    let mut fw = BufWriter::new(File::create(&output_file).expect("Failed to create output file"));
    let mut fw_detail =
        BufWriter::new(File::create(&detail_file).expect("Failed to create detail file"));
    let mut fw_lengths =
        BufWriter::new(File::create(&lengths_file).expect("Failed to create lengths file"));

    // Write header for detail file
    writeln!(
        fw_detail,
        "sequence_id\torientation\tindex\ttype\tlength\tis_flank\tsequence"
    )
    .unwrap();

    for result in &results {
        let decomposition = &result.decomposition;
        let cut_seq = &result.cut_sequence;
        let was_reversed = result.was_reversed;

        // Calculate statistics for internal monomers
        let all_monomer_lengths: Vec<usize> = decomposition
            .iter()
            .filter(|m| m.starts_with(cut_seq))
            .map(|m| m.len())
            .collect();

        let flank_threshold = if !all_monomer_lengths.is_empty() {
            let avg: f64 =
                all_monomer_lengths.iter().sum::<usize>() as f64 / all_monomer_lengths.len() as f64;
            (avg * 0.7) as usize
        } else {
            result.period * 50 / 100
        };

        // Count internal monomers (excluding flanks)
        let mut internal_count = 0;
        let mut internal_lengths: Vec<usize> = Vec::new();
        for (i, m) in decomposition.iter().enumerate() {
            if m.starts_with(cut_seq) {
                if i == decomposition.len() - 1 && m.len() < flank_threshold {
                    continue;
                }
                internal_count += 1;
                internal_lengths.push(m.len());
            }
        }

        let orientation = if was_reversed { "rev" } else { "fwd" };

        let header_info = if !internal_lengths.is_empty() {
            let min_len = *internal_lengths.iter().min().unwrap();
            let max_len = *internal_lengths.iter().max().unwrap();
            let avg_len: f64 =
                internal_lengths.iter().sum::<usize>() as f64 / internal_lengths.len() as f64;
            format!(
                "{} cut={} orientation={} n_monomers={} range={}-{} avg={:.1}",
                result.header, cut_seq, orientation, internal_count, min_len, max_len, avg_len
            )
        } else {
            format!(
                "{} cut={} orientation={} n_monomers=0",
                result.header, cut_seq, orientation
            )
        };

        // Write decomposed FASTA
        writeln!(fw, ">{}", header_info).unwrap();
        writeln!(fw, "{}", decomposition.join(" ")).unwrap();

        // Write lengths file
        writeln!(fw_lengths, ">{}", header_info).unwrap();
        let lengths: Vec<String> = decomposition.iter().map(|m| m.len().to_string()).collect();
        writeln!(fw_lengths, "{}", lengths.join(" ")).unwrap();

        // Write detailed monomer information
        for (i, monomer) in decomposition.iter().enumerate() {
            let (piece_type, is_flank) = if i == 0 && !monomer.starts_with(cut_seq) {
                ("LEFT_FLANK", "TRUE")
            } else if i == decomposition.len() - 1 && monomer.len() < flank_threshold {
                ("RIGHT_FLANK", "TRUE")
            } else {
                ("MONOMER", "FALSE")
            };

            writeln!(
                fw_detail,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                result.header,
                orientation,
                i,
                piece_type,
                monomer.len(),
                is_flank,
                monomer
            )
            .unwrap();
        }
    }

    eprintln!("Processing complete!");
    eprintln!("  {} sequences processed", total);
    eprintln!(
        "  {} total monomers",
        results.iter().map(|r| r.decomposition.len()).sum::<usize>()
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_process_array() {
        let result = process_array(
            "test",
            "ACGTACGTACGTACGT",
            &None,
            100,
            false,
        );

        assert_eq!(result.header, "test");
        assert!(result.decomposition.len() >= 2);

        // Verify reconstruction
        let reconstructed: String = result.decomposition.join("");
        assert_eq!(reconstructed, "ACGTACGTACGTACGT");
    }
}
