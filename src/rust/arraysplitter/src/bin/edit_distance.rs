//! Edit Distance Calculator for Monomer TSV files
//!
//! Calculates edit distance between consecutive monomers within each sequence.
//! Architecture: Reader -> Workers -> Writer (bounded channels, constant memory)

use clap::Parser;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::sync::mpsc::{self, Receiver, SyncSender};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::thread;
use triple_accel::levenshtein_exp;

/// Calculate edit distance between consecutive monomers
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input monomers TSV file
    #[arg(short, long)]
    input: String,

    /// Output TSV file with edit distances
    #[arg(short, long)]
    output: String,

    /// Number of worker threads [default: number of CPUs]
    #[arg(short, long)]
    threads: Option<usize>,

    /// Maximum sequence length for edit distance calculation [default: 10000]
    /// Sequences longer than this will have "TOO_LONG" instead of edit distance
    #[arg(short, long, default_value = "10000")]
    max_length: usize,

    /// Verbose output
    #[arg(short, long)]
    verbose: bool,
}

/// A monomer record from the TSV
#[derive(Debug, Clone)]
struct MonomerRecord {
    sequence_id: String,
    orientation: String,
    index: usize,
    piece_type: String,
    length: usize,
    is_flank: String,
    sequence: String,
}

/// A group of monomers for one sequence
struct SequenceGroup {
    sequence_id: String,
    monomers: Vec<MonomerRecord>,
}

/// Output record with edit distance
struct OutputRecord {
    sequence_id: String,
    records: Vec<(MonomerRecord, EditDistResult)>, // (monomer, edit_distance_to_prev)
}

/// Parse a TSV line into a MonomerRecord
fn parse_line(line: &str) -> Option<MonomerRecord> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 7 {
        return None;
    }

    Some(MonomerRecord {
        sequence_id: fields[0].to_string(),
        orientation: fields[1].to_string(),
        index: fields[2].parse().unwrap_or(0),
        piece_type: fields[3].to_string(),
        length: fields[4].parse().unwrap_or(0),
        is_flank: fields[5].to_string(),
        sequence: fields[6].to_string(),
    })
}

/// Calculate Levenshtein edit distance between two sequences
fn edit_distance(s1: &str, s2: &str) -> usize {
    levenshtein_exp(s1.as_bytes(), s2.as_bytes()) as usize
}

/// Reader thread: reads TSV and groups by sequence_id
fn reader_thread(
    input_path: String,
    tx: SyncSender<Option<SequenceGroup>>,
    num_workers: usize,
    total_count: Arc<AtomicUsize>,
) {
    let file = File::open(&input_path).expect("Failed to open input file");
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    // Skip header
    let _header = lines.next();

    let mut current_id: Option<String> = None;
    let mut current_group: Vec<MonomerRecord> = Vec::new();
    let mut count = 0;

    for line in lines {
        let line = line.expect("Failed to read line");
        if line.trim().is_empty() {
            continue;
        }

        if let Some(record) = parse_line(&line) {
            // Check if we're starting a new sequence
            if current_id.as_ref() != Some(&record.sequence_id) {
                // Send previous group if exists
                if let Some(id) = current_id.take() {
                    if !current_group.is_empty() {
                        tx.send(Some(SequenceGroup {
                            sequence_id: id,
                            monomers: std::mem::take(&mut current_group),
                        }))
                        .expect("Failed to send group");
                        count += 1;
                    }
                }
                current_id = Some(record.sequence_id.clone());
            }
            current_group.push(record);
        }
    }

    // Send last group
    if let Some(id) = current_id {
        if !current_group.is_empty() {
            tx.send(Some(SequenceGroup {
                sequence_id: id,
                monomers: current_group,
            }))
            .expect("Failed to send last group");
            count += 1;
        }
    }

    total_count.store(count, Ordering::SeqCst);

    // Send termination signals for all workers
    for _ in 0..num_workers {
        tx.send(None).expect("Failed to send termination");
    }
}

/// Edit distance result - either a value, too long, or first monomer
#[derive(Debug, Clone)]
enum EditDistResult {
    Value(usize),
    TooLong,
    First,
}

/// Process a sequence group: calculate edit distances
fn process_group(group: SequenceGroup, max_length: usize) -> OutputRecord {
    // Sort by index to ensure correct order
    let mut monomers = group.monomers;
    monomers.sort_by_key(|m| m.index);

    // Calculate edit distances
    let mut results: Vec<(MonomerRecord, EditDistResult)> = Vec::with_capacity(monomers.len());

    for i in 0..monomers.len() {
        let ed = if i > 0 {
            // Skip edit distance for very long sequences (O(n*m) complexity)
            if monomers[i].sequence.len() > max_length || monomers[i - 1].sequence.len() > max_length {
                EditDistResult::TooLong
            } else {
                EditDistResult::Value(edit_distance(&monomers[i - 1].sequence, &monomers[i].sequence))
            }
        } else {
            EditDistResult::First // First monomer has no previous
        };
        results.push((monomers[i].clone(), ed));
    }

    OutputRecord {
        sequence_id: group.sequence_id,
        records: results,
    }
}

/// Writer thread: receives results and writes to file
fn writer_thread(
    rx: Receiver<Option<OutputRecord>>,
    output_path: String,
    total_count: Arc<AtomicUsize>,
    num_workers: usize,
    verbose: bool,
) {
    let mut fw = BufWriter::new(File::create(&output_path).expect("Failed to create output file"));

    // Write header with new column (edit_dist before sequence for easier viewing)
    writeln!(
        fw,
        "sequence_id\torientation\tindex\ttype\tlength\tis_flank\tedit_dist_to_prev\tsequence"
    )
    .unwrap();

    let mut processed = 0;
    let mut finished_workers = 0;
    let mut total_monomers = 0;
    let mut total_edit_distance: usize = 0;
    let mut edit_distance_count = 0;
    let mut too_long_count = 0;

    // Statistics per sequence
    let mut sequence_stats: HashMap<String, Vec<usize>> = HashMap::new();

    loop {
        match rx.recv() {
            Ok(Some(result)) => {
                processed += 1;
                total_monomers += result.records.len();

                for (monomer, ed) in &result.records {
                    let ed_str = match ed {
                        EditDistResult::Value(d) => {
                            total_edit_distance += *d;
                            edit_distance_count += 1;
                            sequence_stats
                                .entry(result.sequence_id.clone())
                                .or_default()
                                .push(*d);
                            d.to_string()
                        }
                        EditDistResult::TooLong => {
                            too_long_count += 1;
                            "TOO_LONG".to_string()
                        }
                        EditDistResult::First => "NA".to_string(),
                    };

                    writeln!(
                        fw,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        monomer.sequence_id,
                        monomer.orientation,
                        monomer.index,
                        monomer.piece_type,
                        monomer.length,
                        monomer.is_flank,
                        ed_str,
                        monomer.sequence
                    )
                    .unwrap();
                }

                // Progress output
                let total = total_count.load(Ordering::SeqCst);
                if total > 0 && processed % 100 == 0 {
                    eprint!(
                        "\rProcessed: {}/{} ({:.1}%)",
                        processed,
                        total,
                        (processed as f64 / total as f64) * 100.0
                    );
                }
            }
            Ok(None) => {
                finished_workers += 1;
                if finished_workers >= num_workers {
                    break;
                }
            }
            Err(_) => break,
        }
    }

    let total = total_count.load(Ordering::SeqCst);
    eprintln!("\rProcessed: {}/{} (100.0%)", processed, total);
    eprintln!("Total monomers: {}", total_monomers);
    eprintln!("Total edit distance pairs: {}", edit_distance_count);
    if too_long_count > 0 {
        eprintln!("Skipped (too long): {}", too_long_count);
    }

    if edit_distance_count > 0 {
        let mean_ed = total_edit_distance as f64 / edit_distance_count as f64;
        eprintln!("Mean edit distance: {:.2}", mean_ed);

        if verbose {
            // Calculate variance
            let mut variance_sum = 0.0;
            for eds in sequence_stats.values() {
                for &ed in eds {
                    variance_sum += (ed as f64 - mean_ed).powi(2);
                }
            }
            let variance = variance_sum / edit_distance_count as f64;
            let std_dev = variance.sqrt();
            eprintln!("Std deviation: {:.2}", std_dev);
            eprintln!("CV: {:.2}%", (std_dev / mean_ed) * 100.0);
        }
    }
}

fn main() {
    let args = Args::parse();

    let num_workers = args.threads.unwrap_or_else(num_cpus::get);
    eprintln!("Using {} worker threads", num_workers);

    // Check input file
    if !Path::new(&args.input).exists() {
        eprintln!("Error: Input file '{}' not found", args.input);
        std::process::exit(1);
    }

    eprintln!("Input: {}", args.input);
    eprintln!("Output: {}", args.output);
    eprintln!("Max sequence length: {}", args.max_length);

    // Channels with bounded capacity (controls memory)
    let (input_tx, input_rx) = mpsc::sync_channel::<Option<SequenceGroup>>(num_workers * 2);
    let (output_tx, output_rx) = mpsc::sync_channel::<Option<OutputRecord>>(num_workers * 2);

    let total_count = Arc::new(AtomicUsize::new(0));
    let verbose = args.verbose;
    let max_length = args.max_length;

    // Spawn reader thread
    let reader_total = Arc::clone(&total_count);
    let input_path = args.input.clone();
    let reader_handle = thread::spawn(move || {
        reader_thread(input_path, input_tx, num_workers, reader_total);
    });

    // Spawn worker threads
    let input_rx = Arc::new(std::sync::Mutex::new(input_rx));
    let mut worker_handles = Vec::new();

    for _ in 0..num_workers {
        let rx = Arc::clone(&input_rx);
        let tx = output_tx.clone();

        let handle = thread::spawn(move || {
            loop {
                let group = {
                    let rx = rx.lock().unwrap();
                    rx.recv()
                };

                match group {
                    Ok(Some(group)) => {
                        let result = process_group(group, max_length);
                        tx.send(Some(result)).expect("Failed to send result");
                    }
                    Ok(None) | Err(_) => {
                        tx.send(None).expect("Failed to send termination");
                        break;
                    }
                }
            }
        });
        worker_handles.push(handle);
    }

    // Drop extra sender so writer knows when to stop
    drop(output_tx);

    // Spawn writer thread
    let writer_total = Arc::clone(&total_count);
    let output_path = args.output.clone();
    let writer_handle = thread::spawn(move || {
        writer_thread(output_rx, output_path, writer_total, num_workers, verbose);
    });

    // Wait for completion
    reader_handle.join().expect("Reader thread panicked");
    for handle in worker_handles {
        handle.join().expect("Worker thread panicked");
    }
    writer_handle.join().expect("Writer thread panicked");

    eprintln!("Done!");
}
