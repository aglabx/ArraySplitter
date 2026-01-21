//! ArraySplitter CLI - De novo decomposition of satellite DNA arrays into monomers
//!
//! Architecture: Reader -> Workers -> Writer (bounded channels, constant memory)

use clap::Parser;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::sync::mpsc::{self, SyncSender, Receiver};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::thread;

use arraysplitter_rs::{
    decompose_array, decompose_array_with_cuts, apply_all_heuristics,
    parse_cuts, get_revcomp,
    decompose::is_canonical_orientation,
};

/// De novo decomposition of satellite DNA arrays into monomers
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input FASTA file
    #[arg(short, long)]
    input: String,

    /// Output prefix
    #[arg(short, long)]
    output: String,

    /// Number of worker threads [default: number of CPUs]
    #[arg(short, long)]
    threads: Option<usize>,

    /// Predefined cut sequences (comma-separated)
    #[arg(short, long)]
    cuts: Option<String>,

    /// Depth for hint discovery [default: 100]
    #[arg(short, long, default_value = "100")]
    depth: usize,

    /// Verbose output
    #[arg(short, long)]
    verbose: bool,
}

/// Input record for workers
struct InputRecord {
    id: String,
    sequence: String,
}

/// Output record from workers
struct OutputRecord {
    header: String,
    decomposition: Vec<String>,
    cut_sequence: String,
    was_reversed: bool,
    period: usize,
}

/// Process a single array
fn process_array(
    id: &str,
    array: &str,
    predefined_cuts: &Option<Vec<String>>,
    depth: usize,
    verbose: bool,
) -> OutputRecord {
    let is_canonical = is_canonical_orientation(array);
    let (working_array, was_reversed) = if is_canonical {
        (array.to_string(), false)
    } else {
        (get_revcomp(array), true)
    };

    // Adaptive depth for long sequences
    let effective_depth = if working_array.len() > 100_000 {
        depth.min(30)
    } else if working_array.len() > 50_000 {
        depth.min(50)
    } else {
        depth
    };

    let result = if let Some(cuts) = predefined_cuts {
        decompose_array_with_cuts(&working_array, cuts, verbose)
    } else {
        decompose_array(&working_array, effective_depth, None, verbose)
    };

    let decomposition = apply_all_heuristics(&result.monomers, &result.cut_sequence, verbose);

    OutputRecord {
        header: id.to_string(),
        decomposition,
        cut_sequence: result.cut_sequence,
        was_reversed,
        period: result.period,
    }
}

/// Reader thread: reads FASTA and sends records to channel
fn reader_thread(
    input_path: String,
    tx: SyncSender<Option<InputRecord>>,
    num_workers: usize,
    total_count: Arc<AtomicUsize>,
) {
    let file = File::open(&input_path).expect("Failed to open input file");
    let reader = BufReader::new(file);

    let mut current_id = String::new();
    let mut current_seq = String::new();
    let mut count = 0;

    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        let trimmed = line.trim();

        if trimmed.starts_with('>') {
            // Send previous record if exists
            if !current_id.is_empty() && !current_seq.is_empty() {
                tx.send(Some(InputRecord {
                    id: current_id.clone(),
                    sequence: current_seq.clone(),
                })).expect("Failed to send record");
                count += 1;
            }

            // Parse new header
            let header = &trimmed[1..];
            let parts: Vec<&str> = header.splitn(2, |c| c == ' ' || c == '\t').collect();
            current_id = parts[0].to_string();
            current_seq = String::new();
        } else if !trimmed.is_empty() {
            current_seq.push_str(trimmed);
        }
    }

    // Send last record
    if !current_id.is_empty() && !current_seq.is_empty() {
        tx.send(Some(InputRecord {
            id: current_id,
            sequence: current_seq,
        })).expect("Failed to send last record");
        count += 1;
    }

    total_count.store(count, Ordering::SeqCst);

    // Send termination signals for all workers
    for _ in 0..num_workers {
        tx.send(None).expect("Failed to send termination");
    }
}

/// Writer thread: receives results and writes to files
fn writer_thread(
    rx: Receiver<Option<OutputRecord>>,
    output_prefix: String,
    total_count: Arc<AtomicUsize>,
    num_workers: usize,
) {
    let output_file = format!("{}.decomposed.fasta", output_prefix);
    let detail_file = format!("{}.monomers.tsv", output_prefix);
    let lengths_file = format!("{}.lengths", output_prefix);

    let mut fw = BufWriter::new(File::create(&output_file).expect("Failed to create output file"));
    let mut fw_detail = BufWriter::new(File::create(&detail_file).expect("Failed to create detail file"));
    let mut fw_lengths = BufWriter::new(File::create(&lengths_file).expect("Failed to create lengths file"));

    writeln!(fw_detail, "sequence_id\torientation\tindex\ttype\tlength\tis_flank\tsequence").unwrap();

    let mut processed = 0;
    let mut finished_workers = 0;
    let mut total_monomers = 0;

    loop {
        match rx.recv() {
            Ok(Some(result)) => {
                processed += 1;
                total_monomers += result.decomposition.len();

                let decomposition = &result.decomposition;
                let cut_seq = &result.cut_sequence;
                let orientation = if result.was_reversed { "rev" } else { "fwd" };

                // Calculate flank threshold
                let all_monomer_lengths: Vec<usize> = decomposition
                    .iter()
                    .filter(|m| m.starts_with(cut_seq))
                    .map(|m| m.len())
                    .collect();

                let flank_threshold = if !all_monomer_lengths.is_empty() {
                    let avg: f64 = all_monomer_lengths.iter().sum::<usize>() as f64
                        / all_monomer_lengths.len() as f64;
                    (avg * 0.7) as usize
                } else {
                    result.period * 50 / 100
                };

                // Count internal monomers
                let mut internal_lengths: Vec<usize> = Vec::new();
                for (i, m) in decomposition.iter().enumerate() {
                    if m.starts_with(cut_seq) {
                        if i == decomposition.len() - 1 && m.len() < flank_threshold {
                            continue;
                        }
                        internal_lengths.push(m.len());
                    }
                }

                let header_info = if !internal_lengths.is_empty() {
                    let min_len = *internal_lengths.iter().min().unwrap();
                    let max_len = *internal_lengths.iter().max().unwrap();
                    let avg_len: f64 = internal_lengths.iter().sum::<usize>() as f64
                        / internal_lengths.len() as f64;
                    format!(
                        "{} cut={} orientation={} n_monomers={} range={}-{} avg={:.1}",
                        result.header, cut_seq, orientation, internal_lengths.len(),
                        min_len, max_len, avg_len
                    )
                } else {
                    format!("{} cut={} orientation={} n_monomers=0",
                        result.header, cut_seq, orientation)
                };

                // Write FASTA
                writeln!(fw, ">{}", header_info).unwrap();
                writeln!(fw, "{}", decomposition.join(" ")).unwrap();

                // Write lengths
                writeln!(fw_lengths, ">{}", header_info).unwrap();
                let lengths: Vec<String> = decomposition.iter().map(|m| m.len().to_string()).collect();
                writeln!(fw_lengths, "{}", lengths.join(" ")).unwrap();

                // Write detail TSV
                for (i, monomer) in decomposition.iter().enumerate() {
                    let (piece_type, is_flank) = if i == 0 && !monomer.starts_with(cut_seq) {
                        ("LEFT_FLANK", "TRUE")
                    } else if i == decomposition.len() - 1 && monomer.len() < flank_threshold {
                        ("RIGHT_FLANK", "TRUE")
                    } else {
                        ("MONOMER", "FALSE")
                    };

                    writeln!(
                        fw_detail, "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        result.header, orientation, i, piece_type,
                        monomer.len(), is_flank, monomer
                    ).unwrap();
                }

                // Progress output
                let total = total_count.load(Ordering::SeqCst);
                if total > 0 && processed % 100 == 0 {
                    eprint!("\rProcessed: {}/{} ({:.1}%)",
                        processed, total, (processed as f64 / total as f64) * 100.0);
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

    // Setup output prefix
    let mut output_prefix = args.output.clone();
    if output_prefix.ends_with(".fasta") {
        output_prefix = output_prefix[..output_prefix.len() - 6].to_string();
    } else if output_prefix.ends_with(".fa") {
        output_prefix = output_prefix[..output_prefix.len() - 3].to_string();
    }

    eprintln!("Input: {}", args.input);
    eprintln!("Output prefix: {}", output_prefix);

    let predefined_cuts: Option<Vec<String>> = args.cuts.as_ref().map(|s| parse_cuts(s));
    if let Some(ref cuts) = predefined_cuts {
        eprintln!("Using predefined cuts: {:?}", cuts);
    }

    // Channels with bounded capacity (controls memory)
    let (input_tx, input_rx) = mpsc::sync_channel::<Option<InputRecord>>(num_workers * 2);
    let (output_tx, output_rx) = mpsc::sync_channel::<Option<OutputRecord>>(num_workers * 2);

    let total_count = Arc::new(AtomicUsize::new(0));
    let depth = args.depth;
    let verbose = args.verbose;

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
        let cuts = predefined_cuts.clone();

        let handle = thread::spawn(move || {
            loop {
                let record = {
                    let rx = rx.lock().unwrap();
                    rx.recv()
                };

                match record {
                    Ok(Some(input)) => {
                        let result = process_array(
                            &input.id,
                            &input.sequence,
                            &cuts,
                            depth,
                            verbose,
                        );
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
    let writer_handle = thread::spawn(move || {
        writer_thread(output_rx, output_prefix, writer_total, num_workers);
    });

    // Wait for completion
    reader_handle.join().expect("Reader thread panicked");
    for handle in worker_handles {
        handle.join().expect("Worker thread panicked");
    }
    writer_handle.join().expect("Writer thread panicked");

    eprintln!("Done!");
}
