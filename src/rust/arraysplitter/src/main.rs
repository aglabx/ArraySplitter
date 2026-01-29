//! ArraySplitter CLI - De novo decomposition of satellite DNA arrays into monomers
//!
//! Architecture: Reader -> Workers -> Writer (bounded channels, constant memory)

use clap::Parser;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::sync::mpsc::{self, SyncSender, Receiver};
use std::sync::atomic::{AtomicUsize, AtomicU64, Ordering};
use std::sync::Arc;
use std::thread;
use std::time::Instant;
use sysinfo::{System, Pid, ProcessRefreshKind};

use arraysplitter_rs::{
    decompose_array, decompose_array_with_cuts, decompose_array_autocorr,
    apply_all_heuristics, parse_cuts, get_revcomp,
    decompose::is_canonical_orientation,
    recursive_hor::{RecursiveResult, decompose_hors_to_base, DEFAULT_MIN_SUBMONOMER_LEN, DEFAULT_AUTOCORR_THRESHOLD},
};

/// Global maximum monomer length for ED calculation (set from CLI)
static MAX_ED_LEN: AtomicUsize = AtomicUsize::new(10000);

/// Calculate edit distance between two sequences
/// Returns usize::MAX/2 if either string is too long to avoid quadratic blowup
fn edit_distance(s1: &str, s2: &str) -> usize {
    let len1 = s1.len();
    let len2 = s2.len();
    let max_ed_len = MAX_ED_LEN.load(Ordering::Relaxed);

    // Skip ED for very long monomers
    if len1 > max_ed_len || len2 > max_ed_len {
        return usize::MAX / 2;
    }

    if len1 == 0 { return len2; }
    if len2 == 0 { return len1; }

    let s1_chars: Vec<char> = s1.chars().collect();
    let s2_chars: Vec<char> = s2.chars().collect();

    let mut prev_row: Vec<usize> = (0..=len2).collect();
    let mut curr_row: Vec<usize> = vec![0; len2 + 1];

    for i in 1..=len1 {
        curr_row[0] = i;
        for j in 1..=len2 {
            let cost = if s1_chars[i - 1] == s2_chars[j - 1] { 0 } else { 1 };
            curr_row[j] = (prev_row[j] + 1)
                .min(curr_row[j - 1] + 1)
                .min(prev_row[j - 1] + cost);
        }
        std::mem::swap(&mut prev_row, &mut curr_row);
    }

    prev_row[len2]
}

/// Compute consensus sequence and Phred-like quality string from monomers.
///
/// Quality encodes the fraction of monomers supporting the consensus base:
/// Q = -10 * log10(1 - support_fraction), capped at 93.
/// Encoded as ASCII char = Q + 33 (same as Illumina Phred+33).
///
/// Examples: 100% support → '~' (Q=93), 90% → '+' (Q=10), 50% → '$' (Q=3)
/// Compute consensus from aligned monomers.
/// Returns (consensus_with_case, iupac, quality_digits):
/// - consensus_with_case: majority-vote base; UPPERCASE if support ≥ 80%, lowercase if variable
/// - iupac: IUPAC ambiguity code (bases with frequency ≥ 20% included)
/// - quality_digits: digit 0-9 per position representing support fraction (9 = 100%, 5 = 50-59%)
fn compute_consensus(monomers: &[&str]) -> (String, String, String) {
    if monomers.is_empty() {
        return (String::new(), String::new(), String::new());
    }

    // Use median length as consensus length
    let mut lengths: Vec<usize> = monomers.iter().map(|m| m.len()).collect();
    lengths.sort();
    let consensus_len = lengths[lengths.len() / 2];

    let mut consensus = String::with_capacity(consensus_len);
    let mut iupac = String::with_capacity(consensus_len);
    let mut quality = String::with_capacity(consensus_len);

    for pos in 0..consensus_len {
        let mut counts = [0u32; 4]; // A, C, G, T
        let mut total = 0u32;

        for monomer in monomers {
            let bytes = monomer.as_bytes();
            if pos < bytes.len() {
                match bytes[pos] {
                    b'A' | b'a' => { counts[0] += 1; total += 1; },
                    b'C' | b'c' => { counts[1] += 1; total += 1; },
                    b'G' | b'g' => { counts[2] += 1; total += 1; },
                    b'T' | b't' => { counts[3] += 1; total += 1; },
                    _ => {},
                }
            }
        }

        if total == 0 {
            consensus.push('N');
            iupac.push('N');
            quality.push('0');
            continue;
        }

        // Find most common base for consensus
        let (best_idx, &best_count) = counts.iter().enumerate()
            .max_by_key(|(_, &c)| c)
            .unwrap();

        let support = best_count as f64 / total as f64;

        // Consensus: always uppercase (quality digits show conservation)
        let base_upper = match best_idx {
            0 => 'A', 1 => 'C', 2 => 'G', 3 => 'T', _ => 'N',
        };
        consensus.push(base_upper);

        // Quality: digit 0-9 representing support fraction
        // 9 = 90-100%, 8 = 80-89%, ..., 0 = 0-9%
        let digit = (support * 10.0).min(9.99) as u8;
        quality.push((b'0' + digit) as char);

        // IUPAC: include bases with frequency ≥ 20%
        let threshold = (total as f64 * 0.2) as u32;
        let present: Vec<usize> = (0..4)
            .filter(|&i| counts[i] >= threshold.max(1))
            .collect();

        let iupac_char = match present.as_slice() {
            [0] => 'A', [1] => 'C', [2] => 'G', [3] => 'T',
            [0, 2] => 'R', // A+G purine
            [1, 3] => 'Y', // C+T pyrimidine
            [1, 2] => 'S', // G+C strong
            [0, 3] => 'W', // A+T weak
            [2, 3] => 'K', // G+T keto
            [0, 1] => 'M', // A+C amino
            [1, 2, 3] => 'B', // not A
            [0, 2, 3] => 'D', // not C
            [0, 1, 3] => 'H', // not G
            [0, 1, 2] => 'V', // not T
            [0, 1, 2, 3] => 'N',
            _ => 'N',
        };
        iupac.push(iupac_char);
    }

    (consensus, iupac, quality)
}

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

    /// Show statistics and outliers after processing
    #[arg(short, long)]
    stats: bool,

    /// Number of outliers to show [default: 20]
    #[arg(long, default_value = "20")]
    top_outliers: usize,

    /// Maximum monomer length for ED calculation [default: 10000]
    #[arg(long, default_value = "10000")]
    max_ed_len: usize,

    /// Decomposition method: autocorr (new pipeline), classic (old), both (run both for comparison)
    #[arg(long, default_value = "autocorr")]
    method: String,
}

/// Input record for workers
struct InputRecord {
    id: String,
    sequence: String,
}

/// Pre-computed per-monomer info
struct MonomerOutput {
    sequence: String,
    piece_type: String,  // "monomer" or "flank"
    source: String,
    ed_tmpl: Option<usize>,
    ed_prev: Option<usize>,
    ed_next: Option<usize>,
}

/// Output record from workers (all computation done)
struct OutputRecord {
    header: String,
    monomers: Vec<MonomerOutput>,
    cut_sequence: String,
    was_reversed: bool,
    period: usize,
    array_len: usize,
    autocorr_period: Option<usize>,
    autocorr_value: Option<f64>,
    autocorr_excess: Option<f64>,
    method: String,
    // Pre-computed statistics:
    n_monomers: usize,
    ed_per_bp: f64,
    cv: f64,
    mean_ed_tmpl: f64,
    mean_ed_prev: f64,
    mean_ed_next: f64,
    // Pre-computed consensus:
    consensus_seq: String,
    iupac_str: String,
    quality_str: String,
    // Recursive HOR decomposition:
    recursive_result: RecursiveResult,
}

/// Process a single array
fn process_array(
    id: &str,
    array: &str,
    predefined_cuts: &Option<Vec<String>>,
    depth: usize,
    verbose: bool,
    method: &str,
) -> OutputRecord {
    let array_len = array.len();

    // Normalize to uppercase for consistent processing
    let array_upper = array.to_uppercase();

    let is_canonical = is_canonical_orientation(&array_upper);
    let (working_array, was_reversed) = if is_canonical {
        (array_upper, false)
    } else {
        (get_revcomp(&array_upper), true)
    };

    let (decomposition, cut_sequence, period, autocorr_period, autocorr_value, monomer_sources) =
        match method {
            "autocorr" => {
                let result = decompose_array_autocorr(&working_array, verbose);
                // No heuristics for autocorr method (it handles its own splitting)
                (
                    result.monomers,
                    result.cut_sequence,
                    result.period,
                    result.autocorr_period,
                    result.autocorr_value,
                    result.monomer_sources,
                )
            }
            _ => {
                // Classic method
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

                let decomposition = apply_all_heuristics(
                    &result.monomers,
                    &result.cut_sequence,
                    &result.alternative_anchors,
                    verbose,
                );
                (
                    decomposition,
                    result.cut_sequence,
                    result.period,
                    None,
                    None,
                    Vec::new(),
                )
            }
        };

    // --- Compute all per-monomer stats in the worker thread ---

    // Determine flank threshold
    let flank_threshold = if period > 0 {
        (period as f64 * 0.7) as usize
    } else if !decomposition.is_empty() {
        let avg: f64 = decomposition.iter().map(|m| m.len()).sum::<usize>() as f64
            / decomposition.len() as f64;
        (avg * 0.7) as usize
    } else {
        0
    };

    // Classify monomer vs flank
    let piece_types: Vec<&str> = decomposition.iter().enumerate().map(|(i, monomer)| {
        let source = if i < monomer_sources.len() { monomer_sources[i].as_str() } else { "" };
        if source.contains("flank") {
            "flank"
        } else if source == "no_period" || source == "no_anchor" {
            "monomer"
        } else if method == "autocorr" {
            if monomer.len() < flank_threshold && (i == 0 || i == decomposition.len() - 1) {
                "flank"
            } else {
                "monomer"
            }
        } else {
            // Classic method
            if monomer.len() < flank_threshold && (i == 0 || i == decomposition.len() - 1) {
                "flank"
            } else {
                "monomer"
            }
        }
    }).collect();

    // Select template monomer (median length among non-flanks)
    let mut non_flank_lengths: Vec<(usize, usize)> = decomposition.iter().enumerate()
        .filter(|(i, _)| piece_types[*i] == "monomer")
        .map(|(i, m)| (i, m.len()))
        .collect();
    non_flank_lengths.sort_by_key(|(_, len)| *len);

    let template_idx = if !non_flank_lengths.is_empty() {
        Some(non_flank_lengths[non_flank_lengths.len() / 2].0)
    } else {
        None
    };
    let template_monomer = template_idx.map(|idx| decomposition[idx].as_str());

    // Compute EDs
    let mut ed_tmpls: Vec<Option<usize>> = vec![None; decomposition.len()];
    let mut ed_prevs: Vec<Option<usize>> = vec![None; decomposition.len()];
    let mut ed_nexts: Vec<Option<usize>> = vec![None; decomposition.len()];

    for (i, monomer) in decomposition.iter().enumerate() {
        if piece_types[i] != "monomer" { continue; }

        if let Some(tmpl) = template_monomer {
            if template_idx == Some(i) {
                ed_tmpls[i] = Some(0);
            } else {
                ed_tmpls[i] = Some(edit_distance(tmpl, monomer));
            }
        }
        if i > 0 && piece_types[i - 1] == "monomer" {
            ed_prevs[i] = Some(edit_distance(&decomposition[i - 1], monomer));
        }
        if i + 1 < decomposition.len() && piece_types[i + 1] == "monomer" {
            ed_nexts[i] = Some(edit_distance(monomer, &decomposition[i + 1]));
        }
    }

    // Compute array-level statistics
    let n_monomers = piece_types.iter().filter(|&&t| t == "monomer").count();
    let ed_tmpl_values: Vec<f64> = ed_tmpls.iter().filter_map(|e| e.map(|v| v as f64)).collect();
    let ed_prev_values: Vec<f64> = ed_prevs.iter().filter_map(|e| e.map(|v| v as f64)).collect();
    let ed_next_values: Vec<f64> = ed_nexts.iter().filter_map(|e| e.map(|v| v as f64)).collect();
    let monomer_len_values: Vec<f64> = decomposition.iter().enumerate()
        .filter(|(i, _)| piece_types[*i] == "monomer")
        .map(|(_, m)| m.len() as f64)
        .collect();

    let mean_ed_tmpl = if !ed_tmpl_values.is_empty() {
        ed_tmpl_values.iter().sum::<f64>() / ed_tmpl_values.len() as f64
    } else { 0.0 };
    let mean_ed_prev = if !ed_prev_values.is_empty() {
        ed_prev_values.iter().sum::<f64>() / ed_prev_values.len() as f64
    } else { 0.0 };
    let mean_ed_next = if !ed_next_values.is_empty() {
        ed_next_values.iter().sum::<f64>() / ed_next_values.len() as f64
    } else { 0.0 };
    let mean_monomer_len = if !monomer_len_values.is_empty() {
        monomer_len_values.iter().sum::<f64>() / monomer_len_values.len() as f64
    } else { 0.0 };
    let ed_per_bp = if mean_monomer_len > 0.0 { mean_ed_prev / mean_monomer_len } else { 0.0 };

    let cv = if monomer_len_values.len() > 1 && mean_monomer_len > 0.0 {
        let variance: f64 = monomer_len_values.iter()
            .map(|&x| (x - mean_monomer_len).powi(2))
            .sum::<f64>() / monomer_len_values.len() as f64;
        variance.sqrt() / mean_monomer_len
    } else { 0.0 };

    // Compute consensus (all monomers except flanks)
    let consensus_monomers: Vec<&str> = decomposition.iter().enumerate()
        .filter(|(i, _)| piece_types[*i] == "monomer")
        .map(|(_, m)| m.as_str())
        .collect();

    let (consensus_seq, iupac_str, quality_str) = if consensus_monomers.len() >= 2 {
        compute_consensus(&consensus_monomers)
    } else {
        (String::new(), String::new(), String::new())
    };

    // Recursive HOR decomposition: decompose each HOR monomer into base-level monomers
    let hor_monomer_seqs: Vec<String> = decomposition.iter()
        .enumerate()
        .filter(|(i, _)| piece_types[*i] == "monomer")
        .map(|(_, m)| m.clone())
        .collect();

    let recursive_result = if !hor_monomer_seqs.is_empty() {
        decompose_hors_to_base(
            &hor_monomer_seqs,
            DEFAULT_MIN_SUBMONOMER_LEN,
            DEFAULT_AUTOCORR_THRESHOLD,
        )
    } else {
        RecursiveResult {
            base_monomers: Vec::new(),
            max_depth: 0,
            n_expected: 0,
            median_period: 0,
            mean_autocorr: 0.0,
            consensus_seq: String::new(),
            iupac_str: String::new(),
            quality_str: String::new(),
        }
    };

    // Build MonomerOutput vec
    let monomers: Vec<MonomerOutput> = decomposition.into_iter().enumerate().map(|(i, seq)| {
        let source = if i < monomer_sources.len() {
            monomer_sources[i].clone()
        } else if method == "classic" {
            if piece_types[i] == "flank" { "classic".to_string() } else { "anchor_graph".to_string() }
        } else {
            "-".to_string()
        };
        MonomerOutput {
            sequence: seq,
            piece_type: piece_types[i].to_string(),
            source,
            ed_tmpl: ed_tmpls[i],
            ed_prev: ed_prevs[i],
            ed_next: ed_nexts[i],
        }
    }).collect();

    OutputRecord {
        header: id.to_string(),
        monomers,
        cut_sequence,
        was_reversed,
        period,
        array_len,
        autocorr_period,
        autocorr_value,
        autocorr_excess: autocorr_value.map(|v| v - 0.25), // excess over uniform random
        method: method.to_string(),
        n_monomers,
        ed_per_bp,
        cv,
        mean_ed_tmpl,
        mean_ed_prev,
        mean_ed_next,
        consensus_seq,
        iupac_str,
        quality_str,
        recursive_result,
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

/// Per-array statistics for --stats output
struct ArrayStats {
    id: String,
    n_monomers: usize,
    lengths: Vec<usize>,
    eds: Vec<usize>,
}

/// Calculate robust CV (MAD/median)
fn robust_cv(vals: &[usize]) -> f64 {
    if vals.len() < 3 {
        return f64::INFINITY;
    }
    let mut sorted: Vec<usize> = vals.to_vec();
    sorted.sort();
    let median = sorted[sorted.len() / 2];
    if median == 0 {
        return f64::INFINITY;
    }
    let mut deviations: Vec<usize> = sorted.iter()
        .map(|&v| if v > median { v - median } else { median - v })
        .collect();
    deviations.sort();
    let mad = deviations[deviations.len() / 2];
    mad as f64 / median as f64
}

/// Sort TSV file by chromosome and position (column 1 only)
/// Header is preserved at top
fn sort_tsv_file(path: &str) {
    use std::process::Command;

    // Use external sort: keep header, then sort rest by genomic position (version sort)
    // Note: using bash for $'\t' syntax support
    let script = format!(
        r#"head -1 "{path}" > "{path}.sorted" && \
           tail -n +2 "{path}" | sort -t$'\t' -k1,1V >> "{path}.sorted" && \
           /bin/mv -f "{path}.sorted" "{path}""#,
        path = path
    );

    let result = Command::new("bash")
        .arg("-c")
        .arg(&script)
        .output();

    match result {
        Ok(output) => {
            if !output.status.success() {
                eprintln!("Warning: Failed to sort {}: {}",
                    path, String::from_utf8_lossy(&output.stderr));
            }
        }
        Err(e) => {
            eprintln!("Warning: Failed to run sort on {}: {}", path, e);
        }
    }
}

/// Sort TSV file with type priority (for hors/monomers files)
/// Order: pred_array -> flank -> monomer/base_monomer (by idx) -> array -> consensus
/// Sorted by: 1) array_id (genomic order), 2) type priority, 3) idx
fn sort_tsv_with_type_priority(path: &str) {
    use std::process::Command;

    // AWK script to add type priority column, sort, then remove it
    // Type priorities: pred_array=0, flank=1, monomer=2, base_monomer=2, array=3, consensus=4
    // Note: using bash for $'\t' syntax support
    let script = format!(
        r#"head -1 "{path}" > "{path}.sorted" && \
           tail -n +2 "{path}" | awk -F'\t' 'BEGIN{{OFS="\t"}} {{
               if ($2=="pred_array") p=0;
               else if ($2=="flank") p=1;
               else if ($2=="monomer" || $2=="base_monomer") p=2;
               else if ($2=="array") p=3;
               else if ($2=="consensus") p=4;
               else p=5;
               print p, $0
           }}' | sort -t$'\t' -k2,2V -k1,1n -k4,4n | cut -f2- >> "{path}.sorted" && \
           /bin/mv -f "{path}.sorted" "{path}""#,
        path = path
    );

    let result = Command::new("bash")
        .arg("-c")
        .arg(&script)
        .output();

    match result {
        Ok(output) => {
            if !output.status.success() {
                eprintln!("Warning: Failed to sort {}: {}",
                    path, String::from_utf8_lossy(&output.stderr));
            }
        }
        Err(e) => {
            eprintln!("Warning: Failed to run sort on {}: {}", path, e);
        }
    }
}

/// Writer thread: receives results and writes to files (streaming, no memory accumulation)
fn writer_thread(
    rx: Receiver<Option<OutputRecord>>,
    output_prefix: String,
    total_count: Arc<AtomicUsize>,
    num_workers: usize,
    collect_stats: bool,
    top_outliers: usize,
) {
    let output_file = format!("{}.decomposed.fasta", output_prefix);
    let hors_file = format!("{}.hors.tsv", output_prefix);
    let monomers_file = format!("{}.monomers.tsv", output_prefix);
    let lengths_file = format!("{}.lengths", output_prefix);
    let summary_file = format!("{}.summary.tsv", output_prefix);

    let mut fw = BufWriter::new(File::create(&output_file).expect("Failed to create output file"));
    let mut fw_hors = BufWriter::new(File::create(&hors_file).expect("Failed to create HORs file"));
    let mut fw_monomers = BufWriter::new(File::create(&monomers_file).expect("Failed to create monomers file"));
    let mut fw_lengths = BufWriter::new(File::create(&lengths_file).expect("Failed to create lengths file"));
    let mut fw_summary = BufWriter::new(File::create(&summary_file).expect("Failed to create summary file"));

    // HORs TSV header (HOR-level decomposition)
    writeln!(fw_hors, "array_id\ttype\tidx\tlength\tsource\ted_tmpl\ted_prev\ted_next\tperiod\tautocorr\tn_expected\ted_per_bp\tcv\tcut_sequence\torientation\tsequence").unwrap();

    // Monomers TSV header (base-level monomers - unified format with parent_idx)
    writeln!(fw_monomers, "array_id\ttype\tidx\tlength\tsource\ted_tmpl\ted_prev\ted_next\tperiod\tautocorr\tn_expected\ted_per_bp\tcv\tcut_sequence\torientation\tparent_idx\tsequence").unwrap();

    // Summary TSV header (one row per array with HOR and monomer statistics + consensus)
    writeln!(fw_summary, "array_id\tarray_length\torientation\tmethod\thor_period\thor_autocorr\thor_n_monomers\thor_mean_ed_tmpl\thor_mean_ed_prev\thor_cv\thor_consensus\thor_iupac\thor_quality\tmono_period\tmono_autocorr\tmono_n_monomers\tmono_mean_ed_tmpl\tmono_mean_ed_prev\tmono_cv\tmono_consensus\tmono_iupac\tmono_quality\tcut_sequence").unwrap();

    let mut processed = 0;
    let mut finished_workers = 0;
    let mut total_monomers = 0;

    // Collect per-array stats if requested (only stats, not full results)
    let mut all_stats: Vec<ArrayStats> = Vec::new();

    // Streaming: write each result immediately as it arrives
    loop {
        match rx.recv() {
            Ok(Some(result)) => {
                processed += 1;

                // Progress output
                let total = total_count.load(Ordering::SeqCst);
                if total > 0 {
                    eprint!("\rProcessing: {}/{} ({:.1}%)",
                        processed, total, (processed as f64 / total as f64) * 100.0);
                }

                // Write this result immediately
                total_monomers += result.monomers.len();

                let cut_seq = &result.cut_sequence;
                let orientation = if result.was_reversed { "rev" } else { "fwd" };

                let autocorr_str = result.autocorr_value
                    .map(|v| format!("{:.4}", v))
                    .unwrap_or_else(|| "-".to_string());
                let n_expected_str = if result.period > 0 {
                    format!("{}", result.array_len / result.period)
                } else {
                    "-".to_string()
                };

                // Write pred_array row
                let excess_str = result.autocorr_excess
                    .map(|v| format!("{:.4}", v))
                    .unwrap_or_else(|| "0".to_string());
                writeln!(
                    fw_hors, "{}\tpred_array\t{}\t{}\t{}\t{}\t0\t0\t{}\t{}\t{}\t0\t0\t-\t{}\t-",
                    result.header,
                    n_expected_str,
                    result.array_len,
                    result.method,
                    excess_str,
                    result.autocorr_period.map(|p| p.to_string()).unwrap_or_else(|| "0".to_string()),
                    autocorr_str,
                    n_expected_str,
                    orientation,
                ).unwrap();

                // Write individual monomer/flank rows
                for (i, mono) in result.monomers.iter().enumerate() {
                    let ed_tmpl_str = mono.ed_tmpl.map(|v| v.to_string()).unwrap_or_else(|| "0".to_string());
                    let ed_prev_str = mono.ed_prev.map(|v| v.to_string()).unwrap_or_else(|| "0".to_string());
                    let ed_next_str = mono.ed_next.map(|v| v.to_string()).unwrap_or_else(|| "0".to_string());

                    let mono_len = mono.sequence.len() as f64;
                    let mono_ed_per_bp = mono.ed_prev
                        .map(|ed| if mono_len > 0.0 { ed as f64 / mono_len } else { 0.0 })
                        .unwrap_or(0.0);
                    let size_dev = if result.period > 0 {
                        (mono_len - result.period as f64).abs() / result.period as f64
                    } else { 0.0 };

                    writeln!(
                        fw_hors, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t1\t{:.4}\t{:.4}\t{}\t{}\t{}",
                        result.header, mono.piece_type, i, mono.sequence.len(),
                        mono.source, ed_tmpl_str, ed_prev_str, ed_next_str,
                        result.period,
                        mono.ed_tmpl.map(|ed| ed as f64 / mono_len).unwrap_or(0.0),
                        mono_ed_per_bp,
                        size_dev,
                        cut_seq,
                        orientation,
                        mono.sequence
                    ).unwrap();
                }

                // Write array summary row
                writeln!(
                    fw_hors, "{}\tarray\t{}\t{}\t{}\t{:.1}\t{:.1}\t{:.1}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{}\t{}\t-",
                    result.header,
                    result.n_monomers,
                    result.array_len,
                    result.method,
                    result.mean_ed_tmpl,
                    result.mean_ed_prev,
                    result.mean_ed_next,
                    result.period,
                    autocorr_str,
                    result.n_monomers,
                    result.ed_per_bp,
                    result.cv,
                    cut_seq,
                    orientation,
                ).unwrap();

                // Write consensus row
                if !result.consensus_seq.is_empty() {
                    let conserved = result.quality_str.chars().filter(|c| *c >= '8').count();
                    let variable = result.quality_str.chars().filter(|c| *c < '8').count();
                    let ambiguous = result.iupac_str.chars()
                        .filter(|c| !matches!(c, 'A' | 'C' | 'G' | 'T'))
                        .count();
                    writeln!(
                        fw_hors, "{}\tconsensus\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{}\t{}\t{}",
                        result.header,
                        result.n_monomers,
                        result.consensus_seq.len(),
                        result.quality_str,
                        conserved,
                        variable,
                        ambiguous,
                        result.period,
                        autocorr_str,
                        result.n_monomers,
                        result.ed_per_bp,
                        result.cv,
                        result.iupac_str,
                        orientation,
                        result.consensus_seq,
                    ).unwrap();
                }

                // Write pred_array line for monomers file (submonomer-level statistics)
                let rec = &result.recursive_result;
                if !rec.base_monomers.is_empty() {
                    // Compute mean ed_per_bp and cv for submonomers
                    let ed_per_bp_values: Vec<f64> = rec.base_monomers.iter()
                        .map(|m| m.ed_per_bp)
                        .collect();
                    let mean_ed_per_bp = if !ed_per_bp_values.is_empty() {
                        ed_per_bp_values.iter().sum::<f64>() / ed_per_bp_values.len() as f64
                    } else { 0.0 };

                    let cv_values: Vec<f64> = rec.base_monomers.iter()
                        .map(|m| m.cv)
                        .collect();
                    let mean_cv = if !cv_values.is_empty() {
                        cv_values.iter().sum::<f64>() / cv_values.len() as f64
                    } else { 0.0 };

                    writeln!(
                        fw_monomers, "{}\tpred_array\t{}\t{}\trecursive\t0\t0\t0\t{}\t{:.4}\t{}\t{:.4}\t{:.4}\t-\t{}\t-\t-",
                        result.header,
                        rec.n_expected,
                        result.array_len,
                        rec.median_period,
                        rec.mean_autocorr,
                        rec.n_expected,
                        mean_ed_per_bp,
                        mean_cv,
                        orientation,
                    ).unwrap();
                }

                // Write base-level monomers from recursive HOR decomposition (unified format)
                for base_mono in &rec.base_monomers {
                    let ed_tmpl_str = base_mono.ed_tmpl.map(|v| v.to_string()).unwrap_or_else(|| "0".to_string());
                    let ed_prev_str = base_mono.ed_prev.map(|v| v.to_string()).unwrap_or_else(|| "0".to_string());
                    let ed_next_str = base_mono.ed_next.map(|v| v.to_string()).unwrap_or_else(|| "0".to_string());

                    // Determine type: "base_monomer" if has parent, "monomer" if single (no further decomp)
                    let mono_type = if rec.base_monomers.len() == 1 && base_mono.source == "base" {
                        "monomer"
                    } else {
                        "base_monomer"
                    };

                    writeln!(
                        fw_monomers, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t1\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}",
                        result.header,
                        mono_type,
                        base_mono.global_idx,
                        base_mono.sequence.len(),
                        base_mono.source,
                        ed_tmpl_str,
                        ed_prev_str,
                        ed_next_str,
                        base_mono.period,
                        base_mono.autocorr,
                        base_mono.ed_per_bp,
                        base_mono.cv,
                        cut_seq,
                        orientation,
                        base_mono.hor_idx,
                        base_mono.sequence
                    ).unwrap();
                }

                // Write summary line (one row per array with HOR and monomer statistics + consensus)
                {
                    // Compute monomer-level statistics
                    let mono_mean_ed_tmpl: f64 = {
                        let values: Vec<f64> = rec.base_monomers.iter()
                            .filter_map(|m| m.ed_tmpl.map(|v| v as f64))
                            .collect();
                        if !values.is_empty() { values.iter().sum::<f64>() / values.len() as f64 } else { 0.0 }
                    };
                    let mono_mean_ed_prev: f64 = {
                        let values: Vec<f64> = rec.base_monomers.iter()
                            .filter_map(|m| m.ed_prev.map(|v| v as f64))
                            .collect();
                        if !values.is_empty() { values.iter().sum::<f64>() / values.len() as f64 } else { 0.0 }
                    };
                    let mono_cv: f64 = {
                        let values: Vec<f64> = rec.base_monomers.iter().map(|m| m.cv).collect();
                        if !values.is_empty() { values.iter().sum::<f64>() / values.len() as f64 } else { 0.0 }
                    };

                    // Use "-" for empty consensus
                    let hor_consensus = if result.consensus_seq.is_empty() { "-".to_string() } else { result.consensus_seq.clone() };
                    let hor_iupac = if result.iupac_str.is_empty() { "-".to_string() } else { result.iupac_str.clone() };
                    let hor_quality = if result.quality_str.is_empty() { "-".to_string() } else { result.quality_str.clone() };
                    let mono_consensus = if rec.consensus_seq.is_empty() { "-".to_string() } else { rec.consensus_seq.clone() };
                    let mono_iupac = if rec.iupac_str.is_empty() { "-".to_string() } else { rec.iupac_str.clone() };
                    let mono_quality = if rec.quality_str.is_empty() { "-".to_string() } else { rec.quality_str.clone() };

                    writeln!(
                        fw_summary, "{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{:.2}\t{:.2}\t{:.4}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{:.2}\t{:.2}\t{:.4}\t{}\t{}\t{}\t{}",
                        result.header,
                        result.array_len,
                        orientation,
                        result.method,
                        result.period,
                        result.autocorr_value.unwrap_or(0.0),
                        result.n_monomers,
                        result.mean_ed_tmpl,
                        result.mean_ed_prev,
                        result.cv,
                        hor_consensus,
                        hor_iupac,
                        hor_quality,
                        rec.median_period,
                        rec.mean_autocorr,
                        rec.n_expected,
                        mono_mean_ed_tmpl,
                        mono_mean_ed_prev,
                        mono_cv,
                        mono_consensus,
                        mono_iupac,
                        mono_quality,
                        cut_seq,
                    ).unwrap();
                }

                // Write FASTA
                let header_info = format!(
                    "{} cut={} orientation={} n_monomers={} period={} method={}",
                    result.header, cut_seq, orientation, result.n_monomers, result.period, result.method
                );
                writeln!(fw, ">{}", header_info).unwrap();
                let seqs: Vec<&str> = result.monomers.iter().map(|m| m.sequence.as_str()).collect();
                writeln!(fw, "{}", seqs.join(" ")).unwrap();

                // Write lengths
                writeln!(fw_lengths, ">{}", header_info).unwrap();
                let lengths: Vec<String> = result.monomers.iter().map(|m| m.sequence.len().to_string()).collect();
                writeln!(fw_lengths, "{}", lengths.join(" ")).unwrap();

                // Collect stats for this array (only if requested)
                if collect_stats {
                    let mut array_lengths: Vec<usize> = Vec::new();
                    let mut array_eds: Vec<usize> = Vec::new();
                    for mono in &result.monomers {
                        if mono.piece_type == "monomer" {
                            array_lengths.push(mono.sequence.len());
                            if let Some(ed) = mono.ed_prev {
                                array_eds.push(ed);
                            }
                        }
                    }

                    if !array_lengths.is_empty() {
                        all_stats.push(ArrayStats {
                            id: result.header.clone(),
                            n_monomers: array_lengths.len(),
                            lengths: array_lengths,
                            eds: array_eds,
                        });
                    }
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

    // Flush all writers before sorting
    drop(fw);
    drop(fw_hors);
    drop(fw_monomers);
    drop(fw_lengths);
    drop(fw_summary);

    let total = total_count.load(Ordering::SeqCst);
    eprintln!("\rProcessed: {}/{} (100.0%)", processed, total);
    eprintln!("Total monomers: {}", total_monomers);

    // Sort all TSV files by genomic position
    eprintln!("Sorting output files...");
    sort_tsv_with_type_priority(&hors_file);     // With type order: pred_array -> monomer -> array -> consensus
    sort_tsv_with_type_priority(&monomers_file); // With type order: pred_array -> base_monomer
    sort_tsv_file(&summary_file);                // Simple sort by array_id (one row per array)
    eprintln!("Sorting complete.");

    // Output statistics if requested
    if collect_stats && !all_stats.is_empty() {
        eprintln!("\n=== Statistics ===");

        let mut stats_with_rcv: Vec<(String, usize, usize, f64, f64)> = all_stats.iter().map(|s| {
            let rcv = robust_cv(&s.lengths);
            let median_len = {
                let mut sorted = s.lengths.clone();
                sorted.sort();
                sorted[sorted.len() / 2]
            };
            let mean_ed = if s.eds.is_empty() {
                0.0
            } else {
                s.eds.iter().sum::<usize>() as f64 / s.eds.len() as f64
            };
            (s.id.clone(), s.n_monomers, median_len, rcv, mean_ed)
        }).collect();

        let excellent = stats_with_rcv.iter().filter(|s| s.3 <= 0.05).count();
        let moderate = stats_with_rcv.iter().filter(|s| s.3 > 0.05 && s.3 <= 0.20).count();
        let poor = stats_with_rcv.iter().filter(|s| s.3 > 0.20).count();
        let total_arrays = stats_with_rcv.len();

        eprintln!("Arrays: {}", total_arrays);
        eprintln!("\nQuality (Robust CV = MAD/median):");
        eprintln!("  Excellent (RCV ≤ 5%):  {:>4} ({:>2}%)", excellent, 100 * excellent / total_arrays);
        eprintln!("  Moderate  (5-20%):     {:>4} ({:>2}%)", moderate, 100 * moderate / total_arrays);
        eprintln!("  Poor      (> 20%):     {:>4} ({:>2}%)", poor, 100 * poor / total_arrays);

        stats_with_rcv.sort_by(|a, b| b.3.partial_cmp(&a.3).unwrap_or(std::cmp::Ordering::Equal));

        eprintln!("\nTop {} outliers (highest RCV):", top_outliers.min(stats_with_rcv.len()));
        eprintln!("{:<55} {:>5} {:>7} {:>7} {:>8}", "Array", "N", "Median", "RCV", "MeanED");
        for (id, n, median, rcv, mean_ed) in stats_with_rcv.iter().take(top_outliers) {
            eprintln!("{:<55} {:>5} {:>7} {:>7.3} {:>8.1}", id, n, median, rcv, mean_ed);
        }

        stats_with_rcv.sort_by(|a, b| a.3.partial_cmp(&b.3).unwrap_or(std::cmp::Ordering::Equal));
        eprintln!("\nTop {} best (lowest RCV):", top_outliers.min(stats_with_rcv.len()));
        eprintln!("{:<55} {:>5} {:>7} {:>7} {:>8}", "Array", "N", "Median", "RCV", "MeanED");
        for (id, n, median, rcv, mean_ed) in stats_with_rcv.iter().take(top_outliers) {
            eprintln!("{:<55} {:>5} {:>7} {:>7.3} {:>8.1}", id, n, median, rcv, mean_ed);
        }
    }
}

fn main() {
    let start_time = Instant::now();
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

    // Validate method
    let method = args.method.clone();
    if !["autocorr", "classic", "both"].contains(&method.as_str()) {
        eprintln!("Error: --method must be 'autocorr', 'classic', or 'both'");
        std::process::exit(1);
    }
    eprintln!("Method: {}", method);

    // Set global MAX_ED_LEN from CLI
    MAX_ED_LEN.store(args.max_ed_len, Ordering::Relaxed);
    eprintln!("Max ED length: {} bp", args.max_ed_len);

    // Resource monitoring: track peak memory and CPU
    let peak_memory_kb = Arc::new(AtomicU64::new(0));
    let total_cpu_samples = Arc::new(AtomicU64::new(0));
    let total_cpu_sum = Arc::new(AtomicU64::new(0)); // sum of cpu_usage * 100 (fixed-point)
    let monitor_running = Arc::new(std::sync::atomic::AtomicBool::new(true));

    // Spawn memory/CPU monitor thread
    let peak_mem = Arc::clone(&peak_memory_kb);
    let cpu_samples = Arc::clone(&total_cpu_samples);
    let cpu_sum = Arc::clone(&total_cpu_sum);
    let monitor_flag = Arc::clone(&monitor_running);
    let monitor_handle = thread::spawn(move || {
        let mut sys = System::new();
        let pid = Pid::from_u32(std::process::id());

        // Initial refresh to set baseline for cpu_usage()
        sys.refresh_processes_specifics(
            ProcessRefreshKind::new().with_memory().with_cpu()
        );
        thread::sleep(std::time::Duration::from_millis(200));

        while monitor_flag.load(Ordering::Relaxed) {
            sys.refresh_processes_specifics(
                ProcessRefreshKind::new().with_memory().with_cpu()
            );

            if let Some(process) = sys.process(pid) {
                let mem_kb = process.memory() / 1024;
                let current_peak = peak_mem.load(Ordering::Relaxed);
                if mem_kb > current_peak {
                    peak_mem.store(mem_kb, Ordering::Relaxed);
                }

                // Track CPU usage (as percentage * 100 for precision)
                let cpu = process.cpu_usage();
                if cpu > 0.0 {
                    cpu_sum.fetch_add((cpu * 100.0) as u64, Ordering::Relaxed);
                    cpu_samples.fetch_add(1, Ordering::Relaxed);
                }
            }
            thread::sleep(std::time::Duration::from_millis(500));
        }
    });

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
        let worker_method = method.clone();

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
                            &worker_method,
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
    let collect_stats = args.stats;
    let top_outliers = args.top_outliers;
    let writer_handle = thread::spawn(move || {
        writer_thread(output_rx, output_prefix, writer_total, num_workers, collect_stats, top_outliers);
    });

    // Wait for completion
    reader_handle.join().expect("Reader thread panicked");
    for handle in worker_handles {
        handle.join().expect("Worker thread panicked");
    }
    writer_handle.join().expect("Writer thread panicked");

    // Stop memory monitor
    monitor_running.store(false, Ordering::Relaxed);
    monitor_handle.join().expect("Monitor thread panicked");

    // Calculate final stats
    let elapsed = start_time.elapsed();
    let peak_mem = peak_memory_kb.load(Ordering::Relaxed);

    // Compute average CPU usage from monitor samples
    let samples = total_cpu_samples.load(Ordering::Relaxed);
    let cpu_sum_val = total_cpu_sum.load(Ordering::Relaxed);
    let cpu_usage = if samples > 0 {
        cpu_sum_val as f64 / samples as f64 / 100.0  // Convert back from fixed-point
    } else {
        0.0
    };

    // Format elapsed time
    let total_secs = elapsed.as_secs();
    let hours = total_secs / 3600;
    let minutes = (total_secs % 3600) / 60;
    let seconds = total_secs % 60;
    let millis = elapsed.subsec_millis();

    eprintln!("\n=== Resource Usage ===");
    if hours > 0 {
        eprintln!("Elapsed time:  {}h {}m {}.{:03}s", hours, minutes, seconds, millis);
    } else if minutes > 0 {
        eprintln!("Elapsed time:  {}m {}.{:03}s", minutes, seconds, millis);
    } else {
        eprintln!("Elapsed time:  {}.{:03}s", seconds, millis);
    }

    // Format memory
    if peak_mem > 1024 * 1024 {
        eprintln!("Peak memory:   {:.2} GB", peak_mem as f64 / 1024.0 / 1024.0);
    } else if peak_mem > 1024 {
        eprintln!("Peak memory:   {:.2} MB", peak_mem as f64 / 1024.0);
    } else {
        eprintln!("Peak memory:   {} KB", peak_mem);
    }

    eprintln!("CPU usage:     {:.1}%", cpu_usage);
    eprintln!("Threads:       {}", num_workers);
    eprintln!("Done!");
}
