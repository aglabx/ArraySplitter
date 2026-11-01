use std::collections::{BTreeMap, HashMap};
use std::env;
use std::fs;
use std::io::{BufRead, BufReader};
use std::process::Command;

#[allow(dead_code)]
#[derive(Debug, Clone)]
struct AsplitRow {
    array_id: String,
    length: usize,
    hor_period: usize,
    hor_n: usize,
    mono_period: usize,
    mono_n: usize,
    mono_ed_tmpl: String,
    status: String,
}

fn read_asplit_summary(path: &str) -> Vec<AsplitRow> {
    let f = fs::File::open(path).unwrap_or_else(|e| panic!("Failed to open {}: {}", path, e));
    let reader = BufReader::new(f);
    let mut lines = reader.lines();
    let header_line = lines.next().expect("Empty summary.tsv").expect("Read err");
    let headers: Vec<&str> = header_line.split('\t').collect();
    let mut col = HashMap::new();
    for (i, &h) in headers.iter().enumerate() {
        col.insert(h.to_string(), i);
    }

    let mut out = Vec::new();
    for line in lines {
        let l = line.expect("Read err");
        if l.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = l.split('\t').collect();
        let get = |k: &str| -> &str {
            col.get(k).and_then(|&idx| parts.get(idx)).copied().unwrap_or("")
        };
        let array_id = get("array_id").to_string();
        let length = get("array_length").parse::<usize>().unwrap_or(0);
        let hor_period = get("hor_period").parse::<usize>().unwrap_or(0);
        let hor_n = get("hor_n_monomers").parse::<usize>().unwrap_or(0);
        let mono_period = get("mono_period").parse::<usize>().unwrap_or(0);
        let mono_n = get("mono_n_monomers").parse::<usize>().unwrap_or(0);
        let mono_ed_tmpl = get("mono_mean_ed_tmpl").to_string();
        let status = get("status").to_string();
        out.push(AsplitRow {
            array_id,
            length,
            hor_period,
            hor_n,
            mono_period,
            mono_n,
            mono_ed_tmpl,
            status,
        });
    }
    out
}

/// Parse FasTAN .1ano file via ONEview binary; returns map array_id -> (dominant_period, total_units, covered_bp)
fn parse_fastan_ano(oneview_bin: &str, ano_path: &str) -> HashMap<String, (usize, usize, usize)> {
    let output = Command::new(oneview_bin)
        .arg(ano_path)
        .output()
        .unwrap_or_else(|e| panic!("Failed to run ONEview on {}: {}", ano_path, e));
    if !output.status.success() {
        panic!("ONEview failed on {}: {}", ano_path, String::from_utf8_lossy(&output.stderr));
    }
    let txt = String::from_utf8_lossy(&output.stdout);

    let mut names: Vec<String> = Vec::new();
    let mut cur_aid: Option<String> = None;
    let mut cur_beg: usize = 0;
    let mut cur_end: usize = 0;
    let mut cur_unit: Option<usize> = None;
    let mut cur_ncuts: usize = 0;

    // array_id -> vec of (span, period, n_units)
    let mut intervals: HashMap<String, Vec<(usize, usize, usize)>> = HashMap::new();

    let flush = |aid: &Option<String>,
                 beg: usize,
                 end: usize,
                 unit: Option<usize>,
                 ncuts: usize,
                 map: &mut HashMap<String, Vec<(usize, usize, usize)>>| {
        if let (Some(a), Some(u)) = (aid, unit) {
            if end > beg {
                let span = end - beg;
                let n_u = if ncuts > 0 { ncuts } else { span / u.max(1) };
                map.entry(a.clone()).or_default().push((span, u, n_u));
            }
        }
    };

    for line in txt.lines() {
        if let Some(rest) = line.strip_prefix("S ") {
            let mut sp = rest.splitn(2, ' ');
            let _len = sp.next();
            if let Some(name) = sp.next() {
                let id = name.split_whitespace().next().unwrap_or(name).to_string();
                names.push(id);
            }
        } else if let Some(rest) = line.strip_prefix("M ") {
            flush(&cur_aid, cur_beg, cur_end, cur_unit, cur_ncuts, &mut intervals);
            let f: Vec<&str> = rest.split_whitespace().collect();
            if f.len() >= 3 {
                let idx = f[0].parse::<usize>().unwrap_or(0);
                cur_aid = names.get(idx).cloned();
                cur_beg = f[1].parse::<usize>().unwrap_or(0);
                cur_end = f[2].parse::<usize>().unwrap_or(0);
                cur_unit = None;
                cur_ncuts = 0;
            }
        } else if let Some(rest) = line.strip_prefix("L ") {
            let f: Vec<&str> = rest.split_whitespace().collect();
            if f.len() >= 2 {
                cur_unit = f[1].parse::<usize>().ok();
            }
        } else if let Some(rest) = line.strip_prefix("P ") {
            let f: Vec<&str> = rest.split_whitespace().collect();
            if !f.is_empty() {
                cur_ncuts = f[0].parse::<usize>().unwrap_or(0);
            }
        }
    }
    flush(&cur_aid, cur_beg, cur_end, cur_unit, cur_ncuts, &mut intervals);

    let mut out = HashMap::new();
    for (aid, list) in intervals {
        let mut by_period: BTreeMap<usize, usize> = BTreeMap::new();
        let mut total_units = 0usize;
        let mut total_span = 0usize;
        for &(span, u, n_u) in &list {
            *by_period.entry(u).or_insert(0) += span;
            total_units += n_u;
            total_span += span;
        }
        let dom_p = by_period
            .into_iter()
            .max_by_key(|&(_, sp)| sp)
            .map(|(u, _)| u)
            .unwrap_or(0);
        out.insert(aid, (dom_p, total_units, total_span));
    }
    out
}

/// Parse centroAnno directory (`*_decomposedResult.csv` and `*_horDecomposedResult.csv`).
/// Returns map array_id -> (mono_period, n_monomers, hor_period_opt)
fn parse_centroanno_dir(dir: &str) -> HashMap<String, (usize, usize, usize)> {
    let mut mono_lengths: HashMap<String, Vec<usize>> = HashMap::new();
    let mut span_by_aid: HashMap<String, BTreeMap<usize, usize>> = HashMap::new();

    let entries = fs::read_dir(dir).unwrap_or_else(|e| panic!("centroAnno dir {}: {}", dir, e));
    {
        for entry in entries.flatten() {
            let path = entry.path();
            let Some(fname) = path.file_name().and_then(|s| s.to_str()) else {
                continue;
            };
            if fname.ends_with("_decomposedResult.csv") && !fname.ends_with("_horDecomposedResult.csv") {
                if let Ok(f) = fs::File::open(&path) {
                    let reader = BufReader::new(f);
                    for line in reader.lines().flatten() {
                        if line.starts_with('#') || line.trim().is_empty() {
                            continue;
                        }
                        let p: Vec<&str> = line.split(',').collect();
                        if p.len() >= 6 {
                            let aid = p[0].to_string();
                            let mlen = p[5].trim().parse::<usize>().unwrap_or(0);
                            if mlen > 0 {
                                mono_lengths.entry(aid).or_default().push(mlen);
                            }
                        }
                    }
                }
            } else if fname.ends_with("_horDecomposedResult.csv") {
                if let Ok(f) = fs::File::open(&path) {
                    let reader = BufReader::new(f);
                    for line in reader.lines().flatten() {
                        if line.starts_with('#') || line.trim().is_empty() {
                            continue;
                        }
                        let p: Vec<&str> = line.split(',').collect();
                        if p.len() >= 6 {
                            let aid = p[0].to_string();
                            let s = p[2].trim().parse::<usize>().unwrap_or(0);
                            let e = p[3].trim().parse::<usize>().unwrap_or(0);
                            let hl = p[5].trim().parse::<usize>().unwrap_or(0);
                            let sp = e.saturating_sub(s).max(hl);
                            if hl > 0 && sp > 0 {
                                *span_by_aid.entry(aid).or_default().entry(hl).or_insert(0) += sp;
                            }
                        }
                    }
                }
            }
        }
    }

    let mut hor_periods: HashMap<String, usize> = HashMap::new();
    for (aid, m) in span_by_aid {
        if let Some((hl, _)) = m.into_iter().max_by_key(|&(_, sp)| sp) {
            hor_periods.insert(aid, hl);
        }
    }

    let mut out = HashMap::new();
    for (aid, mut lens) in mono_lengths {
        let n = lens.len();
        lens.sort_unstable();
        let mut counts: BTreeMap<usize, usize> = BTreeMap::new();
        for &l in &lens {
            *counts.entry(l).or_insert(0) += 1;
        }
        let mode_len = counts
            .into_iter()
            .max_by_key(|&(_, c)| c)
            .map(|(l, _)| l)
            .unwrap_or_else(|| lens[n / 2]);
        let hp = hor_periods.get(&aid).copied().unwrap_or(0);
        out.insert(aid, (mode_len, n, hp));
    }
    out
}

/// Parse /usr/bin/time -v file -> (wall_seconds, max_rss_mib). A missing file or field is an error.
fn parse_time_file(path: &str) -> (f64, f64) {
    let txt = fs::read_to_string(path).unwrap_or_else(|e| panic!("time file {}: {}", path, e));
    let mut wall_s: Option<f64> = None;
    let mut rss_kb: Option<f64> = None;
    for line in txt.lines() {
        let l = line.trim();
        if l.starts_with("Elapsed (wall clock) time") {
            let val = l.split_whitespace().last().unwrap_or_else(|| panic!("{}: empty wall-clock line", path));
            let mut sec = 0.0;
            for part in val.split(':') {
                sec = sec * 60.0 + part.parse::<f64>().unwrap_or_else(|_| panic!("{}: bad wall clock {}", path, val));
            }
            wall_s = Some(sec);
        } else if l.starts_with("Maximum resident set size") {
            let val = l.split_whitespace().last().unwrap_or_else(|| panic!("{}: empty RSS line", path));
            rss_kb = Some(val.parse::<f64>().unwrap_or_else(|_| panic!("{}: bad RSS {}", path, val)));
        }
    }
    (
        wall_s.unwrap_or_else(|| panic!("{}: no wall-clock line", path)),
        rss_kb.unwrap_or_else(|| panic!("{}: no RSS line", path)) / 1024.0,
    )
}

/// Sum the `calls` column (6th) of data/fields_<h>.tsv. A missing file or bad number is an error.
fn sum_nhmmer_from_fields_tsv(tsv_path: &str) -> usize {
    let f = fs::File::open(tsv_path).unwrap_or_else(|e| panic!("{}: {}", tsv_path, e));
    let mut total = 0usize;
    for (i, line) in BufReader::new(f).lines().enumerate() {
        let line = line.unwrap_or_else(|e| panic!("{}: {}", tsv_path, e));
        if i == 0 || line.trim().is_empty() {
            continue;
        }
        let p: Vec<&str> = line.split('\t').collect();
        let c = p.get(5).unwrap_or_else(|| panic!("{} line {}: fewer than 6 columns", tsv_path, i + 1));
        total += c.trim().parse::<usize>().unwrap_or_else(|_| panic!("{} line {}: bad calls {}", tsv_path, i + 1, c));
    }
    total
}

/// One field as seen by one tool, with the same definitions for every tool (Table 5).
struct FieldCall {
    max_period: usize, // longest period the tool reports for the field
    has_113: bool,     // a 108-116-bp unit is reported
    nested: bool,      // a period >= 300 bp is reported above a 108-116-bp unit (two levels)
    n113: usize,       // units reported at the 108-116-bp level
}

fn in113(p: usize) -> bool {
    (108..=116).contains(&p)
}

fn asplit_calls(rows: &[AsplitRow]) -> HashMap<String, FieldCall> {
    rows.iter()
        .filter(|r| r.status != "no_period")
        .map(|r| {
            let has = in113(r.mono_period);
            (
                r.array_id.clone(),
                FieldCall {
                    max_period: r.hor_period.max(r.mono_period),
                    has_113: has,
                    nested: has && r.hor_period >= 300,
                    n113: if has { r.mono_n } else { 0 },
                },
            )
        })
        .collect()
}

/// FasTAN reports one level per interval; the field's period is the unit covering most of it.
fn fastan_calls(m: &HashMap<String, (usize, usize, usize)>) -> HashMap<String, FieldCall> {
    m.iter()
        .map(|(a, &(p, n, _))| {
            (a.clone(), FieldCall { max_period: p, has_113: in113(p), nested: false, n113: if in113(p) { n } else { 0 } })
        })
        .collect()
}

/// centroAnno: modal monomer length and the HOR length covering most of the field.
fn centro_calls(m: &HashMap<String, (usize, usize, usize)>) -> HashMap<String, FieldCall> {
    m.iter()
        .map(|(a, &(mono, n, hor))| {
            let has = in113(mono);
            (
                a.clone(),
                FieldCall { max_period: mono.max(hor), has_113: has, nested: has && hor >= 300, n113: if has { n } else { 0 } },
            )
        })
        .collect()
}

/// One Table 5 row; `out` collects the TSV rows for table5.tsv.
fn table5_row(label: &str, ids: &[String], calls: &HashMap<String, FieldCall>, time_path: &str, nhmmer: usize, out: &mut Vec<String>) {
    let got: Vec<&FieldCall> = ids.iter().filter_map(|id| calls.get(id)).collect();
    let no_call = ids.len() - got.len();
    let ge300 = got.iter().filter(|c| c.max_period >= 300).count();
    let u113 = got.iter().filter(|c| c.has_113).count();
    let nested = got.iter().filter(|c| c.nested).count();
    let n113: usize = got.iter().map(|c| c.n113).sum();
    let (ws, rss) = parse_time_file(time_path);
    let row = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}\t{:.2}\t{:.1}",
        label, ids.len(), no_call, ge300, u113, nested, n113, 100.0 * n113 as f64 / nhmmer as f64, ws, rss
    );
    println!("{}", row);
    out.push(row);
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!("usage: cat113_constellation_eval <workdir of run_cat113.sh> <ONEview>");
        std::process::exit(2);
    }
    let work = args[1].as_str();
    let oneview = args[2].as_str();
    let bench_dir = work.to_string();

    // 1. Evaluate 7-array panel across all tools
    let new_cat7_d32k = read_asplit_summary(&format!(
        "{}/asplit_new/cat7_d32k.summary.tsv",
        work
    ));
    let new_cat7_d5k = read_asplit_summary(&format!(
        "{}/asplit_new/cat7_d5k.summary.tsv",
        work
    ));
    let v182_cat7_d32k = read_asplit_summary(&format!(
        "{}/asplit_v182/cat7_d32k.summary.tsv",
        work
    ));
    let v182_cat7_d5k = read_asplit_summary(&format!(
        "{}/asplit_v182/cat7_d5k.summary.tsv",
        work
    ));
    let v182_cat7_rt03 = read_asplit_summary(&format!(
        "{}/asplit_v182/cat7_d32k_rt0.3.summary.tsv",
        work
    ));
    let fastan_cat7 = parse_fastan_ano(oneview, &format!("{}/fastan/cat7.1ano", bench_dir));
    let centro_cat7 = parse_centroanno_dir(&format!("{}/centroanno/cat7", bench_dir));

    println!("## CAT_7_ARRAY_PANEL_COMPARISON");
    println!("array_id\tlen_bp\tv182_d5k(H/M)\tv182_d32k(H/M)\tv182_rt0.3(H/M)\tconstel_d5k(H/M)\tconstel_d32k(H/M)\tconstel_n_mono\tFasTAN(P/n)\tcentroAnno(M/H)");

    let map_by_id = |v: &[AsplitRow]| -> HashMap<String, AsplitRow> {
        v.iter().map(|r| (r.array_id.clone(), r.clone())).collect()
    };
    let m_new32 = map_by_id(&new_cat7_d32k);
    let m_new5 = map_by_id(&new_cat7_d5k);
    let m_182_32 = map_by_id(&v182_cat7_d32k);
    let m_182_5 = map_by_id(&v182_cat7_d5k);
    let m_182_rt = map_by_id(&v182_cat7_rt03);

    let labels = [
        "E1_113_p7279",
        "F1_113_p2647",
        "D4_113_p3533",
        "A3_113_p15165",
        "E2_113_p11906",
        "F2_113_p3122",
        "A3_483_control",
    ];
    for &lbl in &labels {
        let r32 = &m_new32[lbl];
        let r5 = &m_new5[lbl];
        let o32 = &m_182_32[lbl];
        let o5 = &m_182_5[lbl];
        let ort = &m_182_rt[lbl];
        let ft = fastan_cat7.get(lbl).copied().unwrap_or((0, 0, 0));
        let ca = centro_cat7.get(lbl).copied().unwrap_or((0, 0, 0));
        println!(
            "{}\t{}\t{}/{}\t{}/{}\t{}/{}\t{}/{}\t{}/{}\t{}\t{}/{}\t{}/{}",
            lbl,
            r32.length,
            o5.hor_period,
            o5.mono_period,
            o32.hor_period,
            o32.mono_period,
            ort.hor_period,
            ort.mono_period,
            r5.hor_period,
            r5.mono_period,
            r32.hor_period,
            r32.mono_period,
            r32.mono_n,
            ft.0,
            ft.1,
            ca.0,
            ca.2
        );
    }

    // 2. Evaluate genome-wide cat_pat1.0 (41 fields) and cat_mat1.0 (51 fields)
    let new_pat_32 = read_asplit_summary(&format!("{}/asplit_new/pat_d32k.summary.tsv", bench_dir));
    let new_pat_5 = read_asplit_summary(&format!("{}/asplit_new/pat_d5k.summary.tsv", bench_dir));
    let old_pat_32 = read_asplit_summary(&format!("{}/asplit_v182/pat_d32k.summary.tsv", bench_dir));
    let old_pat_5 = read_asplit_summary(&format!("{}/asplit_v182/pat_d5k.summary.tsv", bench_dir));

    let new_mat_32 = read_asplit_summary(&format!("{}/asplit_new/mat_d32k.summary.tsv", bench_dir));
    let new_mat_5 = read_asplit_summary(&format!("{}/asplit_new/mat_d5k.summary.tsv", bench_dir));
    let old_mat_32 = read_asplit_summary(&format!("{}/asplit_v182/mat_d32k.summary.tsv", bench_dir));
    let old_mat_5 = read_asplit_summary(&format!("{}/asplit_v182/mat_d5k.summary.tsv", bench_dir));

    let pat_ids: Vec<String> = new_pat_32.iter().map(|r| r.array_id.clone()).collect();
    let mat_ids: Vec<String> = new_mat_32.iter().map(|r| r.array_id.clone()).collect();

    let nhmmer_pat_in_fields = sum_nhmmer_from_fields_tsv(&format!(
        "{}/../data/fields_pat.tsv",
        env!("CARGO_MANIFEST_DIR")
    ));
    let nhmmer_mat_in_fields = sum_nhmmer_from_fields_tsv(&format!(
        "{}/../data/fields_mat.tsv",
        env!("CARGO_MANIFEST_DIR")
    ));

    // 2-3. Table 5: every tool with the same field definitions (struct FieldCall)
    let fastan_pat = parse_fastan_ano(oneview, &format!("{}/fastan/pat.1ano", bench_dir));
    let fastan_mat = parse_fastan_ano(oneview, &format!("{}/fastan/mat.1ano", bench_dir));
    let centro_pat = parse_centroanno_dir(&format!("{}/centroanno/pat", bench_dir));
    let centro_mat = parse_centroanno_dir(&format!("{}/centroanno/mat", bench_dir));

    println!(
        "\n## TABLE5 (nhmmer calls in fields: pat {}, mat {}; also written to {}/table5.tsv)",
        nhmmer_pat_in_fields, nhmmer_mat_in_fields, work
    );
    let header = "run\tfields\tno_call\tperiod_ge300\tunit_108_116\tnested_ge300_over_108_116\tunits_108_116\tpct_of_nhmmer\twall_s\trss_mib";
    println!("{}", header);
    let mut rows: Vec<String> = Vec::new();
    for (hap, ids, nh, runs) in [
        ("pat", &pat_ids, nhmmer_pat_in_fields, [(&old_pat_5, "v182_single_anchor_D5000", "asplit_v182_pat_d5k"), (&old_pat_32, "v182_single_anchor_D32000", "asplit_v182_pat_d32k"),
            (&new_pat_5, "constellation_D5000", "asplit_new_pat_d5k"), (&new_pat_32, "constellation_D32000", "asplit_new_pat_d32k")]),
        ("mat", &mat_ids, nhmmer_mat_in_fields, [(&old_mat_5, "v182_single_anchor_D5000", "asplit_v182_mat_d5k"), (&old_mat_32, "v182_single_anchor_D32000", "asplit_v182_mat_d32k"),
            (&new_mat_5, "constellation_D5000", "asplit_new_mat_d5k"), (&new_mat_32, "constellation_D32000", "asplit_new_mat_d32k")]),
    ] {
        for (r, tag, tf) in runs {
            table5_row(&format!("{}_ArraySplitter_{}", hap, tag), ids, &asplit_calls(r), &format!("{}/logs/{}.time", bench_dir, tf), nh, &mut rows);
        }
        let (fm, cm) = if hap == "pat" { (&fastan_pat, &centro_pat) } else { (&fastan_mat, &centro_mat) };
        table5_row(&format!("{}_FasTAN", hap), ids, &fastan_calls(fm), &format!("{}/logs/fastan_{}.time", bench_dir, hap), nh, &mut rows);
        table5_row(&format!("{}_centroAnno", hap), ids, &centro_calls(cm), &format!("{}/logs/centro_{}.time", bench_dir, hap), nh, &mut rows);
        let c113: Vec<String> = ids.iter().filter(|id| centro_calls(cm).get(*id).map_or(false, |c| c.has_113))
            .map(|id| id.split('|').next().unwrap_or(id).to_string()).collect();
        println!("# {} centroAnno fields with a 108-116-bp unit: {}", hap, c113.join(", "));
    }
    let mut body = String::from(header);
    body.push('\n');
    for r in &rows {
        body.push_str(r);
        body.push('\n');
    }
    fs::write(format!("{}/table5.tsv", work), body).unwrap_or_else(|e| panic!("write table5.tsv: {}", e));
}

