//! I/O utilities for reading and writing FASTA, TSV, and other formats

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// A FASTA record
#[derive(Debug, Clone)]
pub struct FastaRecord {
    pub id: String,
    pub description: String,
    pub sequence: String,
}

impl FastaRecord {
    pub fn new(id: &str, description: &str, sequence: &str) -> Self {
        Self {
            id: id.to_string(),
            description: description.to_string(),
            sequence: sequence.to_string(),
        }
    }

    /// Get full header (id + description)
    pub fn header(&self) -> String {
        if self.description.is_empty() {
            self.id.clone()
        } else {
            format!("{} {}", self.id, self.description)
        }
    }
}

/// Read all records from a FASTA file
pub fn read_fasta<P: AsRef<Path>>(path: P) -> std::io::Result<Vec<FastaRecord>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let mut records = Vec::new();
    let mut current_id = String::new();
    let mut current_desc = String::new();
    let mut current_seq = String::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();

        if trimmed.starts_with('>') {
            // Save previous record if exists
            if !current_id.is_empty() {
                records.push(FastaRecord::new(&current_id, &current_desc, &current_seq));
            }

            // Parse new header
            let header = &trimmed[1..];
            let parts: Vec<&str> = header.splitn(2, |c| c == ' ' || c == '\t').collect();
            current_id = parts[0].to_string();
            current_desc = parts.get(1).unwrap_or(&"").to_string();
            current_seq = String::new();
        } else if !trimmed.is_empty() {
            current_seq.push_str(trimmed);
        }
    }

    // Don't forget the last record
    if !current_id.is_empty() {
        records.push(FastaRecord::new(&current_id, &current_desc, &current_seq));
    }

    Ok(records)
}

/// Iterator over FASTA records (memory-efficient for large files)
pub struct FastaReader<R: BufRead> {
    reader: R,
    buffer: String,
    current_id: String,
    current_desc: String,
}

impl<R: BufRead> FastaReader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            buffer: String::new(),
            current_id: String::new(),
            current_desc: String::new(),
        }
    }
}

impl<R: BufRead> Iterator for FastaReader<R> {
    type Item = std::io::Result<FastaRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut sequence = String::new();
        let mut found_header = false;

        loop {
            self.buffer.clear();
            match self.reader.read_line(&mut self.buffer) {
                Ok(0) => {
                    // EOF
                    if !self.current_id.is_empty() && !sequence.is_empty() {
                        let record = FastaRecord::new(&self.current_id, &self.current_desc, &sequence);
                        self.current_id.clear();
                        return Some(Ok(record));
                    }
                    return None;
                }
                Ok(_) => {
                    let trimmed = self.buffer.trim();

                    if trimmed.starts_with('>') {
                        // Return previous record if exists
                        if found_header || !self.current_id.is_empty() {
                            let record = FastaRecord::new(&self.current_id, &self.current_desc, &sequence);

                            // Parse new header
                            let header = &trimmed[1..];
                            let parts: Vec<&str> = header.splitn(2, |c| c == ' ' || c == '\t').collect();
                            self.current_id = parts[0].to_string();
                            self.current_desc = parts.get(1).unwrap_or(&"").to_string();

                            if !record.id.is_empty() {
                                return Some(Ok(record));
                            }
                        } else {
                            // First header
                            let header = &trimmed[1..];
                            let parts: Vec<&str> = header.splitn(2, |c| c == ' ' || c == '\t').collect();
                            self.current_id = parts[0].to_string();
                            self.current_desc = parts.get(1).unwrap_or(&"").to_string();
                            found_header = true;
                        }
                    } else if !trimmed.is_empty() {
                        sequence.push_str(trimmed);
                    }
                }
                Err(e) => return Some(Err(e)),
            }
        }
    }
}

/// Write FASTA records to a file
pub fn write_fasta<P: AsRef<Path>>(path: P, records: &[FastaRecord], line_width: usize) -> std::io::Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    for record in records {
        writeln!(writer, ">{}", record.header())?;

        // Write sequence with line wrapping
        for chunk in record.sequence.as_bytes().chunks(line_width) {
            writeln!(writer, "{}", std::str::from_utf8(chunk).unwrap())?;
        }
    }

    Ok(())
}

/// Monomer record for TSV output
#[derive(Debug, Clone)]
pub struct MonomerRecord {
    pub array_id: String,
    pub index: usize,
    pub monomer_type: String,  // "left_flank", "monomer", "right_flank"
    pub length: usize,
    pub orientation: String,   // "fwd" or "rev"
    pub sequence: String,
}

/// Write monomers TSV file (compatible with Python version)
pub fn write_monomers_tsv<P: AsRef<Path>>(path: P, monomers: &[MonomerRecord]) -> std::io::Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    // Header
    writeln!(writer, "array_id\tindex\ttype\tlength\torientation")?;

    for m in monomers {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}",
            m.array_id, m.index, m.monomer_type, m.length, m.orientation
        )?;
    }

    Ok(())
}

/// Write lengths file (simple format: header\tlength)
pub fn write_lengths<P: AsRef<Path>>(path: P, lengths: &[(String, usize)]) -> std::io::Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    for (header, length) in lengths {
        writeln!(writer, "{}\t{}", header, length)?;
    }

    Ok(())
}

/// Parse predefined cuts from comma-separated string
pub fn parse_cuts(cuts_str: &str) -> Vec<String> {
    cuts_str
        .split(',')
        .map(|s| s.trim().to_uppercase())
        .filter(|s| !s.is_empty())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_fasta_reader() {
        let fasta = ">seq1 description\nACGT\nACGT\n>seq2\nTTTT\n";
        let cursor = Cursor::new(fasta);
        let reader = FastaReader::new(BufReader::new(cursor));

        let records: Vec<_> = reader.collect::<Result<Vec<_>, _>>().unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id, "seq1");
        assert_eq!(records[0].description, "description");
        assert_eq!(records[0].sequence, "ACGTACGT");
        assert_eq!(records[1].id, "seq2");
        assert_eq!(records[1].sequence, "TTTT");
    }

    #[test]
    fn test_parse_cuts() {
        let cuts = parse_cuts("ATG, CGCG, ttaggg");
        assert_eq!(cuts, vec!["ATG", "CGCG", "TTAGGG"]);
    }
}
