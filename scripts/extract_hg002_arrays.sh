#!/usr/bin/env bash
# Extract HG002 maternal satellite arrays from T2T assembly
# Requires: samtools or bedtools

set -euo pipefail

ASSEMBLY="${1:-hg002v1.1.mat.fasta}"
BED="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/test_data/hg002_arrays.bed"
OUTPUT="${2:-hg002v1.1.mat.10kb.fasta}"

if [[ ! -f "${ASSEMBLY}" ]]; then
    echo "Usage: $0 <path_to_hg002v1.1.mat.fasta> [output.fasta]" >&2
    echo "Download assembly from: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/021/951/015/GCA_021951015.1_HG002.mat.cur.20211005/" >&2
    exit 1
fi

echo "Extracting 1,088 arrays using bedtools getfasta..."
bedtools getfasta -fi "${ASSEMBLY}" -bed "${BED}" -name -fo "${OUTPUT}"
echo "Saved to ${OUTPUT}"
