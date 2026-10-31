#!/usr/bin/env bash
# Monomer calls of the cat 113-bp satellite from the profile search (2026-10-30; PI: can alphabet annotate the cat
# monomers and find chromosome-specific sets from them; alphabet's own cut found none).
#
#   cen113_monomer_table.sh NHMMER_WORKDIR TWOBIT CHROM_LIST MODEL MIN_COLS OUTDIR
#
# NHMMER_WORKDIR is the WORKDIR of cen113_profile_search.sh (chrom/<chrom>.tbl). One call per nhmmer hit of MODEL
# whose alignment covers >= MIN_COLS model columns, oriented to the model strand (the 113-column E1 model: its
# columns 1-113 are one monomer, so a call is one monomer in one frame). Hits of the other model are not calls.
# OUTDIR/monomers.tsv: id, chrom, start (1-based), end, strand, hmm_from, hmm_to, bits, seq.
# OUTDIR/monomers.fa: the same sequences, id = chrom:start-end:strand.
# OUTDIR/coverage.tsv: per chromosome, hits of MODEL, hits kept, hits below MIN_COLS (named counts, not dropped
# silently). A chromosome in CHROM_LIST without a .tbl stops the run.
set -euo pipefail
NW=${1:?nhmmer workdir}; TB=${2:?2bit}; CL=${3:?chrom list}; MODEL=${4:?model name}; MINC=${5:?min model columns}
O=${6:?outdir}
command -v twoBitToFa > /dev/null || { echo "cen113_monomer_table: twoBitToFa not on PATH" >&2; exit 1; }
mkdir -p "$O"
: > "$O/calls.bed"
printf 'chrom\tmodel_hits\tkept\tbelow_min_cols\n' > "$O/coverage.tsv"
while read -r c; do
  [ -n "$c" ] || continue
  [ -s "$NW/chrom/$c.tbl" ] || { echo "cen113_monomer_table: $NW/chrom/$c.tbl missing" >&2; exit 1; }
  # tblout: 1 target, 3 query, 5 hmm_from, 6 hmm_to, 7 ali_from, 8 ali_to, 12 strand, 14 bits
  awk -v c="$c" -v m="$MODEL" -v minc="$MINC" -v bed="$O/calls.bed" -v cov="$O/coverage.tsv" '
    !/^#/ && $1 != "alphabet_guard" && $3 == m {
      n++; if ($6 - $5 + 1 < minc) {low++; next}
      s = ($12 == "+") ? $7 : $8; e = ($12 == "+") ? $8 : $7
      printf "%s\t%d\t%d\t%s:%d-%d:%s|%d|%d|%s\t0\t%s\n", c, s - 1, e, c, s, e, $12, $5, $6, $14, $12 >> bed; k++}
    END {printf "%s\t%d\t%d\t%d\n", c, n, k, low >> cov}' "$NW/chrom/$c.tbl"
done < "$CL"
sort -k1,1 -k2,2n "$O/calls.bed" > "$O/calls.sorted.bed"
twoBitToFa -bed="$O/calls.sorted.bed" "$TB" "$O/calls.fwd.fa"
# one line per record, reverse-complement the minus-strand calls, split the name back into columns
awk 'function rc(s,   r, i, b) {r = ""; for (i = length(s); i > 0; i--) {b = substr(s, i, 1)
                                     r = r (b == "A" ? "T" : b == "C" ? "G" : b == "G" ? "C" : b == "T" ? "A" : "N")}; return r}
     function out() {if (h == "") return; s = toupper(s); split(h, f, "|"); split(f[1], g, ":"); split(g[2], se, "-")
                     if (g[3] == "-") s = rc(s)
                     printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", f[1], g[1], se[1], se[2], g[3], f[2], f[3], f[4], s}
     /^>/ {out(); h = substr($0, 2); s = ""; next} {s = s $0} END {out()}' "$O/calls.fwd.fa" \
  | { printf 'id\tchrom\tstart\tend\tstrand\thmm_from\thmm_to\tbits\tseq\n'; cat; } > "$O/monomers.tsv"
nb=$(wc -l < "$O/calls.sorted.bed"); nm=$(($(wc -l < "$O/monomers.tsv") - 1))
[ "$nb" -eq "$nm" ] || { echo "cen113_monomer_table: $nb calls but $nm sequences" >&2; exit 1; }
dup=$(awk 'NR > 1 {print $1}' "$O/monomers.tsv" | sort | uniq -d | wc -l)
[ "$dup" -eq 0 ] || { echo "cen113_monomer_table: $dup duplicated ids" >&2; exit 1; }
awk -F'\t' 'NR > 1 {printf ">%s\n%s\n", $1, $9}' "$O/monomers.tsv" > "$O/monomers.fa"
/bin/rm -f "$O/calls.fwd.fa" "$O/calls.bed"
echo "cen113_monomer_table: $nm monomers"
cat "$O/coverage.tsv"
