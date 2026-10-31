#!/usr/bin/env bash
# Second pass of the cat 113-bp monomer calls with a profile per chromosome (2026-10-30,
# EXP-cat-113-monomer-letters step 2: the first pass misses one or two HOR positions per copy, 150-220-bp gaps
# that weak hits show to be 113-bp monomers).
#
#   cen113_second_pass.sh MONOMERS_TSV ARRAYS_BED TWOBIT OUTDIR [THREADS]
#
# Per chromosome with >= 10 distinct full-frame calls (hmm_from <= 3, hmm_to >= 110 of the E1 profile): up to 200
# of them (every n-th of the distinct sequences sorted, deterministic), MAFFT --auto, hmmbuild --dna; nhmmer
# (E <= 1e-5) over that chromosome's arrays (ARRAYS_BED, +- 200 bp). A hit covering >= 80% of the model's columns
# is a call (pass 1 kept >= 90 of 113 columns, the same share).
# OUTDIR/monomers.tsv: pass-1 calls plus pass-2 calls that overlap no pass-1 call by >= 50% of their length;
# column pass = 1 or 2, hmm columns are of the model that made the call. OUTDIR/added.tsv: per chromosome, pass-2
# hits, calls, added. Chromosomes without a profile are listed there with reason, not skipped silently.
set -euo pipefail
M=${1:?monomers.tsv}; AB=${2:?arrays.bed}; TB=${3:?2bit}; O=${4:?outdir}; T=${5:-8}
for t in mafft hmmbuild nhmmer twoBitToFa bedtools; do command -v $t > /dev/null || { echo "cen113_second_pass: $t not on PATH" >&2; exit 1; }; done
mkdir -p "$O/chrom"
printf 'chrom\tdistinct_full_frame\tprofile\tpass2_hits\tpass2_calls\tadded\n' > "$O/added.tsv"
: > "$O/pass2.tsv"
for c in $(cut -f1 "$AB" | sort -u); do
  W="$O/chrom/$c"; mkdir -p "$W"
  awk -F'\t' -v c="$c" 'NR > 1 && $2 == c && $6 <= 3 && $7 >= 110 {print $9}' "$M" | sort -u > "$W/full.txt"
  n=$(wc -l < "$W/full.txt")
  if [ "$n" -lt 10 ]; then printf '%s\t%d\tnone: < 10 distinct full-frame calls\t-\t-\t0\n' "$c" "$n" >> "$O/added.tsv"; continue; fi
  awk -v n="$n" -v k=200 'BEGIN {step = (n > k) ? n / k : 1} {i++; if (i >= nxt) {printf ">s%d\n%s\n", i, $0; nxt += step}}' "$W/full.txt" \
    | awk '/^>/ {m++} m <= 200' > "$W/sample.fa"
  mafft --auto --thread "$T" --quiet "$W/sample.fa" > "$W/aln.fa"
  hmmbuild --dna -n "cat113_$c" "$W/model.hmm" "$W/aln.fa" > "$W/hmmbuild.log"
  L=$(awk '$1 == "LENG" {print $2; exit}' "$W/model.hmm")
  awk -F'\t' -v OFS='\t' -v c="$c" '$1 == c {s = $2 - 200; if (s < 0) s = 0; print $1, s, $3 + 200, $1 "_" s}' "$AB" > "$W/regions.bed"
  { echo ">alphabet_guard"; awk 'BEGIN {srand(1); for (i = 0; i < 1000; i++) printf "%s", substr("ACGT", int(rand() * 4) + 1, 1); print ""}'
    twoBitToFa -bed="$W/regions.bed" "$TB" stdout; } > "$W/regions.fa"
  [ "$(grep -c '^>' "$W/regions.fa")" -eq $(($(wc -l < "$W/regions.bed") + 1)) ] || { echo "cen113_second_pass: region cut failed on $c" >&2; exit 1; }
  nhmmer --dna --cpu "$T" --noali -E 1e-5 --tblout "$W/pass2.tbl" "$W/model.hmm" "$W/regions.fa" > "$W/nhmmer.log" \
    || { echo "cen113_second_pass: nhmmer failed on $c" >&2; exit 1; }
  # region name is <chrom>_<0-based offset>; tblout 5 hmm_from, 6 hmm_to, 7/8 ali, 12 strand, 14 bits
  awk -v OFS='\t' -v L="$L" -v c="$c" '!/^#/ && $1 != "alphabet_guard" {h++; if ($6 - $5 + 1 < 0.8 * L) next
      off = substr($1, length(c) + 2) + 0; s = ($12 == "+") ? $7 : $8; e = ($12 == "+") ? $8 : $7
      print c, off + s - 1, off + e, $12, $5, $6, $14}
      END {print h + 0 > "/dev/stderr"}' "$W/pass2.tbl" 2> "$W/hits.n" | sort -k2,2n > "$W/pass2.bed"
  # keep pass-2 calls that overlap no pass-1 call by >= 50% of their own length
  awk -F'\t' -v OFS='\t' -v c="$c" 'NR > 1 && $2 == c {print $2, $3 - 1, $4}' "$M" | sort -k2,2n > "$W/pass1.bed"
  if [ -s "$W/pass2.bed" ]; then
    bedtools intersect -a "$W/pass2.bed" -b "$W/pass1.bed" -wao \
      | awk -F'\t' -v OFS='\t' '{k = $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7; ov[k] += $NF; if (!(k in seen)) {seen[k] = 1; ord[++n] = k}}
          END {for (i = 1; i <= n; i++) {split(ord[i], f, OFS); if (ov[ord[i]] < 0.5 * (f[3] - f[2])) print ord[i]}}' > "$W/added.bed"
  else : > "$W/added.bed"; fi
  cat "$W/added.bed" >> "$O/pass2.tsv"
  printf '%s\t%d\tcat113_%s (LENG %s)\t%s\t%d\t%d\n' "$c" "$n" "$c" "$L" "$(cat "$W/hits.n")" "$(wc -l < "$W/pass2.bed")" "$(wc -l < "$W/added.bed")" >> "$O/added.tsv"
done
# sequences of the added calls, oriented to their model, then the merged table
if [ -s "$O/pass2.tsv" ]; then
  awk -F'\t' -v OFS='\t' '{print $1, $2, $3, $1 ":" $2 + 1 "-" $3 ":" $4 "|" $5 "|" $6 "|" $7}' "$O/pass2.tsv" | sort -k1,1 -k2,2n > "$O/pass2.named.bed"
  twoBitToFa -bed="$O/pass2.named.bed" "$TB" "$O/pass2.fwd.fa"
  [ "$(grep -c '^>' "$O/pass2.fwd.fa")" -eq "$(wc -l < "$O/pass2.named.bed")" ] || { echo "cen113_second_pass: pass-2 sequence count mismatch" >&2; exit 1; }
fi
{ printf 'id\tchrom\tstart\tend\tstrand\thmm_from\thmm_to\tbits\tseq\tpass\n'
  awk -F'\t' -v OFS='\t' 'NR > 1 {print $0, 1}' "$M"
  if [ -s "$O/pass2.tsv" ]; then
    awk 'function rc(s,   r, i, b) {r = ""; for (i = length(s); i > 0; i--) {b = substr(s, i, 1)
                                         r = r (b == "A" ? "T" : b == "C" ? "G" : b == "G" ? "C" : b == "T" ? "A" : "N")}; return r}
         function out() {if (h == "") return; s = toupper(s); split(h, f, "|"); split(f[1], g, ":"); split(g[2], se, "-")
                         if (g[3] == "-") s = rc(s)
                         printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t2\n", f[1], g[1], se[1], se[2], g[3], f[2], f[3], f[4], s}
         /^>/ {out(); h = substr($0, 2); s = ""; next} {s = s $0} END {out()}' "$O/pass2.fwd.fa"
  fi; } > "$O/monomers.tsv"
dup=$(awk 'NR > 1 {print $1}' "$O/monomers.tsv" | sort | uniq -d | wc -l)
[ "$dup" -eq 0 ] || { echo "cen113_second_pass: $dup duplicated ids" >&2; exit 1; }
awk -F'\t' 'NR > 1 {printf ">%s\n%s\n", $1, $9}' "$O/monomers.tsv" > "$O/monomers.fa"
column -t -s $'\t' "$O/added.tsv"
echo "cen113_second_pass: $(($(wc -l < "$O/monomers.tsv") - 1)) monomers ($(awk -F'\t' 'NR > 1 && $10 == 2' "$O/monomers.tsv" | wc -l) added)"
