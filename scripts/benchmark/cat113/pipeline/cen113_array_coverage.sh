#!/usr/bin/env bash
# How much of each cat 113-bp array the monomer calls cover, and where they miss (2026-10-30,
# EXP-cat-113-monomer-letters step 2).
#
#   cen113_array_coverage.sh NHMMER_WORKDIR MONOMERS_TSV OTHER_MODEL MIN_CALLS OUTDIR
#
# Arrays: hits of both profiles (NHMMER_WORKDIR/chrom/*.hits.tsv) merged when <= 5 kb apart, as in
# cen113_profile_search.sh. Calls: MONOMERS_TSV of cen113_monomer_table.sh.
# OUTDIR/groups.tsv: every group, with span, calls, bp covered by calls, fraction covered, and kept = yes when it
#   has >= MIN_CALLS calls (groups with fewer stay in the table, marked no).
# OUTDIR/arrays.bed: the kept groups (array id = chrom:start-end, 1-based).
# OUTDIR/exceptions.tsv: gaps >= 150 bp between calls inside a kept array, with the bp of OTHER_MODEL hits in them.
set -euo pipefail
NW=${1:?nhmmer workdir}; M=${2:?monomers.tsv}; OM=${3:?other model}; MINC=${4:?min calls}; O=${5:?outdir}
command -v bedtools > /dev/null || { echo "cen113_array_coverage: bedtools not on PATH" >&2; exit 1; }
mkdir -p "$O"
ls "$NW"/chrom/*.hits.tsv > /dev/null
cat "$NW"/chrom/*.hits.tsv | awk -F'\t' -v OFS='\t' '{print $1, $2 - 1, $3, $4}' | sort -k1,1 -k2,2n > "$O/hits.bed"
awk -F'\t' -v OFS='\t' -v m="$OM" '$4 == m' "$O/hits.bed" > "$O/other_hits.bed"
[ -s "$O/other_hits.bed" ] || { echo "cen113_array_coverage: no hits of $OM" >&2; exit 1; }
awk -F'\t' -v OFS='\t' 'NR > 1 {print $2, $3 - 1, $4, $1}' "$M" | sort -k1,1 -k2,2n > "$O/calls.bed"
bedtools merge -i "$O/hits.bed" -d 5000 -c 4 -o count > "$O/groups.bed"
# coverage columns appended: n calls overlapping, bp covered, group length, fraction
bedtools coverage -a "$O/groups.bed" -b "$O/calls.bed" \
  | awk -F'\t' -v OFS='\t' -v minc="$MINC" '
      BEGIN {print "array", "chrom", "start", "end", "span", "hits", "calls", "bp_called", "fraction_called", "kept"}
      {print $1 ":" $2 + 1 "-" $3, $1, $2 + 1, $3, $3 - $2, $4, $5, $6, sprintf("%.3f", $8), ($5 >= minc ? "yes" : "no")}' \
  > "$O/groups.tsv"
awk -F'\t' -v OFS='\t' 'NR > 1 && $10 == "yes" {print $2, $3 - 1, $4, $1}' "$O/groups.tsv" > "$O/arrays.bed"
bedtools subtract -a "$O/arrays.bed" -b "$O/calls.bed" | awk -F'\t' -v OFS='\t' '$3 - $2 >= 150' > "$O/gaps.bed"
{ printf 'array\tchrom\tgap_start\tgap_end\tgap_length\tother_model_bp\tother_model_fraction\n'
  if [ -s "$O/gaps.bed" ]; then
    bedtools coverage -a "$O/gaps.bed" -b "$O/other_hits.bed" \
      | awk -F'\t' -v OFS='\t' '{print $4, $1, $2 + 1, $3, $3 - $2, $6, sprintf("%.3f", $8)}'
  fi; } > "$O/exceptions.tsv"
echo "cen113_array_coverage: $(($(wc -l < "$O/groups.tsv") - 1)) groups, $(wc -l < "$O/arrays.bed") kept, $(($(wc -l < "$O/exceptions.tsv") - 1)) gaps >= 150 bp"
