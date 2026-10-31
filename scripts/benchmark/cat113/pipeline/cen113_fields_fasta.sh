#!/usr/bin/env bash
# The cat 113-bp satellite fields (arrays) as FASTA (2026-10-31, PI «а ты можешь мне дать поля 113 как фасту»).
#
#   cen113_fields_fasta.sh GROUPS_TSV MONOMERS_TSV TWOBIT ASSEMBLY_REPORT OUT_FA
#
# GROUPS_TSV: groups.tsv of cen113_array_coverage.sh; the kept arrays (>= its MIN_CALLS calls) are written.
# Sequence is the assembly's forward strand. Header:
#   >NAME|ACCESSION:START-END|LENGTHbp|calls=N|calls_minus=M
# NAME = chromosome name from the NCBI assembly report; calls_minus = how many of the array's monomer calls lie on the
# minus strand (the monomer frame of the E1 profile reads them reverse-complemented). OUT_FA.tsv: the same as a table.
# A chromosome missing from the report, or a record count that differs from the kept arrays, stops the run.
set -euo pipefail
G=${1:?groups.tsv}; M=${2:?monomers.tsv}; TB=${3:?2bit}; R=${4:?assembly report}; O=${5:?out fasta}
command -v twoBitToFa > /dev/null || { echo "cen113_fields_fasta: twoBitToFa not on PATH" >&2; exit 1; }
grep -v '^#' "$R" | awk -F'\t' '$2 == "assembled-molecule" {print $5 "\t" $1}' > "$O.names"
# kept arrays with their minus-strand call counts
awk -F'\t' -v OFS='\t' '
  FILENAME == ARGV[1] {nm[$1] = $2; next}
  FILENAME == ARGV[2] {if (FNR == 1) next; n++; c[n] = $2; s[n] = $3; e[n] = $4; st[n] = $5; next}
  FNR == 1 {next}
  $10 == "yes" {if (!($2 in nm)) {print "cen113_fields_fasta: " $2 " not in the assembly report" > "/dev/stderr"; exit 2}
                minus = 0; calls = 0
                for (i = 1; i <= n; i++) if (c[i] == $2 && s[i] >= $3 && e[i] <= $4) {calls++; if (st[i] == "-") minus++}
                print $2, $3 - 1, $4, nm[$2] "|" $2 ":" $3 "-" $4 "|" $5 "bp|calls=" calls "|calls_minus=" minus}' \
  "$O.names" "$M" "$G" > "$O.bed"
twoBitToFa -bed="$O.bed" "$TB" "$O"
k=$(wc -l < "$O.bed"); f=$(grep -c '^>' "$O")
[ "$k" -eq "$f" ] || { echo "cen113_fields_fasta: $k arrays but $f records" >&2; exit 1; }
{ printf 'name\taccession\tstart\tend\tlength\tcalls\tcalls_minus\n'
  awk -F'\t' -v OFS='\t' '{split($4, h, "|"); split(h[2], a, ":"); split(a[2], se, "-"); sub("bp", "", h[3]); sub("calls=", "", h[4]); sub("calls_minus=", "", h[5])
                           print h[1], a[1], se[1], se[2], h[3], h[4], h[5]}' "$O.bed"; } > "$O.tsv"
/bin/rm -f "$O.names" "$O.bed"
echo "cen113_fields_fasta: $f fields, $(grep -v '^>' "$O" | tr -d '\n' | wc -c) bp -> $O"
