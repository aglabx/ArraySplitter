#!/usr/bin/env bash
# TRF on chosen tandem arrays of a tree 2bit, to see which period TRF finds where FasTAN reported a long one
# (2026-10-30, PI «давай TRF запустим на этих массивах»: the cat arrays carrying the 113-bp centromeric unit, which
# FasTAN reported with periods of 2.6-15 kb; CLM-cat-locus-is-the-pericentromeric-483-satellite).
#
#   trf_on_arrays.sh ARRAYS.bed TWOBIT WORKDIR
#
# ARRAYS.bed: seq, start (0-based), end, label. Each array is cut with twoBitToFa and run through TRF with satellome
# 1.6.x's own parameters (src/satellome/constants.py: 2 5 7 80 10 50 2000; -l 10 allows arrays up to 10 Mb),
# or TRF_PARAMS.
# WORKDIR/<label>.fa, <label>.fa.<params>.dat, and on stdout one row per TRF repeat:
# label, start, end (in the array), length, period, copies, consensus size, % match, % indel, score, consensus.
# A failed cut or a wrong-length sequence stops the run. TRF 4.10.0 with permissive parameters sometimes exits 255
# with no .dat (cat, 2026-10-30: 2 3 3 75 20 30 on some 100-kb windows); such an array is named in WORKDIR/FAILED.tsv
# with TRF's exit code, the run goes on, and the script exits 3 at the end if any failed.
set -euo pipefail
BED=${1:?arrays bed}; TB=${2:?2bit}; W=${3:?workdir}
read -r -a P <<< "${TRF_PARAMS:-2 5 7 80 10 50 2000}"   # TRF_PARAMS: another set, e.g. the most permissive TRF documents
[ "${#P[@]}" -eq 7 ] || { echo "trf_on_arrays: TRF_PARAMS must hold 7 numbers" >&2; exit 1; }
command -v twoBitToFa > /dev/null || { echo "trf_on_arrays: twoBitToFa not on PATH" >&2; exit 1; }
TRF=${TRF:-trf}
command -v "$TRF" > /dev/null || { echo "trf_on_arrays: $TRF not on PATH" >&2; exit 1; }
mkdir -p "$W"
printf 'label\ttrf_exit\n' > "$W/FAILED.tsv"
nfail=0
printf 'label\tstart\tend\tlength\tperiod\tcopies\tconsensus_size\tpct_match\tpct_indel\tscore\tconsensus\n'
while IFS=$'\t' read -r seq start end label; do
  [ -n "$seq" ] || continue
  fa=$W/$label.fa
  twoBitToFa -seq="$seq" -start="$start" -end="$end" "$TB" "$fa" || { echo "trf_on_arrays: cut failed for $label" >&2; exit 1; }
  n=$(grep -v '^>' "$fa" | tr -d '\n' | wc -c)
  [ "$n" -eq $((end - start)) ] || { echo "trf_on_arrays: $label is $n bp, expected $((end - start))" >&2; exit 1; }
  rc=0; (cd "$W" && "$TRF" "$label.fa" "${P[@]}" -l 10 -d -h > "$label.trf.log" 2>&1) || rc=$?   # TRF's exit code is not a status
  dat=$W/$label.fa.$(IFS=.; echo "${P[*]}").dat
  if [ ! -s "$dat" ]; then
    echo "trf_on_arrays: no TRF output for $label (exit $rc; see $W/$label.trf.log)" >&2
    printf '%s\t%s\n' "$label" "$rc" >> "$W/FAILED.tsv"; nfail=$((nfail + 1)); continue
  fi
  awk -v l="$label" 'NF >= 14 && $1 ~ /^[0-9]+$/ && $2 ~ /^[0-9]+$/ {
      printf "%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", l, $1, $2, $2 - $1 + 1, $3, $4, $5, $6, $7, $8, $14}' "$dat"
done < "$BED"
[ "$nfail" -eq 0 ] || { echo "trf_on_arrays: $nfail array(s) without TRF output, listed in $W/FAILED.tsv" >&2; exit 3; }
