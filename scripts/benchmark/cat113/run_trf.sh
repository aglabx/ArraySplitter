#!/usr/bin/env bash
# The TRF columns of Table 4: TRF (paper: v4.10.0) on each array of the 7-array panel, at two parameter sets.
#
#   run_trf.sh WORKDIR INPUTS TRF
#
#   satellome    2 5 7 80 10 50 2000   (the parameters satellome 1.6.x runs TRF with)
#   permissive   2 3 3 75 20 30 2000
# Both with -l 10 -d -h. Per set: WORKDIR/trf_<set>/<label>.fa and .dat, WORKDIR/trf_<set>.tsv (one row per TRF
# repeat) and WORKDIR/trf_summary_<set>.tsv (trf_period_summary.py: union coverage of periods 100-126, 200-252
# and 440-530 bp per array; Table 4 reports the 100-126-bp share). TRF 4.10.0 sometimes exits with no .dat
# (permissive set on the 1.02-Mb A3_483_control): such an array is named in WORKDIR/trf_<set>/FAILED.tsv, the
# summary prints it as no_trf_rows (never 0 %), and the script exits 3 at the end.
set -euo pipefail
W=${1:?workdir}; IN=${2:?inputs dir}; TRF=${3:?trf binary}
H=$(cd "$(dirname "$0")" && pwd)
[ -x "$TRF" ] || { echo "run_trf: $TRF not executable" >&2; exit 1; }
(cd "$IN" && md5sum -c "$H/data/inputs.md5") > /dev/null || { echo "run_trf: inputs differ from data/inputs.md5" >&2; exit 1; }
declare -A PARAMS=([satellome]="2 5 7 80 10 50 2000" [permissive]="2 3 3 75 20 30 2000")
for set in satellome permissive; do
  read -r -a P <<< "${PARAMS[$set]}"
  D=$W/trf_$set; mkdir -p "$D"
  printf 'label\ttrf_exit\n' > "$D/FAILED.tsv"
  printf 'label\tstart\tend\tlength\tperiod\tcopies\tconsensus_size\tpct_match\tpct_indel\tscore\tconsensus\n' > "$W/trf_$set.tsv"
  while IFS=$'\t' read -r acc s e label; do
    [ -n "$acc" ] || continue
    awk -v l="$label" '/^>/ {keep = (substr($0, 2) == l)} keep' "$IN/cat7.fa" > "$D/$label.fa"
    n=$(grep -v '^>' "$D/$label.fa" | tr -d '\n' | wc -c)
    [ "$n" -eq $((e - s)) ] || { echo "run_trf: $label is $n bp, expected $((e - s))" >&2; exit 1; }
    rc=0; (cd "$D" && "$TRF" "$label.fa" "${P[@]}" -l 10 -d -h > "$label.trf.log" 2>&1) || rc=$?   # TRF's exit code is not a status
    dat=$D/$label.fa.$(IFS=.; echo "${P[*]}").dat
    if [ ! -s "$dat" ]; then
      echo "run_trf: no TRF output for $label at $set (exit $rc; see $D/$label.trf.log)" >&2
      printf '%s\t%s\n' "$label" "$rc" >> "$D/FAILED.tsv"; continue
    fi
    awk -v l="$label" 'NF >= 14 && $1 ~ /^[0-9]+$/ && $2 ~ /^[0-9]+$/ {
        printf "%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", l, $1, $2, $2 - $1 + 1, $3, $4, $5, $6, $7, $8, $14}' "$dat" >> "$W/trf_$set.tsv"
  done < "$H/data/cat7_arrays.bed"
  python3 -I "$H/trf_period_summary.py" "$W/trf_$set.tsv" "$H/data/cat7_arrays.bed" 100-126,200-252,440-530 > "$W/trf_summary_$set.tsv"
done
echo "run_trf: $W/trf_summary_satellome.tsv and $W/trf_summary_permissive.tsv"
nf=$(awk 'FNR > 1' "$W"/trf_*/FAILED.tsv | wc -l)
[ "$nf" -eq 0 ] || { echo "run_trf: $nf array/parameter pair(s) without TRF output, named in $W/trf_*/FAILED.tsv" >&2; exit 3; }
