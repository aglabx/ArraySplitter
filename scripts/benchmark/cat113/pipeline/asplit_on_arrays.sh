#!/usr/bin/env bash
# ArraySplitter (a given build) on the arrays cut by trf_on_arrays.sh (2026-10-30, PI «А если на аррейсплиттер
# последней версии?»: the cat catalogue was decomposed by 1.7.4 inside satellome 1.17.2; the latest tag is v1.8.2).
#
#   asplit_on_arrays.sh FA_DIR LABELS BIN WORKDIR [ARRAYSPLITTER OPTIONS...]
#
# FA_DIR/<label>.fa for each label in LABELS (one per line; trf_on_arrays.sh writes them). They are joined into
# WORKDIR/arrays.fa with the label as the header, BIN runs once with -t 4 and the options given, and stdout gets one
# row per array from its summary.tsv: label, length, hor_period, hor_n_monomers, mono_period, mono_n_monomers,
# mono_mean_ed_tmpl, mono_consensus. An array missing from the summary is named and the script exits 3.
set -euo pipefail
FD=${1:?fa dir}; LIST=${2:?labels}; BIN=${3:?arraysplitter binary}; W=${4:?workdir}; shift 4
[ -x "$BIN" ] || { echo "asplit_on_arrays: $BIN not executable" >&2; exit 1; }
mkdir -p "$W"
: > "$W/arrays.fa"
while read -r label; do
  [ -n "$label" ] || continue
  [ -s "$FD/$label.fa" ] || { echo "asplit_on_arrays: $FD/$label.fa missing" >&2; exit 1; }
  { echo ">$label"; grep -v '^>' "$FD/$label.fa"; } >> "$W/arrays.fa"
done < "$LIST"
"$BIN" --version > "$W/binary.version"
md5sum "$BIN" > "$W/binary.md5"
"$BIN" -i "$W/arrays.fa" -o "$W/out" -t 4 "$@" > "$W/run.log" 2>&1 || { echo "asplit_on_arrays: $BIN failed (see $W/run.log)" >&2; exit 1; }
[ -s "$W/out.summary.tsv" ] || { echo "asplit_on_arrays: no $W/out.summary.tsv" >&2; exit 1; }
printf 'label\tlength\thor_period\thor_n_monomers\tmono_period\tmono_n_monomers\tmono_mean_ed_tmpl\tmono_consensus\n'
awk -F'\t' 'NR == FNR {want[$1] = 1; order[++n] = $1; next}
    FNR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
    ($1 in want) {row[$1] = $1 "\t" $c["array_length"] "\t" $c["hor_period"] "\t" $c["hor_n_monomers"] "\t" \
        $c["mono_period"] "\t" $c["mono_n_monomers"] "\t" $c["mono_mean_ed_tmpl"] "\t" $c["mono_consensus"]}
    END {miss = 0; for (i = 1; i <= n; i++) {if (order[i] == "") continue
           if (order[i] in row) print row[order[i]]; else {print "asplit_on_arrays: " order[i] " not in summary" > "/dev/stderr"; miss++}}
         if (miss) exit 3}' "$LIST" "$W/out.summary.tsv"
