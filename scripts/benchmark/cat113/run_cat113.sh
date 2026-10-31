#!/usr/bin/env bash
# Every run behind the paper's §3.5 (Tables 4-5, Figure 2), in one workdir.
#
#   run_cat113.sh WORKDIR INPUTS ASPLIT ASPLIT_V182 FASTGA_DIR FASTAN CENTROANNO
#
# INPUTS      directory with cat7.fa, pat.fa, mat.fa (fetch_inputs.sh)
# ASPLIT      the arraysplitter binary under test (v1.9.0)
# ASPLIT_V182 the single-anchor baseline (v1.8.2, commit 433722a)
# FASTGA_DIR  FASTGA build with FAtoGDB and ONEview (FasTAN's input and output tools)
# FASTAN      FasTAN binary (paper: commit 963ac7f), CENTROANNO centroAnno binary (paper: commit 0b1f7dc)
#
# Lanes (4 threads for ArraySplitter and FasTAN, centroAnno defaults), each timed with GNU time -v:
#   asplit_new/{cat7,pat,mat}_{d5k,d32k}    ASPLIT, default cap and --max-period 32000
#   asplit_v182/{cat7,pat,mat}_{d5k,d32k}   ASPLIT_V182, the same
#   asplit_v182/cat7_d32k_rt0.3             ASPLIT_V182 --max-period 32000 --recursion-termination 0.3
#   fastan/{cat7,pat,mat}.1ano              FasTAN -ampv -T4
#   centroanno/{cat7,pat,mat}               centroAnno -x anno-asm
# Then the evaluator (eval/, Tables 4-5 -> eval_summary.txt) and the Figure 2 script (figures/).
# WORKDIR/binaries.md5 records every binary; WORKDIR/MANIFEST.md5 pins every deterministic output.
# Any failed lane stops the run.
set -euo pipefail
W=${1:?workdir}; IN=${2:?inputs dir}; NEW=${3:?arraysplitter}; OLD=${4:?arraysplitter v1.8.2}
FG=${5:?FASTGA dir}; FT=${6:?FasTAN}; CA=${7:?centroAnno}
H=$(cd "$(dirname "$0")" && pwd)
TIME=/usr/bin/time
for x in "$NEW" "$OLD" "$FG/FAtoGDB" "$FG/ONEview" "$FT" "$CA" "$TIME"; do
  [ -x "$x" ] || { echo "run_cat113: $x not executable" >&2; exit 1; }
done
(cd "$IN" && md5sum -c "$H/data/inputs.md5") > /dev/null || { echo "run_cat113: inputs differ from data/inputs.md5" >&2; exit 1; }
mkdir -p "$W"/{asplit_new,asplit_v182,fastan,centroanno,logs,figures}
md5sum "$NEW" "$OLD" "$FT" "$CA" "$FG/FAtoGDB" "$FG/ONEview" > "$W/binaries.md5"
{ "$NEW" --version; "$OLD" --version; } > "$W/binaries.version"

asplit() { # bin lane input tag [options]
  local bin=$1 lane=$2 in=$3 tag=$4; shift 4
  "$TIME" -v -o "$W/logs/${lane}_$tag.time" "$bin" -i "$IN/$in.fa" -o "$W/$lane/$tag" -t 4 "$@" \
    > "$W/logs/${lane}_$tag.log" 2>&1 || { echo "run_cat113: $lane $tag failed (logs/${lane}_$tag.log)" >&2; exit 1; }
  [ -s "$W/$lane/$tag.summary.tsv" ] || { echo "run_cat113: no $lane/$tag.summary.tsv" >&2; exit 1; }
}
for p in cat7 pat mat; do
  asplit "$NEW" asplit_new "$p" "${p}_d5k"
  asplit "$NEW" asplit_new "$p" "${p}_d32k" --max-period 32000
  asplit "$OLD" asplit_v182 "$p" "${p}_d5k"
  asplit "$OLD" asplit_v182 "$p" "${p}_d32k" --max-period 32000
done
asplit "$OLD" asplit_v182 cat7 cat7_d32k_rt0.3 --max-period 32000 --recursion-termination 0.3

for p in cat7 pat mat; do
  cp "$IN/$p.fa" "$W/fastan/$p.fa"
  "$FG/FAtoGDB" "$W/fastan/$p.fa" > "$W/logs/fatogdb_$p.log" 2>&1 || { echo "run_cat113: FAtoGDB $p failed" >&2; exit 1; }
  "$TIME" -v -o "$W/logs/fastan_$p.time" "$FT" -ampv -T4 -o"$W/fastan/$p.1ano" "$W/fastan/$p.1gdb" \
    > "$W/logs/fastan_$p.log" 2>&1 || { echo "run_cat113: FasTAN $p failed" >&2; exit 1; }
  "$TIME" -v -o "$W/logs/centro_$p.time" "$CA" "$IN/$p.fa" -o "$W/centroanno/$p" -x anno-asm \
    > "$W/logs/centro_$p.log" 2>&1 || { echo "run_cat113: centroAnno $p failed" >&2; exit 1; }
done

(cd "$H/eval" && cargo build --release --quiet)
md5sum "$H/eval/target/release/cat113_constellation_eval" >> "$W/binaries.md5"
"$H/eval/target/release/cat113_constellation_eval" "$W" "$FG/ONEview" > "$W/eval_summary.txt"
python3 "$H/plot_cat113_fields.py" "$W/asplit_new" "$W/figures" > "$W/logs/figures.log" 2>&1 \
  || { echo "run_cat113: figure script failed (logs/figures.log)" >&2; exit 1; }

(cd "$W" && find asplit_new asplit_v182 centroanno -type f | sort | xargs md5sum; md5sum fastan/*.1ano eval_summary.txt) > "$W/MANIFEST.md5"
echo "run_cat113: done, $(wc -l < "$W/MANIFEST.md5") files in $W/MANIFEST.md5"
