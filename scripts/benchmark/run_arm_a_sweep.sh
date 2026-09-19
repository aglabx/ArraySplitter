#!/usr/bin/env bash
# Arm A parameter sensitivity sweep for ArraySplitter --method autocorr.
#
# Sweeps four CLI levers introduced in ArraySplitter commit 911c389 on the
# 12-fixture sim panel + the 60-Mb zfinch panel. Each parameter is swept
# independently (others held at default).
#
# Designed to run from ~/asplit_bench on aglab0; rsync this script there.
#
# Defaults reproduced (must produce byte-identical output to no-flag run):
#   --excess-floor 0.05  --recursion-termination 0.5
#   --multiplet-factor 1.5  --period-finder mode-dependent

set -eu

OUT=ablation_arm_a
mkdir -p "$OUT"
ASPLIT=~/asplit_bench/arraysplitter_repo/src/rust/arraysplitter/target/release/arraysplitter

SIM=~/asplit_bench/sim_fixtures.fasta
ZF=~/asplit_bench/zf.fa

# Helper: one (panel × params) cell.  $1=panel-label, $2=fasta, $3..=flags
run_cell() {
    local label="$1"; shift
    local fa="$1"; shift
    local tag="$1"; shift
    local out="$OUT/${label}_${tag}"
    if [ -f "${out}.summary.tsv" ]; then
        echo "skip ${label} ${tag} (already)"
        return
    fi
    echo "--- ${label} ${tag} ---"
    /usr/bin/time -f 'real %es' \
        "$ASPLIT" -i "$fa" -o "$out" --max-period 100000 -t 4 --method autocorr "$@" \
        2>&1 | tail -2
}

echo "##### sim sweep #####"
# Default (sanity)
run_cell sim "$SIM" default
# excess-floor: 0.02, 0.05 (default), 0.10, 0.20
for v in 0.02 0.10 0.20; do
    run_cell sim "$SIM" "excess${v}" --excess-floor "$v"
done
# recursion-termination: 0.3, 0.5 (default), 0.7
for v in 0.3 0.7; do
    run_cell sim "$SIM" "recurse${v}" --recursion-termination "$v"
done
# multiplet-factor: 1.1, 1.3, 1.5 (default), 2.0
for v in 1.1 1.3 2.0; do
    run_cell sim "$SIM" "mult${v}" --multiplet-factor "$v"
done
# period-finder: raw, refined, mode-dependent (default)
for v in raw refined; do
    run_cell sim "$SIM" "finder-${v}" --period-finder "$v"
done

echo "##### zfinch sweep #####"
run_cell zf "$ZF" default
for v in 0.02 0.10 0.20; do
    run_cell zf "$ZF" "excess${v}" --excess-floor "$v"
done
for v in 0.3 0.7; do
    run_cell zf "$ZF" "recurse${v}" --recursion-termination "$v"
done
for v in 1.1 1.3 2.0; do
    run_cell zf "$ZF" "mult${v}" --multiplet-factor "$v"
done
for v in raw refined; do
    run_cell zf "$ZF" "finder-${v}" --period-finder "$v"
done

echo
echo "##### sweep done #####"
ls -lh "$OUT"/*.summary.tsv | head -50
