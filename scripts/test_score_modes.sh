#!/bin/bash
# Test different SCORE_MODE formulas for cut-weighted scoring
# Usage: ./test_score_modes.sh

set -e

WORKDIR="/mnt/data/claude/arraysplitter_rust"
INPUT="/mnt/data/claude/phase1_alpha_extraction/all_alpha_satellite.fasta"
RUST_SRC="$WORKDIR/src/anchor_graph.rs"
RESULTS_DIR="$WORKDIR/score_mode_tests"
ANALYZE="$WORKDIR/scripts/analyze_decomposition.py"

# Score modes to test
MODES="0 1 2 3 4"
MODE_NAMES=("mean_ed" "mean_ed/sqrt(cuts)" "mean_ed/log2(cuts)" "mean_ed/cuts" "mean_ed*(1+10/sqrt(cuts))")

mkdir -p "$RESULTS_DIR"

echo "=== Testing different SCORE_MODE formulas ==="
echo "Input: $INPUT"
echo "Results: $RESULTS_DIR"
echo ""
echo "Modes:"
echo "  0: mean_ed (baseline)"
echo "  1: mean_ed / sqrt(cuts)"
echo "  2: mean_ed / log2(cuts)"
echo "  3: mean_ed / cuts"
echo "  4: mean_ed * (1 + 10/sqrt(cuts))"
echo ""

# Save original file
cp "$RUST_SRC" "$RUST_SRC.backup"

for MODE in $MODES; do
    echo "=========================================="
    echo "Testing SCORE_MODE = $MODE (${MODE_NAMES[$MODE]})"
    echo "=========================================="

    # Update SCORE_MODE in Rust code
    sed -i "s/const SCORE_MODE: usize = [0-9]*;/const SCORE_MODE: usize = $MODE;/" "$RUST_SRC"

    # Verify change
    grep "const SCORE_MODE" "$RUST_SRC"

    # Rebuild
    echo "Building..."
    cd "$WORKDIR"
    ~/.cargo/bin/cargo build --release 2>&1 | tail -2

    # Run analysis
    OUTPUT_PREFIX="$RESULTS_DIR/mode_${MODE}"
    echo "Running analysis -> $OUTPUT_PREFIX"

    START_TIME=$(date +%s)
    ./target/release/arraysplitter_rs -i "$INPUT" -o "$OUTPUT_PREFIX" 2>&1 | tail -5
    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))

    echo "Elapsed: ${ELAPSED}s"

    # Analyze results
    echo "Analyzing results..."
    python3 "$ANALYZE" \
        -m "${OUTPUT_PREFIX}.monomers.tsv" \
        -f "${OUTPUT_PREFIX}.decomposed.fasta" \
        -o "${OUTPUT_PREFIX}.analysis.tsv" 2>&1 | grep -A10 "=== Summary ==="

    echo ""
done

# Restore original file
cp "$RUST_SRC.backup" "$RUST_SRC"
echo "Restored original anchor_graph.rs"

# Summary table
echo ""
echo "=========================================="
echo "SUMMARY: SCORE_MODE comparison"
echo "=========================================="
printf "%-30s | %12s | %12s | %10s | %10s | %12s\n" "Mode" "EXACT" "HALF/DBL" "3RD/3PL" "MISMATCH" "Correct*"
echo "-------------------------------|--------------|--------------|------------|------------|-------------"

for MODE in $MODES; do
    OUTPUT_PREFIX="$RESULTS_DIR/mode_${MODE}"
    if [ -f "${OUTPUT_PREFIX}.analysis.tsv" ]; then
        TOTAL=$(tail -n +2 "${OUTPUT_PREFIX}.analysis.tsv" | wc -l)
        EXACT=$(tail -n +2 "${OUTPUT_PREFIX}.analysis.tsv" | awk -F'\t' '$12=="EXACT"' | wc -l)
        CLOSE=$(tail -n +2 "${OUTPUT_PREFIX}.analysis.tsv" | awk -F'\t' '$12=="CLOSE"' | wc -l)
        HALF=$(tail -n +2 "${OUTPUT_PREFIX}.analysis.tsv" | awk -F'\t' '$12=="HALF/DOUBLE"' | wc -l)
        THIRD=$(tail -n +2 "${OUTPUT_PREFIX}.analysis.tsv" | awk -F'\t' '$12=="THIRD/TRIPLE"' | wc -l)
        MISMATCH=$(tail -n +2 "${OUTPUT_PREFIX}.analysis.tsv" | awk -F'\t' '$12=="MISMATCH"' | wc -l)

        # Count correct HALF/DOUBLE (ratio < 0.6 = found monomers in HOR)
        HALF_CORRECT=$(tail -n +2 "${OUTPUT_PREFIX}.analysis.tsv" | awk -F'\t' '$12=="HALF/DOUBLE" && $11+0 < 0.6' | wc -l)

        CORRECT=$((EXACT + CLOSE + HALF_CORRECT))
        if [ "$TOTAL" -gt 0 ]; then
            EXACT_PCT=$(echo "scale=1; $EXACT * 100 / $TOTAL" | bc)
            MISMATCH_PCT=$(echo "scale=1; $MISMATCH * 100 / $TOTAL" | bc)
            CORRECT_PCT=$(echo "scale=1; $CORRECT * 100 / $TOTAL" | bc)
        fi

        printf "%-30s | %4d (%5s%%) | %4d         | %4d       | %4d (%5s%%) | %4d (%5s%%)\n" \
            "${MODE_NAMES[$MODE]}" \
            "$EXACT" "$EXACT_PCT" "$HALF" "$THIRD" "$MISMATCH" "$MISMATCH_PCT" "$CORRECT" "$CORRECT_PCT"
    fi
done

echo ""
echo "* Correct = EXACT + CLOSE + HALF/DOUBLE(ratio<0.6, i.e. found monomers in HOR)"
echo ""
echo "Done! Results in $RESULTS_DIR"
