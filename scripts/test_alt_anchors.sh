#!/bin/bash
# Test redundant anchors improvement
# Compare current (with alt anchors) vs previous Mode 3 (without alt anchors)
# Usage: ./test_alt_anchors.sh

set -e

WORKDIR="/mnt/data/claude/arraysplitter_rust"
INPUT="/mnt/data/claude/phase1_alpha_extraction/all_alpha_satellite.fasta"
RESULTS_DIR="$WORKDIR/alt_anchor_tests"
ANALYZE="$WORKDIR/scripts/analyze_decomposition.py"
PREV_RESULTS="$WORKDIR/score_mode_tests/mode_3"

mkdir -p "$RESULTS_DIR"

echo "=== Testing Redundant Anchors ==="
echo "Input: $INPUT"
echo "Results: $RESULTS_DIR"
echo ""

# Build
echo "Building..."
cd "$WORKDIR"
~/.cargo/bin/cargo build --release 2>&1 | tail -2

# Run with alt anchors (current code)
OUTPUT_PREFIX="$RESULTS_DIR/with_alt_anchors"
echo "Running with alternative anchors -> $OUTPUT_PREFIX"

START_TIME=$(date +%s)
./target/release/arraysplitter_rs -i "$INPUT" -o "$OUTPUT_PREFIX" -v 2>"$RESULTS_DIR/verbose.log"
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
echo "Elapsed: ${ELAPSED}s"

# Count how many times alt anchors were used
ALT_COUNT=$(grep -c "Alternative anchors" "$RESULTS_DIR/verbose.log" || echo "0")
ED_SPLITS=$(grep -c "ED-split" "$RESULTS_DIR/verbose.log" || echo "0")
echo "Arrays with alternative anchors: $ALT_COUNT"
echo "ED-splits performed: $ED_SPLITS"

# Analyze results
echo ""
echo "Analyzing results..."
python3 "$ANALYZE" \
    -m "${OUTPUT_PREFIX}.monomers.tsv" \
    -f "${OUTPUT_PREFIX}.decomposed.fasta" \
    -o "${OUTPUT_PREFIX}.analysis.tsv" 2>&1 | grep -A10 "=== Summary ==="

# Compare with previous Mode 3 results
echo ""
echo "=========================================="
echo "COMPARISON: With Alt Anchors vs Without"
echo "=========================================="
printf "%-25s | %12s | %12s | %10s | %10s | %12s\n" "Version" "EXACT" "HALF/DBL" "3RD/3PL" "MISMATCH" "Correct*"
echo "--------------------------|--------------|--------------|------------|------------|-------------"

for VERSION in "prev_mode3:$PREV_RESULTS" "with_alt:$OUTPUT_PREFIX"; do
    NAME="${VERSION%%:*}"
    PREFIX="${VERSION##*:}"
    ANALYSIS="${PREFIX}.analysis.tsv"

    if [ -f "$ANALYSIS" ]; then
        TOTAL=$(tail -n +2 "$ANALYSIS" | wc -l)
        EXACT=$(tail -n +2 "$ANALYSIS" | awk -F'\t' '$12=="EXACT"' | wc -l)
        CLOSE=$(tail -n +2 "$ANALYSIS" | awk -F'\t' '$12=="CLOSE"' | wc -l)
        HALF=$(tail -n +2 "$ANALYSIS" | awk -F'\t' '$12=="HALF/DOUBLE"' | wc -l)
        THIRD=$(tail -n +2 "$ANALYSIS" | awk -F'\t' '$12=="THIRD/TRIPLE"' | wc -l)
        MISMATCH=$(tail -n +2 "$ANALYSIS" | awk -F'\t' '$12=="MISMATCH"' | wc -l)
        HALF_CORRECT=$(tail -n +2 "$ANALYSIS" | awk -F'\t' '$12=="HALF/DOUBLE" && $11+0 < 0.6' | wc -l)

        CORRECT=$((EXACT + CLOSE + HALF_CORRECT))
        if [ "$TOTAL" -gt 0 ]; then
            EXACT_PCT=$(echo "scale=1; $EXACT * 100 / $TOTAL" | bc)
            MISMATCH_PCT=$(echo "scale=1; $MISMATCH * 100 / $TOTAL" | bc)
            CORRECT_PCT=$(echo "scale=1; $CORRECT * 100 / $TOTAL" | bc)
        fi

        printf "%-25s | %4d (%5s%%) | %4d         | %4d       | %4d (%5s%%) | %4d (%5s%%)\n" \
            "$NAME" \
            "$EXACT" "$EXACT_PCT" "$HALF" "$THIRD" "$MISMATCH" "$MISMATCH_PCT" "$CORRECT" "$CORRECT_PCT"
    else
        echo "$NAME: analysis file not found ($ANALYSIS)"
    fi
done

echo ""
echo "* Correct = EXACT + CLOSE + HALF/DOUBLE(ratio<0.6)"
echo ""
echo "Done! Results in $RESULTS_DIR"
