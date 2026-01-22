#!/bin/bash
# Test different LAMBDA values for dual ED optimization
# Usage: ./test_lambda_values.sh

set -e

WORKDIR="/mnt/data/claude/arraysplitter_rust"
INPUT="/mnt/data/claude/phase1_alpha_extraction/all_alpha_satellite.fasta"
RUST_SRC="$WORKDIR/src/anchor_graph.rs"
RESULTS_DIR="$WORKDIR/lambda_tests"

# Lambda values to test
LAMBDAS="0.0 0.5 1.0 1.5 2.0 3.0"

mkdir -p "$RESULTS_DIR"

echo "=== Testing different LAMBDA values ==="
echo "Input: $INPUT"
echo "Results: $RESULTS_DIR"
echo ""

# Save original file
cp "$RUST_SRC" "$RUST_SRC.backup"

for LAMBDA in $LAMBDAS; do
    echo "=========================================="
    echo "Testing LAMBDA = $LAMBDA"
    echo "=========================================="

    # Update LAMBDA in Rust code
    sed -i "s/const LAMBDA: f64 = [0-9.]*;/const LAMBDA: f64 = $LAMBDA;/" "$RUST_SRC"

    # Verify change
    grep "const LAMBDA" "$RUST_SRC"

    # Rebuild
    echo "Building..."
    cd "$WORKDIR"
    ~/.cargo/bin/cargo build --release 2>&1 | tail -2

    # Run analysis
    OUTPUT_PREFIX="$RESULTS_DIR/lambda_${LAMBDA//./_}"
    echo "Running analysis -> $OUTPUT_PREFIX"

    START_TIME=$(date +%s)
    ./target/release/arraysplitter_rs -i "$INPUT" -o "$OUTPUT_PREFIX" 2>&1 | tail -5
    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))

    echo "Elapsed: ${ELAPSED}s"

    # Analyze results
    echo "Analyzing results..."
    python3 "$WORKDIR/scripts/analyze_decomposition.py" \
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
echo "SUMMARY: LAMBDA comparison"
echo "=========================================="
echo "LAMBDA | EXACT | CLOSE | HALF/DBL | 3RD/3PL | MISMATCH"
echo "-------|-------|-------|----------|---------|----------"

for LAMBDA in $LAMBDAS; do
    OUTPUT_PREFIX="$RESULTS_DIR/lambda_${LAMBDA//./_}"
    if [ -f "${OUTPUT_PREFIX}.analysis.tsv" ]; then
        # Count matches from analysis file
        TOTAL=$(tail -n +2 "${OUTPUT_PREFIX}.analysis.tsv" | wc -l)
        EXACT=$(tail -n +2 "${OUTPUT_PREFIX}.analysis.tsv" | awk -F'\t' '$12=="EXACT"' | wc -l)
        CLOSE=$(tail -n +2 "${OUTPUT_PREFIX}.analysis.tsv" | awk -F'\t' '$12=="CLOSE"' | wc -l)
        HALF=$(tail -n +2 "${OUTPUT_PREFIX}.analysis.tsv" | awk -F'\t' '$12=="HALF/DOUBLE"' | wc -l)
        THIRD=$(tail -n +2 "${OUTPUT_PREFIX}.analysis.tsv" | awk -F'\t' '$12=="THIRD/TRIPLE"' | wc -l)
        MISMATCH=$(tail -n +2 "${OUTPUT_PREFIX}.analysis.tsv" | awk -F'\t' '$12=="MISMATCH"' | wc -l)

        EXACT_PCT=$(echo "scale=1; $EXACT * 100 / $TOTAL" | bc)
        MISMATCH_PCT=$(echo "scale=1; $MISMATCH * 100 / $TOTAL" | bc)

        printf "%6s | %5s | %5s | %8s | %7s | %s\n" \
            "$LAMBDA" "$EXACT ($EXACT_PCT%)" "$CLOSE" "$HALF" "$THIRD" "$MISMATCH ($MISMATCH_PCT%)"
    fi
done

echo ""
echo "Done! Results in $RESULTS_DIR"
