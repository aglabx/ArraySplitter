#!/usr/bin/env bash
# Test harness for ArraySplitter (Rust engine)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CARGO_DIR="${SCRIPT_DIR}/src/rust/arraysplitter"

echo "========================================="
echo "Running ArraySplitter Rust Test Suite"
echo "========================================="
(cd "${CARGO_DIR}" && cargo test --release)

echo -e "\n========================================="
echo "Running Regression Harness"
echo "========================================="
if [[ -f "${SCRIPT_DIR}/test_data/zebra_finch_satdna.fasta" ]]; then
    bash "${SCRIPT_DIR}/tests/regression/run.sh"
else
    echo "Notice: test_data/zebra_finch_satdna.fasta not found locally; skipping panel regression."
fi

echo -e "\n========================================="
echo "All tests completed successfully!"
echo "========================================="