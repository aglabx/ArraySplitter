#!/usr/bin/env bash
# run_all.sh: Run the complete benchmark and validation suite for ArraySplitter.
#
# Usage:
#   bash scripts/benchmark/run_all.sh [OPTIONS]
#
# Options:
#   --threads <N>       Worker threads for arraysplitter (default: 4)
#   --out-dir <DIR>     Output directory for benchmark results (default: results/benchmark_run)
#   --bin <PATH>        Path to arraysplitter binary (default: auto-detect)
#   --skip-run          Skip decomposition run and reuse existing output prefix
#   --keep              Keep output files on exit (default: keep)
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

THREADS=4
OUT_DIR="${REPO_ROOT}/results/benchmark_run"
BIN=""
SKIP_RUN=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads)
            THREADS="$2"; shift 2 ;;
        --out-dir)
            OUT_DIR="$2"; shift 2 ;;
        --bin)
            BIN="$2"; shift 2 ;;
        --skip-run)
            SKIP_RUN=1; shift ;;
        --keep)
            shift ;; # default behavior
        *)
            echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

# Locate binary
if [[ -z "${BIN}" ]]; then
    if [[ -x "${REPO_ROOT}/src/rust/arraysplitter/target/release/arraysplitter" ]]; then
        BIN="${REPO_ROOT}/src/rust/arraysplitter/target/release/arraysplitter"
    elif command -v arraysplitter >/dev/null 2>&1; then
        BIN="$(command -v arraysplitter)"
    else
        echo "[error] arraysplitter binary not found. Please run: (cd src/rust/arraysplitter && cargo build --release)" >&2
        exit 1
    fi
fi

echo "=== ArraySplitter Benchmark Suite ==="
echo "Repo root : ${REPO_ROOT}"
echo "Binary    : ${BIN}"
echo "Threads   : ${THREADS}"
echo "Out dir   : ${OUT_DIR}"
echo

mkdir -p "${OUT_DIR}"
PREFIX="${OUT_DIR}/zf_run"
ZF_FA="${REPO_ROOT}/test_data/zebra_finch_satdna.fasta"
ZF_GZ="${REPO_ROOT}/test_data/zebra_finch_satdna.fasta.gz"

if [[ ! -f "${ZF_FA}" ]]; then
    if [[ -f "${ZF_GZ}" ]]; then
        echo "[setup] Decompressing ${ZF_GZ}..."
        gzip -dc "${ZF_GZ}" > "${ZF_FA}"
    else
        echo "[error] Test panel not found: ${ZF_GZ}" >&2
        exit 1
    fi
fi

if [[ "${SKIP_RUN}" -eq 0 ]]; then
    echo "[1/6] Running ArraySplitter on zebra finch satDNA panel..."
    "${BIN}" -i "${ZF_FA}" -o "${PREFIX}" -t "${THREADS}" --method autocorr
    echo "      Decomposition complete: $(wc -l < "${PREFIX}.summary.tsv") arrays in summary.tsv"
else
    echo "[1/6] Skipping ArraySplitter run (reusing ${PREFIX})"
fi

echo
echo "[2/6] Validating parent_idx hierarchy joins (check_join.py)..."
python3 "${SCRIPT_DIR}/check_join.py" "${PREFIX}"

echo
echo "[3/6] Null control for cut agreement metric (zf_metric_null.py)..."
python3 "${SCRIPT_DIR}/zf_metric_null.py" "${PREFIX}" "${REPO_ROOT}"

echo
echo "[4/6] Autocorrelation support analysis for HOR calls (zf_hor_support.py)..."
python3 "${SCRIPT_DIR}/zf_hor_support.py" "${ZF_FA}" "${PREFIX}"

echo
echo "[5/6] Sub-HOR structure, duplicates, and edge monomers (zf_subhor_and_edges.py)..."
python3 "${SCRIPT_DIR}/zf_subhor_and_edges.py" "${ZF_FA}" "${PREFIX}"

echo
echo "[6/6] Subsampling aliasing & sampling error audit (alias_check.py + sampling_error.py)..."
python3 "${SCRIPT_DIR}/alias_check.py" "${OUT_DIR}" "${REPO_ROOT}"
python3 "${SCRIPT_DIR}/sampling_error.py" "${OUT_DIR}/alias/alias.fa"

echo
echo "=== All benchmark suite checks completed successfully ==="
