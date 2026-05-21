#!/usr/bin/env bash
# Regression harness: builds the release binary, runs it on the locked
# input, and diffs the 5 output files (md5/size/lines) against the manifest.
#
# Usage:
#   bash tests/regression/run.sh                  # build + run + check
#   bash tests/regression/run.sh --no-build       # skip cargo build
#   bash tests/regression/run.sh --keep           # keep run outputs on disk
#   MANIFEST=tests/regression/<other>.manifest \
#       bash tests/regression/run.sh              # use a non-default manifest

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CARGO_DIR="${REPO_ROOT}/src/rust/arraysplitter"
BIN="${CARGO_DIR}/target/release/arraysplitter"
INPUT="${REPO_ROOT}/test_data/zebra_finch_satdna.fasta"
MANIFEST="${MANIFEST:-${SCRIPT_DIR}/zfinch_iter2.manifest.tsv}"

NO_BUILD=0
KEEP=0
for arg in "$@"; do
    case "${arg}" in
        --no-build) NO_BUILD=1 ;;
        --keep)     KEEP=1 ;;
        *) echo "Unknown arg: ${arg}" >&2; exit 2 ;;
    esac
done

[[ -f "${MANIFEST}" ]] || { echo "Manifest not found: ${MANIFEST}" >&2; exit 2; }
[[ -f "${INPUT}"    ]] || { echo "Input not found:    ${INPUT}"    >&2; exit 2; }

if [[ "${NO_BUILD}" -eq 0 ]]; then
    echo "[build] cargo build --release"
    (cd "${CARGO_DIR}" && cargo build --release --quiet)
fi
[[ -x "${BIN}" ]] || { echo "Binary missing: ${BIN}" >&2; exit 2; }

OUT_DIR="${REPO_ROOT}/results/regression_runs/zfinch_$(date +%Y%m%d_%H%M%S)_$$"
mkdir -p "${OUT_DIR}"
PREFIX="${OUT_DIR}/zfinch"

cleanup() {
    if [[ "${KEEP}" -eq 0 ]]; then
        /bin/rm -rf "${OUT_DIR}"
    else
        echo "[keep] outputs preserved at ${OUT_DIR}"
    fi
}
trap cleanup EXIT

echo "[run] ${BIN} -i <input> -o ${PREFIX} -t 4 --method autocorr"
"${BIN}" -i "${INPUT}" -o "${PREFIX}" -t 4 --method autocorr >/dev/null 2>&1

# md5 implementation differs between macOS (md5 -q) and Linux (md5sum).
if command -v md5 >/dev/null 2>&1; then
    md5_of() { md5 -q "$1"; }
elif command -v md5sum >/dev/null 2>&1; then
    md5_of() { md5sum "$1" | awk '{print $1}'; }
else
    echo "Neither md5 nor md5sum found" >&2; exit 2
fi

# stat byte-size flag differs between BSD (-f%z) and GNU (-c%s).
if stat -f%z /dev/null >/dev/null 2>&1; then
    size_of() { stat -f%z "$1"; }
else
    size_of() { stat -c%s "$1"; }
fi

fail=0
checked=0
printf '\n%-22s %-7s %-7s %-7s\n' "FILE" "MD5" "SIZE" "LINES"
printf '%-22s %-7s %-7s %-7s\n' "----" "---" "----" "-----"

# Iterate manifest rows (skip comments and blank lines).
while IFS=$'\t' read -r ext exp_md5 exp_size exp_lines; do
    [[ -z "${ext}" || "${ext}" =~ ^# ]] && continue
    checked=$((checked + 1))
    f="${PREFIX}${ext}"
    if [[ ! -f "${f}" ]]; then
        printf '%-22s %s\n' "${ext}" "MISSING (file not produced)"
        fail=$((fail + 1))
        continue
    fi
    got_md5=$(md5_of "${f}")
    got_size=$(size_of "${f}")
    got_lines=$(wc -l < "${f}" | tr -d ' ')

    md5_ok="OK";   [[ "${got_md5}"   == "${exp_md5}"   ]] || md5_ok="FAIL"
    size_ok="OK";  [[ "${got_size}"  == "${exp_size}"  ]] || size_ok="FAIL"
    lines_ok="OK"; [[ "${got_lines}" == "${exp_lines}" ]] || lines_ok="FAIL"

    printf '%-22s %-7s %-7s %-7s\n' "${ext}" "${md5_ok}" "${size_ok}" "${lines_ok}"

    if [[ "${md5_ok}" != "OK" ]]; then
        printf '   md5  expected=%s got=%s\n' "${exp_md5}" "${got_md5}"
        fail=$((fail + 1))
    fi
    if [[ "${size_ok}" != "OK" ]]; then
        printf '   size expected=%s got=%s\n' "${exp_size}" "${got_size}"
        fail=$((fail + 1))
    fi
    if [[ "${lines_ok}" != "OK" ]]; then
        printf '   lines expected=%s got=%s\n' "${exp_lines}" "${got_lines}"
        fail=$((fail + 1))
    fi
done < "${MANIFEST}"

echo
if [[ "${fail}" -eq 0 && "${checked}" -gt 0 ]]; then
    echo "PASS — ${checked}/${checked} files match manifest"
    exit 0
else
    echo "FAIL — ${fail} mismatch(es) across ${checked} file(s)"
    exit 1
fi
