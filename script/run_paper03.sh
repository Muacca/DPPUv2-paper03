#!/usr/bin/env bash
# run_paper03.sh
# Run all paper03 and proofs scripts, saving stdout+stderr to the paper data directory.
#
# Usage (from the script/ directory):
#   bash run_paper03.sh                          # run all scripts  -> data/
#   bash run_paper03.sh proofs                   # run proofs/ only -> data/
#   bash run_paper03.sh paper03                  # run paper03/ only -> data/
#   bash run_paper03.sh --output-dir my_logs     # custom output dir
#   bash run_paper03.sh -o my_logs proofs        # combined
#
# Output: <output-dir>/<script_name>.log  (one file per script)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PAPER_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SCRIPTS_DIR="${SCRIPT_DIR}/scripts"

# ── argument parsing ──────────────────────────────────────────────────────────

TARGET="all"
DATA_DIR="${PAPER_ROOT}/data"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --output-dir|-o)
            DATA_DIR="$2"
            shift 2
            ;;
        all|proofs|paper03)
            TARGET="$1"
            shift
            ;;
        *)
            echo "Usage: bash run_paper03.sh [--output-dir DIR] [all|proofs|paper03]" >&2
            exit 1
            ;;
    esac
done

mkdir -p "${DATA_DIR}"

export PYTHONIOENCODING=utf-8
export PYTHONUTF8=1
export DPPU_LOG_DISABLE=1

# ── helper ───────────────────────────────────────────────────────────────────

run_script() {
    local rel_path="$1"          # e.g.  paper03/spin2_schur_degeneracy.py
    local script_file="${SCRIPTS_DIR}/${rel_path}"
    local base
    base="$(basename "${rel_path}" .py)"
    local log_file="${DATA_DIR}/${base}.log"

    echo ">>> ${rel_path}"
    local t0
    t0=$(date +%s)

    if python "${script_file}" > "${log_file}" 2>&1; then
        local elapsed=$(( $(date +%s) - t0 ))
        echo "    OK  (${elapsed}s)  -> ${DATA_DIR}/${base}.log"
    else
        local elapsed=$(( $(date +%s) - t0 ))
        echo "    FAILED  (${elapsed}s)  -> ${DATA_DIR}/${base}.log"
        FAILED_SCRIPTS+=("${rel_path}")
    fi
}

FAILED_SCRIPTS=()

# ── script lists ─────────────────────────────────────────────────────────────

PROOFS_SCRIPTS=(
    proofs/riemann_block_classification.py
    proofs/pontryagin_block_decomposition.py
    proofs/plucker_cancellation_proof.py
    proofs/pfaffian_decomposability.py
    proofs/vt_mode_cancellation.py
    proofs/analytic_proof_P_equals_zero.py
    proofs/analytical_proof.py
    proofs/symbolic_block_all_topologies.py
)

PAPER03_SCRIPTS=(
    paper03/pontryagin_scan_s3.py
    paper03/pontryagin_scan_nil3.py
    paper03/pontryagin_anisotropic_twist.py
    paper03/torsion_mode_analysis.py
    paper03/spin2_schur_degeneracy.py
    paper03/spin2_spectrum_splitting.py
    paper03/spin2_diagonal_hessians.py
    paper03/spin2_offdiag_hessian.py
    paper03/spin1_kinetic_metric.py
    paper03/spin1_null_mode_analysis.py
    paper03/spin1_auxiliary_exclusion.py
    paper03/spin1_appI_projection_artifact.py
    paper03/spin1_symbolic_firstprinciples.py
    paper03/spin1_cross_hessian.py
    paper03/homogeneous_consistency.py
    paper03/eft_invariant_generation.py
    paper03/spin0_torsion_cross_hessian.py
    paper03/spin0_reta_block_eigenvalues.py
    paper03/spin0_rLetaV_hessian.py
)

# ── dispatch ─────────────────────────────────────────────────────────────────

echo "========================================"
echo "  run_paper03.sh  (target: ${TARGET})"
echo "  output dir: ${DATA_DIR}"
echo "========================================"
echo ""

if [[ "${TARGET}" == "all" || "${TARGET}" == "proofs" ]]; then
    echo "--- proofs/ ---"
    for s in "${PROOFS_SCRIPTS[@]}"; do
        run_script "${s}"
    done
    echo ""
fi

if [[ "${TARGET}" == "all" || "${TARGET}" == "paper03" ]]; then
    echo "--- paper03/ ---"
    for s in "${PAPER03_SCRIPTS[@]}"; do
        run_script "${s}"
    done
    echo ""
fi

# ── summary ──────────────────────────────────────────────────────────────────

echo "========================================"
if [[ ${#FAILED_SCRIPTS[@]} -eq 0 ]]; then
    echo "  All scripts completed successfully."
else
    echo "  ${#FAILED_SCRIPTS[@]} script(s) failed:"
    for s in "${FAILED_SCRIPTS[@]}"; do
        echo "    - ${s}"
    done
fi
echo "========================================"
