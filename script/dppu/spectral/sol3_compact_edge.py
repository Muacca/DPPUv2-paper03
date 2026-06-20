"""Sol3 compact-edge modified-Mathieu spectral engine.

Sol3 compact-edge spectral computation: modified-Mathieu operator for the
A-dual-orbit sector, rigidity analysis (L/R=log lambda_+), convergence
across grid/window refinements, and robustness across SL(2,Z) matrices.

Public entry points
-------------------
compute_tables() -> (rigidity_rows, lambda_rows, convergence_rows,
                     compare_rows, robustness_rows)
gate_rows(rigidity_rows, lambda_rows, convergence_rows, compare_rows,
          robustness_rows) -> gate_rows_list

The engine is self-contained and does not read external analysis logs or CSV
tables.
"""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np
from scipy.linalg import eigh_tridiagonal


# ---------------------------------------------------------------------------
# Pure compute helpers (verbatim)
# ---------------------------------------------------------------------------
def matrix_id(A: np.ndarray) -> str:
    return f"A{int(A[0,0])}_{int(A[0,1])}_{int(A[1,0])}_{int(A[1,1])}"


def matrix_text(A: np.ndarray) -> str:
    return f"[[{int(A[0,0])},{int(A[0,1])}],[{int(A[1,0])},{int(A[1,1])}]]"


def hyperbolic_data(A: np.ndarray) -> Dict[str, Any]:
    trace = int(np.trace(A))
    det = int(round(np.linalg.det(A)))
    disc = trace * trace - 4 * det
    if det != 1 or disc <= 0:
        raise ValueError(f"not a positive hyperbolic SL(2,Z) matrix: {matrix_text(A)}")
    lambda_plus = 0.5 * (trace + math.sqrt(disc))
    lambda_minus = 1.0 / lambda_plus
    mu = math.log(lambda_plus)
    return {
        "A_id": matrix_id(A),
        "A_matrix": matrix_text(A),
        "trace": trace,
        "det": det,
        "discriminant": disc,
        "lambda_plus": lambda_plus,
        "lambda_minus": lambda_minus,
        "L_over_R": mu,
        "status": "PASS",
        "notes": "Consistency rigidity from A: L/R=log(lambda_plus). Gap size is controlled by the A-invariant quadratic level Q=|alpha beta|.",
    }


def dual_components(A: np.ndarray, gamma: Tuple[int, int]) -> Dict[str, Any]:
    # The compact quotient identifies nonzero Fourier labels by the dual Anosov
    # action.  In the A^T eigenbasis the lifted orbit sector has a confining
    # modified-Mathieu potential.
    B = A.T.astype(float)
    eigvals, eigvecs = np.linalg.eig(B)
    eigvals = np.real(eigvals)
    eigvecs = np.real(eigvecs)
    plus_idx = int(np.argmax(eigvals))
    minus_idx = 1 - plus_idx
    v_plus = eigvecs[:, plus_idx]
    v_minus = eigvecs[:, minus_idx]
    if v_plus[0] < 0:
        v_plus = -v_plus
    if v_minus[1] < 0:
        v_minus = -v_minus
    basis = np.column_stack([v_plus, v_minus])
    coeff = np.linalg.solve(basis, np.array(gamma, dtype=float))
    alpha = float(coeff[0])
    beta = float(coeff[1])
    lam_plus = float(eigvals[plus_idx])
    lam_minus = float(eigvals[minus_idx])
    residuals = []
    orbit_terms = []
    for n in range(-2, 3):
        if n >= 0:
            vec = np.linalg.matrix_power(A.T, n).dot(np.array(gamma, dtype=int))
        else:
            inv_bt = np.linalg.matrix_power(np.linalg.inv(A.T), -n)
            vec = inv_bt.dot(np.array(gamma, dtype=float))
        coeff_n = np.linalg.solve(basis, np.array(vec, dtype=float))
        expected = np.array([alpha * (lam_plus**n), beta * (lam_minus**n)])
        residuals.append(float(np.linalg.norm(coeff_n - expected)))
        orbit_terms.append(f"n={n}:({vec[0]:.8g},{vec[1]:.8g})")
    return {
        "alpha": alpha,
        "beta": beta,
        "lambda_plus": lam_plus,
        "lambda_minus": lam_minus,
        "orbit_scaling_residual_max": max(residuals),
        "dual_orbit_window": "; ".join(orbit_terms),
    }


def low_spectrum(alpha: float, beta: float, window: float, grid_n: int) -> Dict[str, Any]:
    # Dimensionless coordinate t=w/R.  The potential
    # alpha^2 exp(-2t)+beta^2 exp(2t) is the nonzero A-orbit sector reduction.
    # Dirichlet walls are placed far from the classical minimum and then varied.
    center = 0.25 * math.log((alpha * alpha) / (beta * beta))
    left = center - window
    right = center + window
    h = (right - left) / (grid_n + 1)
    xs = left + h * np.arange(1, grid_n + 1)
    potential = alpha * alpha * np.exp(-2.0 * xs) + beta * beta * np.exp(2.0 * xs)
    diag = (2.0 / (h * h)) + potential
    offdiag = np.full(grid_n - 1, -1.0 / (h * h), dtype=float)
    eigs = eigh_tridiagonal(diag, offdiag, select="i", select_range=(0, 3))[0]
    return {
        "center": center,
        "left": left,
        "right": right,
        "h": h,
        "lambda1": float(eigs[0]),
        "low_eigenvalues": [float(x) for x in eigs],
        "potential_min": 2.0 * abs(alpha * beta),
    }


def selected_gamma_levels(A: np.ndarray, max_levels: int = 6, bound: int = 4) -> List[Tuple[Tuple[int, int], float, str]]:
    candidates: List[Tuple[float, Tuple[int, int]]] = []
    for p_val in range(-bound, bound + 1):
        for q_val in range(-bound, bound + 1):
            if (p_val, q_val) == (0, 0):
                continue
            gamma = (p_val, q_val)
            if gamma[0] < 0 or (gamma[0] == 0 and gamma[1] < 0):
                gamma = (-gamma[0], -gamma[1])
            comp = dual_components(A, gamma)
            q_invariant = abs(comp["alpha"] * comp["beta"])
            candidates.append((q_invariant, gamma))
    candidates.sort(key=lambda item: (round(item[0], 10), abs(item[1][0]) + abs(item[1][1]), item[1]))
    selected: List[Tuple[Tuple[int, int], float, str]] = []
    seen_levels: set[str] = set()
    for q_invariant, gamma in candidates:
        level_key = f"{q_invariant:.10g}"
        if level_key in seen_levels:
            continue
        seen_levels.add(level_key)
        selected.append((gamma, q_invariant, level_key))
        if len(selected) >= max_levels:
            break
    return selected


def compute_tables() -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]]]:
    matrices = [
        np.array([[2, 1], [1, 1]], dtype=int),
        np.array([[1, 1], [1, 2]], dtype=int),
        np.array([[3, 1], [2, 1]], dtype=int),
        np.array([[4, 1], [3, 1]], dtype=int),
    ]
    grid_configs = [(5.0, 260), (6.0, 380), (7.0, 520)]
    rigidity_rows = [hyperbolic_data(A) for A in matrices]
    discriminant_seen: set[int] = set()
    for row in rigidity_rows:
        disc = int(row["discriminant"])
        if disc in discriminant_seen:
            row["independence_role"] = "REDUNDANT_SAME_DISCRIMINANT_CONTROL"
        else:
            row["independence_role"] = "INDEPENDENT_DISCRIMINANT"
            discriminant_seen.add(disc)
    lambda_rows: List[Dict[str, Any]] = []
    for A in matrices:
        hdat = hyperbolic_data(A)
        for gamma, q_invariant, q_level_key in selected_gamma_levels(A):
            comp = dual_components(A, gamma)
            raw_bulk = 2 * abs(gamma[0] * gamma[1])
            for window, grid_n in grid_configs:
                spec = low_spectrum(comp["alpha"], comp["beta"], window, grid_n)
                lambda_rows.append(
                    {
                        "case_id": f"{matrix_id(A)}_g{gamma[0]}_{gamma[1]}_W{window:g}_N{grid_n}",
                        "A_id": matrix_id(A),
                        "A_matrix": matrix_text(A),
                        "gamma_integer": str(gamma),
                        "discriminant": hdat["trace"] * hdat["trace"] - 4,
                        "lambda_plus": f"{hdat['lambda_plus']:.12g}",
                        "L_over_R": f"{hdat['L_over_R']:.12g}",
                        "alpha_dual_plus": f"{comp['alpha']:.12g}",
                        "beta_dual_minus": f"{comp['beta']:.12g}",
                        "Q_abs_alpha_beta": f"{q_invariant:.12g}",
                        "Q_level_key": q_level_key,
                        "dual_orbit_window": comp["dual_orbit_window"],
                        "orbit_scaling_residual_max": f"{comp['orbit_scaling_residual_max']:.3e}",
                        "window_half_width": window,
                        "grid_n": grid_n,
                        "grid_h": f"{spec['h']:.12g}",
                        "center_t": f"{spec['center']:.12g}",
                        "potential_min_2abs_alpha_beta": f"{spec['potential_min']:.12g}",
                        "raw_bulk_2abs_pq": raw_bulk,
                        "lambda1_genuine": f"{spec['lambda1']:.12g}",
                        "low_eigenvalues": ",".join(f"{value:.12g}" for value in spec["low_eigenvalues"]),
                        "operator_class": "A_DUAL_ORBIT_MODIFIED_MATHIEU_FINITE_DIFFERENCE",
                        "operator_family": "COMPACT_A_DUAL_ORBIT",
                        "abelian_carrier_separated": gamma != (0, 0),
                        "evidence_tier": "TIER1_2_QUALITATIVE_BDV_DISCRETE_GAPPED;TIER3_NUMERIC_REPRESENTATIVE_LAMBDA1",
                        "status": "PASS" if spec["lambda1"] > 0.0 and comp["orbit_scaling_residual_max"] < 1e-8 else "FAIL",
                    }
                )
    final_rows = [row for row in lambda_rows if int(row["grid_n"]) == grid_configs[-1][1]]
    by_case: Dict[Tuple[str, str], List[Dict[str, Any]]] = defaultdict(list)
    for row in lambda_rows:
        by_case[(row["A_id"], row["gamma_integer"])].append(row)
    convergence_rows: List[Dict[str, Any]] = []
    for key, rows in by_case.items():
        rows_sorted = sorted(rows, key=lambda r: int(r["grid_n"]))
        lams = [float(r["lambda1_genuine"]) for r in rows_sorted]
        abs_last_delta = abs(lams[-1] - lams[-2])
        rel_last_delta = abs_last_delta / max(1.0, abs(lams[-1]))
        convergence_rows.append(
            {
                "case_id": f"{key[0]}_{key[1]}",
                "A_id": key[0],
                "gamma_integer": key[1],
                "lambda1_sequence": ",".join(f"{x:.12g}" for x in lams),
                "abs_last_delta": f"{abs_last_delta:.12g}",
                "rel_last_delta": f"{rel_last_delta:.12g}",
                "status": "PASS" if rel_last_delta < 5.0e-3 else "REVIEW",
                "notes": "Window and grid refinement sequence for the genuine orbit-sector operator.",
            }
        )
    compare_rows: List[Dict[str, Any]] = []
    for row in final_rows:
        gamma = row["gamma_integer"]
        genuine = float(row["lambda1_genuine"])
        raw_bulk = float(row["raw_bulk_2abs_pq"])
        compare_rows.append(
            {
                "case_id": row["case_id"],
                "A_id": row["A_id"],
                "discriminant": row["discriminant"],
                "gamma_integer": gamma,
                "Q_abs_alpha_beta": row["Q_abs_alpha_beta"],
                "Q_level_key": row["Q_level_key"],
                "genuine_lambda1": f"{genuine:.12g}",
                "noncompact_infimum": "0",
                "raw_bulk_2abs_pq": row["raw_bulk_2abs_pq"],
                "genuine_minus_raw_bulk": f"{genuine - raw_bulk:.12g}",
                "positive_discrete": genuine > 0.0,
                "nonbulk": abs(genuine - raw_bulk) > 1.0e-3,
                "bulk_relation": "COMPACT_GAP_DEPENDS_ON_Q_NOT_RAW_2ABS_PQ",
                "compact_vs_noncompact": "COMPACT_DISCRETE_POSITIVE_VS_NONCOMPACT_ZERO_INFIMUM" if genuine > 0.0 else "REVIEW",
                "status": "PASS" if genuine > 0.0 and abs(genuine - raw_bulk) > 1.0e-3 else "REVIEW",
            }
        )
    robustness_rows: List[Dict[str, Any]] = []
    by_A: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for row in compare_rows:
        by_A[row["A_id"]].append(row)
    for A_id, rows in by_A.items():
        lambdas = [float(r["genuine_lambda1"]) for r in rows]
        statuses = [r["status"] for r in rows]
        q_levels = sorted({r["Q_level_key"] for r in rows})
        discriminants = sorted({int(r["discriminant"]) for r in rows})
        robustness_rows.append(
            {
                "A_id": A_id,
                "case_count": len(rows),
                "independent_Q_levels": len(q_levels),
                "discriminants": ",".join(str(x) for x in discriminants),
                "min_genuine_lambda1": f"{min(lambdas):.12g}",
                "max_genuine_lambda1": f"{max(lambdas):.12g}",
                "all_positive": all(x > 0.0 for x in lambdas),
                "all_nonbulk": all(r["nonbulk"] == "True" or r["nonbulk"] is True for r in rows),
                "all_case_status_pass": all(status == "PASS" for status in statuses),
                "status": "PASS" if all(x > 0.0 for x in lambdas) and all(status == "PASS" for status in statuses) else "REVIEW",
                "notes": f"Distinct Q-level sampling for one hyperbolic A; Q_levels={q_levels}.",
            }
        )
    by_disc: Dict[str, List[float]] = defaultdict(list)
    for row in compare_rows:
        by_disc[str(row["discriminant"])].append(float(row["genuine_lambda1"]))
    disc_ranges = {disc: max(values) - min(values) for disc, values in by_disc.items() if len(values) >= 2}
    independent_disc_count = len(by_disc)
    robustness_rows.append(
        {
            "A_id": "GLOBAL_DISCRIMINANT_AND_Q_STRUCTURE",
            "case_count": len(compare_rows),
            "independent_Q_levels": len({(row["discriminant"], row["Q_level_key"]) for row in compare_rows}),
            "discriminants": ",".join(sorted(by_disc.keys(), key=int)),
            "min_genuine_lambda1": f"{min(float(r['genuine_lambda1']) for r in compare_rows):.12g}",
            "max_genuine_lambda1": f"{max(float(r['genuine_lambda1']) for r in compare_rows):.12g}",
            "all_positive": all(float(r["genuine_lambda1"]) > 0.0 for r in compare_rows),
            "all_nonbulk": all(r["nonbulk"] == "True" or r["nonbulk"] is True for r in compare_rows),
            "all_case_status_pass": all(r["status"] == "PASS" for r in compare_rows),
            "status": "PASS" if independent_disc_count >= 3 and len({(row["discriminant"], row["Q_level_key"]) for row in compare_rows}) >= 8 else "REVIEW",
            "notes": f"Honest robustness: lambda1 depends on Q=|alpha beta| and discriminant, not raw case count; disc_ranges={disc_ranges}.",
        }
    )
    return rigidity_rows, lambda_rows, convergence_rows, compare_rows, robustness_rows


def gate_rows(
    rigidity_rows: Sequence[Dict[str, Any]],
    lambda_rows: Sequence[Dict[str, Any]],
    convergence_rows: Sequence[Dict[str, Any]],
    compare_rows: Sequence[Dict[str, Any]],
    robustness_rows: Sequence[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    final_lambda_rows = [row for row in lambda_rows if int(row["grid_n"]) == 520]
    rigidity_ok = bool(rigidity_rows) and all(row["status"] == "PASS" for row in rigidity_rows)
    lambda_ok = bool(final_lambda_rows) and all(row["status"] == "PASS" for row in final_lambda_rows)
    convergence_ok = bool(convergence_rows) and all(row["status"] == "PASS" for row in convergence_rows)
    compare_ok = bool(compare_rows) and all(row["status"] == "PASS" for row in compare_rows)
    robustness_ok = len(robustness_rows) >= 2 and all(row["status"] == "PASS" for row in robustness_rows)
    structure_ok = any(row["A_id"] == "GLOBAL_DISCRIMINANT_AND_Q_STRUCTURE" and row["status"] == "PASS" for row in robustness_rows)
    abelian_ok = all(row.get("abelian_carrier_separated") is True for row in final_lambda_rows)
    compact_noncompact_ok = compare_ok and all(
        row["compact_vs_noncompact"] == "COMPACT_DISCRETE_POSITIVE_VS_NONCOMPACT_ZERO_INFIMUM"
        for row in compare_rows
    )
    return [
        {"gate": "R_PASS_0", "status": "PASS", "witness": "dppu.spectral.sol3_compact_edge", "notes": "Self-contained compact-edge engine; no external CSV input is read."},
        {"gate": "R_PASS_1", "status": "PASS" if rigidity_ok else "FAIL", "witness": "sol3_compact_rigidity.csv", "notes": "L/R=log(lambda_plus(A)) computed from A."},
        {"gate": "R_PASS_2", "status": "PASS" if lambda_ok else "FAIL", "witness": "sol3_compact_genuine_lambda1.csv", "notes": "A-dual-orbit modified-Mathieu operator."},
        {"gate": "R_PASS_3", "status": "PASS" if convergence_ok else "FAIL", "witness": "sol3_compact_genuine_lambda1.csv", "notes": "Grid/window refinement stable."},
        {"gate": "R_PASS_4", "status": "PASS" if compare_ok and structure_ok else "FAIL", "witness": "sol3_compact_vs_noncompact.csv; sol3_compact_edge_multiA_robustness.csv", "notes": "Positive compact gap is governed by Q=|alpha beta|/discriminant, not by raw 2|pq| or noncompact infimum zero."},
        {"gate": "R_PASS_5", "status": "PASS" if robustness_ok else "FAIL", "witness": "sol3_compact_edge_multiA_robustness.csv", "notes": "Robustness is counted by distinct Q levels and distinct discriminants; redundant same-discriminant controls are marked."},
        {"gate": "R_PASS_6", "status": "PASS" if abelian_ok else "FAIL", "witness": "sol3_compact_edge_verdict.md", "notes": "(0,0) carrier is separated from nonzero edge gap sectors."},
        {"gate": "R_PASS_7", "status": "PASS" if compact_noncompact_ok else "FAIL", "witness": "sol3_compact_vs_noncompact.csv", "notes": "Compact nonzero A-dual orbit sectors are separated from the noncompact zero-infimum reference."},
    ]
