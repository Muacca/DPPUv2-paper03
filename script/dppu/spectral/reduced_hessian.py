"""DBR reduced-Hessian audit on the constraint-surface reduced phase space.

Pure compute module — no I/O, no logging, no side effects on import.

Provides:
  — scalar_axis_hessian_audits  (repair_mod parameter dropped; imports
    build_scalar_axis_projection and dirac_standard from dppu.dirac.scalar_axis)
  — reduced_hessian_pair_audit
  — (helpers: poisson, matrix_rows)

The ``repair_mod`` parameter has been removed.  ``build_scalar_axis_projection``
and ``dirac_standard`` are imported directly from ``dppu.dirac.scalar_axis``.
"""

from __future__ import annotations

from typing import Any, Dict, List, Sequence

import sympy as sp

from dppu.dirac.scalar_axis import build_scalar_axis_projection, dirac_standard


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def poisson(f: sp.Expr, g: sp.Expr, q_vars: Sequence[sp.Symbol], p_vars: Sequence[sp.Symbol]) -> sp.Expr:
    return sp.expand(sum(sp.diff(f, q) * sp.diff(g, p) - sp.diff(f, p) * sp.diff(g, q) for q, p in zip(q_vars, p_vars)))


def matrix_rows(matrix: sp.Matrix) -> List[List[str]]:
    return [[str(sp.simplify(matrix[i, j])) for j in range(matrix.cols)] for i in range(matrix.rows)]


# ---------------------------------------------------------------------------
# Public entry points
# ---------------------------------------------------------------------------

def reduced_hessian_pair_audit(result: Dict[str, Any]) -> Dict[str, Any]:
    q_cfg: List[sp.Symbol] = result["q_cfg"]
    p_cfg: List[sp.Symbol] = result["p_cfg"]
    vars_all = q_cfg + p_cfg
    n_config = len(q_cfg)
    if result["N_phys"] == 0:
        return {
            "basis_type": "constraint_surface_reduced_phase_space",
            "phase_basis": [str(var) for var in vars_all],
            "constraint_surface_dim": 0,
            "physical_symplectic_rank": 0,
            "physical_pair_columns": [],
            "H_pair": [],
            "H_times_pair_zero": True,
            "ker_H_pair_dim": 0,
            "pair_subspace_dim": 0,
            "ImPhi_equals_kerH": True,
            "status": "PASS",
            "notes": "No scalar zero-mode pair in this representative.",
        }

    constraint_matrix = sp.Matrix([[sp.diff(c, var) for var in vars_all] for c in result["constraints"]])
    null_cols = constraint_matrix.nullspace()
    if not null_cols:
        return {
            "basis_type": "constraint_surface_reduced_phase_space",
            "phase_basis": [str(var) for var in vars_all],
            "constraint_surface_dim": 0,
            "physical_symplectic_rank": 0,
            "physical_pair_columns": [],
            "H_pair": [],
            "H_times_pair_zero": False,
            "ker_H_pair_dim": 0,
            "pair_subspace_dim": 0,
            "ImPhi_equals_kerH": False,
            "status": "FAIL",
            "notes": "Constraint-surface nullspace was empty.",
        }

    basis = sp.Matrix.hstack(*null_cols)
    symplectic = sp.zeros(2 * n_config, 2 * n_config)
    for idx in range(n_config):
        symplectic[idx, n_config + idx] = 1
        symplectic[n_config + idx, idx] = -1
    omega_reduced = basis.T * symplectic * basis
    physical_rank = int(omega_reduced.rank())
    hessian = sp.Matrix([[sp.diff(result["H_c"], a, b) for b in vars_all] for a in vars_all]) / 2

    selected_pair = None
    selected_h_pair = None
    for i in range(basis.shape[1]):
        for j in range(i + 1, basis.shape[1]):
            omega_pair = sp.Matrix([[omega_reduced[i, i], omega_reduced[i, j]], [omega_reduced[j, i], omega_reduced[j, j]]])
            if sp.simplify(omega_pair.det()) == 0:
                continue
            sub_basis = basis[:, [i, j]]
            h_pair = sp.simplify(sub_basis.T * hessian * sub_basis)
            if h_pair == sp.zeros(2, 2):
                selected_pair = sub_basis
                selected_h_pair = h_pair
                break
        if selected_pair is not None:
            break

    if selected_pair is None or selected_h_pair is None:
        return {
            "basis_type": "constraint_surface_reduced_phase_space",
            "phase_basis": [str(var) for var in vars_all],
            "constraint_surface_dim": basis.shape[1],
            "physical_symplectic_rank": physical_rank,
            "physical_pair_columns": [],
            "H_pair": [],
            "H_times_pair_zero": False,
            "ker_H_pair_dim": 0,
            "pair_subspace_dim": 0,
            "ImPhi_equals_kerH": False,
            "status": "FAIL",
            "notes": "No nondegenerate zero-Hessian symplectic pair was found.",
        }

    h_times_pair = sp.simplify(hessian * selected_pair)
    h_times_pair_zero = h_times_pair == sp.zeros(h_times_pair.rows, h_times_pair.cols)
    ker_dim = len(selected_h_pair.nullspace())
    pair_dim = selected_pair.rank()
    im_equals_ker = h_times_pair_zero and ker_dim == pair_dim == 2
    return {
        "basis_type": "constraint_surface_reduced_phase_space",
        "phase_basis": [str(var) for var in vars_all],
        "constraint_surface_dim": basis.shape[1],
        "physical_symplectic_rank": physical_rank,
        "physical_pair_columns": matrix_rows(selected_pair),
        "H_pair": matrix_rows(selected_h_pair),
        "H_times_pair_zero": bool(h_times_pair_zero),
        "ker_H_pair_dim": int(ker_dim),
        "pair_subspace_dim": int(pair_dim),
        "ImPhi_equals_kerH": bool(im_equals_ker),
        "status": "PASS" if im_equals_ker and physical_rank == 2 else "FAIL",
        "notes": "One nondegenerate constraint-surface symplectic pair has zero reduced Hamiltonian Hessian.",
    }


def scalar_axis_hessian_audits(topology: str, local: Dict[str, Any], active_axes: Sequence[str]) -> Dict[str, Any]:
    """Run scalar-axis DBR Hessian audits for all three axes.

    The ``repair_mod`` parameter has been dropped; ``build_scalar_axis_projection``
    and ``dirac_standard`` are imported directly from
    ``dppu.dirac.scalar_axis``.
    """
    rows: Dict[str, Any] = {}
    for axis in ("x", "y", "z"):
        proj = build_scalar_axis_projection(local, axis)
        points = [
            ("endpoint_symbolic", {local["N"]: sp.S.One, local["q"]: sp.S.One, proj["k"]: sp.S.One}),
            ("generic_numeric_1", {local["N"]: sp.Integer(2), local["q"]: sp.Integer(3), proj["k"]: sp.Integer(2)}),
            ("generic_numeric_2", {local["N"]: sp.Integer(5), local["q"]: sp.Integer(7), proj["k"]: sp.Integer(3)}),
        ]
        point_rows = []
        for point_id, subs in points:
            result = dirac_standard(proj, local, subs)
            audit = reduced_hessian_pair_audit(result)
            point_rows.append(
                {
                    "point_id": point_id,
                    "N_phys": str(result["N_phys"]),
                    "rank_C": int(result["rank_C"]),
                    "N_constraints": int(result["N_constraints"]),
                    "N_first": int(result["N_first"]),
                    "N_second": int(result["N_second"]),
                    "tertiary_count": int(result.get("tertiary_count", 0)),
                    "hessian_audit": audit,
                    "status": audit["status"],
                }
            )
        active = axis in active_axes
        expected_nphys = "1" if active else "0"
        status = "PASS" if all(row["N_phys"] == expected_nphys and row["status"] == "PASS" for row in point_rows) else "FAIL"
        rows[axis] = {
            "axis": axis,
            "projection": proj["projection"],
            "active_in_K_form": active,
            "expected_N_phys": expected_nphys,
            "basis_type": "q_cfg+p_cfg reduced scalar-axis phase-space",
            "q_cfg": [str(q) for q in proj["q_cfg"]],
            "point_rows": point_rows,
            "status": status,
        }
    return rows
