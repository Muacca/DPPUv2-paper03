"""Scalar-axis projection and Dirac standard reduction for Branch-A locals.

Pure compute module — no I/O, no logging, no side effects on import.

Provides:
  — build_scalar_axis_projection
  — dirac_standard
  — reduced_hamiltonian_status
  — (helpers: poisson, coeff_vector, is_independent)

The helper functions poisson / coeff_vector / is_independent are provided as
module-level functions; they are NOT re-imported from
``dppu.dirac.constraint_closure`` because the constraint_closure versions use
``sp.simplify(result)`` for poisson while these versions use ``sp.expand(...)``
— the implementations are numerically equivalent but differ in intermediate
form.  To stay behavior-preserving (diff-zero gate), the expand-based
implementations are kept here verbatim.
"""

from __future__ import annotations

from typing import Any, Dict, List, Sequence, Tuple

import sympy as sp


# ---------------------------------------------------------------------------
# Local helpers (expand-based variants; see module docstring for rationale)
# ---------------------------------------------------------------------------

def poisson(f: sp.Expr, g: sp.Expr, q_vars: Sequence[sp.Symbol], p_vars: Sequence[sp.Symbol]) -> sp.Expr:
    return sp.expand(
        sum(sp.diff(f, q) * sp.diff(g, p) - sp.diff(f, p) * sp.diff(g, q) for q, p in zip(q_vars, p_vars))
    )


def coeff_vector(expr: sp.Expr, vars_all: Sequence[sp.Symbol]) -> List[sp.Expr]:
    expr = sp.expand(expr)
    return [sp.simplify(sp.diff(expr, var)) for var in vars_all]


def is_independent(expr: sp.Expr, existing: Sequence[sp.Expr], vars_all: Sequence[sp.Symbol]) -> bool:
    expr = sp.expand(expr)
    if expr == sp.S.Zero:
        return False
    old_rows = [coeff_vector(item, vars_all) for item in existing]
    new_rows = old_rows + [coeff_vector(expr, vars_all)]
    old_rank = sp.Matrix(old_rows).rank() if old_rows else 0
    return sp.Matrix(new_rows).rank() > old_rank


# ---------------------------------------------------------------------------
# Main public functions
# ---------------------------------------------------------------------------

def build_scalar_axis_projection(local: Dict[str, Any], axis: str) -> Dict[str, Any]:
    axis_index = {"x": 0, "y": 1, "z": 2}[axis]
    k = sp.Symbol(f"k_{axis}", positive=True)
    n_var, b_var, s_var, e_var, x_var = sp.symbols(f"nS_{axis} BS_{axis} sS_{axis} ES_{axis} XS_{axis}")
    lorentz_pairs: List[Tuple[int, int]] = list(local["lorentz_pairs"])
    k0 = {(a, b): sp.Symbol(f"K0S_{axis}_{a}{b}") for a, b in lorentz_pairs}
    omega = {(a, b): sp.Symbol(f"omS_{axis}_{a}{b}") for a, b in lorentz_pairs}
    domega = {(a, b): sp.Symbol(f"domS_{axis}_{a}{b}") for a, b in lorentz_pairs}

    q_cfg: List[sp.Symbol] = [n_var, b_var, s_var, e_var] + [k0[pair] for pair in lorentz_pairs] + [x_var] + [omega[pair] for pair in lorentz_pairs]
    velocity_by_q = {omega[pair]: domega[pair] for pair in lorentz_pairs}
    sectors: Dict[sp.Symbol, str] = {n_var: "multiplier", b_var: "scalar", s_var: "scalar", e_var: "scalar", x_var: "spin-1"}
    for pair in lorentz_pairs:
        sectors[k0[pair]] = "multiplier"
        sectors[omega[pair]] = "connection"

    direction = [sp.S.Zero, sp.S.Zero, sp.S.Zero]
    direction[axis_index] = sp.S.One
    norm = sp.sqrt(sp.Rational(1, 6))
    h_diag = [(-2 if idx == axis_index else 1) * norm for idx in range(3)]

    subs: Dict[sp.Symbol, sp.Expr] = {local["n_pert"]: n_var, local["s"]: s_var}
    for idx in range(3):
        subs[local["NADM"][idx]] = k * b_var * direction[idx]

    qA00, qA11, qA01, qA02, qA12 = local["qA_vars"]
    subs.update({qA00: e_var * h_diag[0], qA11: e_var * h_diag[1], qA01: sp.S.Zero, qA02: sp.S.Zero, qA12: sp.S.Zero})

    X01, X02, X12 = local["X_vars"]
    for idx, frame_var in enumerate([X01, X02, X12]):
        subs[frame_var] = k * x_var * direction[idx]

    for a, b in lorentz_pairs:
        subs[local["K0_syms"][(a, b)]] = k0[(a, b)]
        for i in range(1, 4):
            subs[local["omega_syms"][(a, b, i)]] = k * omega[(a, b)] * direction[i - 1]
            subs[local["dot_omega_pos_syms"][(a, b, i)]] = k * domega[(a, b)] * direction[i - 1]

    return {
        "case_id": f"Nil3_scalar_axis_{axis}",
        "sector": "scalar",
        "projection": f"nonzero_frame_k_parallel_{axis}",
        "axis": axis,
        "k": k,
        "q_cfg": q_cfg,
        "velocity_by_q": velocity_by_q,
        "sectors": sectors,
        "L_proj": sp.expand(local["L_total_A"].subs(subs)),
        "L_static": sp.expand(local["L_static_A"].subs(subs)),
        "L_kin": sp.expand(local["L_kin_A"].subs(subs)),
    }


def dirac_standard(proj: Dict[str, Any], local: Dict[str, Any], rank_subs: Dict[sp.Symbol, sp.Expr]) -> Dict[str, Any]:
    q_cfg: List[sp.Symbol] = proj["q_cfg"]
    p_cfg = [sp.Symbol(f"p_{q}") for q in q_cfg]
    p_by_q = dict(zip(q_cfg, p_cfg))
    vars_all = q_cfg + p_cfg

    primaries: List[sp.Expr] = []
    for q_var in q_cfg:
        vel = proj["velocity_by_q"].get(q_var)
        momentum = sp.diff(proj["L_proj"], vel) if vel is not None else sp.S.Zero
        primaries.append(sp.expand(p_by_q[q_var] - momentum).subs(rank_subs))

    h_c = sp.expand((-proj["L_static"]).subs(rank_subs))
    c_primary = sp.Matrix([[poisson(a, b, q_cfg, p_cfg) for b in primaries] for a in primaries])
    h_brackets = [sp.expand(poisson(phi, h_c, q_cfg, p_cfg)) for phi in primaries]

    secondaries: List[sp.Expr] = []
    for null_vec in c_primary.nullspace():
        sec = sp.expand(sum(null_vec[i] * h_brackets[i] for i in range(len(q_cfg))))
        if is_independent(sec, primaries + secondaries, vars_all):
            secondaries.append(sec)

    constraints = primaries + secondaries
    c_full = sp.Matrix([[poisson(a, b, q_cfg, p_cfg) for b in constraints] for a in constraints])
    full_h_brackets = [sp.expand(poisson(c, h_c, q_cfg, p_cfg)) for c in constraints]
    tertiary: List[sp.Expr] = []
    for null_vec in c_full.nullspace():
        condition = sp.expand(sum(null_vec[i] * full_h_brackets[i] for i in range(len(constraints))))
        if is_independent(condition, constraints + tertiary, vars_all):
            tertiary.append(condition)

    n_config = len(q_cfg)
    n_constraints = len(constraints)
    rank_c = int(c_full.rank())
    n_first = n_constraints - rank_c
    n_second = rank_c
    n_phys = sp.Rational(2 * n_config - 2 * n_first - n_second, 2)
    return {
        "q_cfg": q_cfg,
        "p_cfg": p_cfg,
        "constraints": constraints,
        "H_c": h_c,
        "N_config": n_config,
        "N_phase": 2 * n_config,
        "N_primary": len(primaries),
        "N_secondary": len(secondaries),
        "N_constraints": n_constraints,
        "rank_C": rank_c,
        "kernel_dim": n_constraints - rank_c,
        "N_first": n_first,
        "N_second": n_second,
        "N_phys": n_phys,
        "tertiary_count": len(tertiary),
    }


def reduced_hamiltonian_status(result: Dict[str, Any]) -> Dict[str, str]:
    q_cfg: List[sp.Symbol] = result["q_cfg"]
    p_cfg: List[sp.Symbol] = result["p_cfg"]
    vars_all = q_cfg + p_cfg
    n_config = len(q_cfg)
    if result["N_phys"] == 0:
        return {
            "physical_pair_count": "0",
            "kinetic_or_hamiltonian_test": "N/A_NO_PHYSICAL_PAIR",
            "eigenvalues_or_signature": "N/A",
            "ghost_status": "N/A_NO_PHYSICAL_PAIR",
            "status": "PASS",
            "notes": "DBR leaves no physical pair in this representative.",
        }

    constraint_matrix = sp.Matrix([[sp.diff(c, var) for var in vars_all] for c in result["constraints"]])
    null_cols = constraint_matrix.nullspace()
    if not null_cols:
        return {
            "physical_pair_count": str(result["N_phys"]),
            "kinetic_or_hamiltonian_test": "CONSTRAINT_SURFACE_NULLSPACE_EMPTY",
            "eigenvalues_or_signature": "UNRESOLVED",
            "ghost_status": "UNRESOLVED",
            "status": "REVIEW_REQUIRED",
            "notes": "No constraint-surface basis was available.",
        }

    basis = sp.Matrix.hstack(*null_cols)
    symplectic = sp.zeros(2 * n_config, 2 * n_config)
    for idx in range(n_config):
        symplectic[idx, n_config + idx] = 1
        symplectic[n_config + idx, idx] = -1
    omega_reduced = basis.T * symplectic * basis
    physical_rank = int(omega_reduced.rank())
    hessian = sp.Matrix([[sp.diff(result["H_c"], a, b) for b in vars_all] for a in vars_all]) / 2

    nondegenerate_signatures: List[str] = []
    for i in range(basis.shape[1]):
        for j in range(i + 1, basis.shape[1]):
            omega_pair = sp.Matrix([[omega_reduced[i, i], omega_reduced[i, j]], [omega_reduced[j, i], omega_reduced[j, j]]])
            if omega_pair.det() == 0:
                continue
            sub_basis = basis[:, [i, j]]
            h_pair = sp.simplify(sub_basis.T * hessian * sub_basis)
            if h_pair == sp.zeros(2, 2):
                nondegenerate_signatures.append("ZERO_QUADRATIC_HAMILTONIAN")
            else:
                nondegenerate_signatures.append("NONZERO_QUADRATIC_HAMILTONIAN:" + sp.sstr(h_pair))

    if physical_rank == 2 and nondegenerate_signatures and all(item == "ZERO_QUADRATIC_HAMILTONIAN" for item in nondegenerate_signatures):
        return {
            "physical_pair_count": str(result["N_phys"]),
            "kinetic_or_hamiltonian_test": "constraint_surface_symplectic_reduction",
            "eigenvalues_or_signature": "reduced_quadratic_hamiltonian_zero",
            "ghost_status": "NO_NEGATIVE_DIRECTION_DETECTED_DEGENERATE_ZERO_MODE",
            "status": "PASS",
            "notes": "The scalar pair is a degenerate zero-mode at quadratic order, not a propagating ghost in this representative.",
        }

    return {
        "physical_pair_count": str(result["N_phys"]),
        "kinetic_or_hamiltonian_test": "constraint_surface_symplectic_reduction",
        "eigenvalues_or_signature": ";".join(nondegenerate_signatures) if nondegenerate_signatures else "UNRESOLVED",
        "ghost_status": "UNRESOLVED_NONZERO_REDUCED_HAMILTONIAN",
        "status": "REVIEW_REQUIRED",
        "notes": f"physical symplectic rank={physical_rank}",
    }
