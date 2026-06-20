"""Per-topology Branch-A local Lagrangian builders (DBR route).

Pure compute module — no I/O, no logging, no side effects on import.

Provides Branch-A first-order local Lagrangian builders for all four
homogeneous topologies (Nil3, S3, Sol3, T3), together with their
in-file helpers: nil3_structure_constants, s3_structure_constants,
sol3_structure_constants, lc_connection_from_structure_constants,
scalar_curvature_scaled, build_full_local_projection,
dirac_primary_nullspace, poisson, coeff_vector, is_independent.

The old builders received ``(nil3_mod, engine)`` / ``(engine, topology)``
module arguments; this module replaces those with direct imports from
``dppu.dirac.constraint_closure``.
"""

from __future__ import annotations

from itertools import combinations
from typing import Any, Dict, List, Sequence, Tuple

import sympy as sp

from dppu.dirac.constraint_closure import (
    build_L_conn_alg,
    eps3,
    lc4,
)


# ---------------------------------------------------------------------------
# Structure-constant helpers
# ---------------------------------------------------------------------------

def nil3_structure_constants(q: sp.Symbol) -> Dict[Tuple[int, int, int], sp.Expr]:
    # Spatial 1-based frame convention: C^3_{12}=-1/q.
    return {
        (3, 1, 2): -1 / q,
        (3, 2, 1): 1 / q,
    }


def s3_structure_constants(q: sp.Symbol) -> Dict[Tuple[int, int, int], sp.Expr]:
    return {
        (u, i, j): sp.simplify(2 * eps3(u, i, j) / q)
        for u in range(1, 4)
        for i in range(1, 4)
        for j in range(1, 4)
        if eps3(u, i, j) != 0
    }


def sol3_structure_constants(q: sp.Symbol) -> Dict[Tuple[int, int, int], sp.Expr]:
    return {
        (2, 1, 2): 1 / q,
        (2, 2, 1): -1 / q,
        (3, 1, 3): -1 / q,
        (3, 3, 1): 1 / q,
    }


def t3_structure_constants(q: sp.Symbol) -> Dict[Tuple[int, int, int], sp.Expr]:  # noqa: ARG001
    return {}


def _get_structure(topology: str, q: sp.Symbol) -> Dict[Tuple[int, int, int], sp.Expr]:
    if topology == "Nil3":
        return nil3_structure_constants(q)
    if topology == "S3":
        return s3_structure_constants(q)
    if topology == "Sol3":
        return sol3_structure_constants(q)
    if topology == "T3":
        return t3_structure_constants(q)
    raise ValueError(f"Unknown topology: {topology!r}")


# ---------------------------------------------------------------------------
# LC connection from structure constants
# ---------------------------------------------------------------------------

def lc_connection_from_structure_constants(
    q: sp.Symbol,
    structure: Dict[Tuple[int, int, int], sp.Expr],
) -> Dict[Tuple[int, int, int], sp.Expr]:
    # Convention: omega_lc(c,e,rho)=-Gamma^c_{rho,e}.
    def C(up: int, i: int, j: int) -> sp.Expr:
        return structure.get((up, i, j), sp.S.Zero)

    def gamma(up: int, i: int, j: int) -> sp.Expr:
        return sp.simplify((C(up, i, j) - C(i, j, up) + C(j, up, i)) / 2)

    return {
        (c, e, rho): sp.simplify(-gamma(c, rho, e))
        for c in range(1, 4)
        for e in range(1, 4)
        for rho in range(1, 4)
    }


def scalar_curvature_scaled(
    q: sp.Symbol,
    structure: Dict[Tuple[int, int, int], sp.Expr],
) -> sp.Expr:
    def C(up: int, i: int, j: int) -> sp.Expr:
        return structure.get((up, i, j), sp.S.Zero)

    def gamma(up: int, i: int, j: int) -> sp.Expr:
        return sp.simplify((C(up, i, j) - C(i, j, up) + C(j, up, i)) / 2)

    scalar = sp.S.Zero
    for upper in range(1, 4):
        for k in range(1, 4):
            for i in range(1, 4):
                for j in range(1, 4):
                    value = sum(
                        gamma(m, j, k) * gamma(upper, i, m)
                        - gamma(m, i, k) * gamma(upper, j, m)
                        for m in range(1, 4)
                    )
                    value -= sum(C(m, i, j) * gamma(upper, m, k) for m in range(1, 4))
                    if upper == i and k == j:
                        scalar += value
    return sp.factor(sp.cancel(scalar * q**2))


# ---------------------------------------------------------------------------
# Core Branch-A builder (Nil3 / S3 / Sol3 / T3 via topology dispatch)
# ---------------------------------------------------------------------------

def build_branchA_local_L2_with_lc(topology: str) -> Dict[str, Any]:
    """Build the Branch-A first-order quadratic local Lagrangian block for *topology*.

    Old signature was ``build_branchA_local_L2_with_lc(engine, topology)`` where
    ``engine`` was a module-object argument.  Here we replace those
    ``engine.eps3`` / ``engine.lc4`` / ``engine.build_L_conn_alg`` calls with
    direct imports from ``dppu.dirac.constraint_closure``.
    """
    N, q = sp.symbols("N q", positive=True)

    def e_bg(a: int, mu: int) -> sp.Expr:
        if a == 0 and mu == 0:
            return N
        if a == mu and a > 0:
            return q
        return sp.S.Zero

    s = sp.Symbol("s")
    n_pert = sp.Symbol("n_pert")
    NADM = [sp.Symbol(f"NADM_{i}") for i in range(1, 4)]
    qA00, qA11, qA01, qA02, qA12 = sp.symbols("qA_00 qA_11 qA_01 qA_02 qA_12")
    qA22 = -(qA00 + qA11)
    qA_comp = {
        (0, 0): qA00,
        (1, 1): qA11,
        (2, 2): qA22,
        (0, 1): qA01,
        (1, 0): qA01,
        (0, 2): qA02,
        (2, 0): qA02,
        (1, 2): qA12,
        (2, 1): qA12,
    }
    X01, X02, X12 = sp.symbols("X_01 X_02 X_12")
    X_comp = {
        (0, 1): X01,
        (1, 0): -X01,
        (0, 2): X02,
        (2, 0): -X02,
        (1, 2): X12,
        (2, 1): -X12,
        (0, 0): sp.S.Zero,
        (1, 1): sp.S.Zero,
        (2, 2): sp.S.Zero,
    }

    lorentz_pairs = list(combinations(range(4), 2))
    K0_syms = {(a, b): sp.Symbol(f"K0_{a}{b}") for a, b in lorentz_pairs}
    omega_syms = {
        (a, b, i): sp.Symbol(f"om_{a}{b}_{i}")
        for a, b in lorentz_pairs
        for i in range(1, 4)
    }
    dot_omega_pos_syms = {
        (a, b, i): sp.Symbol(f"domega_{a}{b}_{i}")
        for a, b in lorentz_pairs
        for i in range(1, 4)
    }
    dot_omega: Dict[Tuple[int, int, int], sp.Expr] = {}
    for a in range(4):
        for b in range(4):
            for i in range(1, 4):
                if a < b:
                    dot_omega[(a, b, i)] = dot_omega_pos_syms[(a, b, i)]
                elif a > b:
                    dot_omega[(a, b, i)] = -dot_omega_pos_syms[(b, a, i)]
                else:
                    dot_omega[(a, b, i)] = sp.S.Zero

    de_A: Dict[Tuple[int, int], sp.Expr] = {}
    for j in range(1, 4):
        de_A[(0, j)] = sp.S.Zero
    for a in range(1, 4):
        for j in range(1, 4):
            ia, ij = a - 1, j - 1
            val = qA_comp.get((ia, ij), sp.S.Zero) + X_comp.get((ia, ij), sp.S.Zero)
            if ia == ij:
                val += s / 3
            de_A[(a, j)] = val

    de_ext_A = dict(de_A)
    de_ext_A[(0, 0)] = n_pert
    for j in range(1, 4):
        de_ext_A[(j, 0)] = q * NADM[j - 1]

    structure = _get_structure(topology, q)
    omega_lc_map = lc_connection_from_structure_constants(q, structure)

    def omega_lc(c: int, e: int, rho: int) -> sp.Expr:
        if c == 0 or e == 0 or rho == 0:
            return sp.S.Zero
        return omega_lc_map.get((c, e, rho), sp.S.Zero)

    def K_sp_up(e: int, d: int, sigma: int) -> sp.Expr:
        if e == d or sigma == 0:
            return sp.S.Zero
        if e < d:
            return omega_syms.get((e, d, sigma), sp.S.Zero)
        return -omega_syms.get((d, e, sigma), sp.S.Zero)

    L_kin_A = sp.S.Zero
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    lc = lc4(a, b, c, d)
                    if lc == 0:
                        continue
                    for i in range(1, 4):
                        for j in range(1, 4):
                            for k in range(1, 4):
                                eijk = eps3(i, j, k)
                                if eijk == 0:
                                    continue
                                de_val = de_A.get((d, k), sp.S.Zero)
                                if de_val == sp.S.Zero:
                                    continue
                                L_kin_A += -2 * lc * eijk * e_bg(c, j) * de_val * dot_omega[(a, b, i)]
    L_kin_A = sp.expand(L_kin_A)

    K_build = {
        (a, b, mu): sp.Symbol(f"K_build_{a}{b}_{mu}")
        for a, b in lorentz_pairs
        for mu in range(4)
    }
    N0_b = sp.Symbol("N0_build", positive=True)
    q0_b = sp.Symbol("q0_build", positive=True)
    L_conn_alg_raw = build_L_conn_alg(N0_b, q0_b, K_build)
    K_subs: Dict[sp.Symbol, sp.Expr] = {N0_b: N + n_pert, q0_b: q}
    for a, b in lorentz_pairs:
        K_subs[K_build[(a, b, 0)]] = K0_syms[(a, b)]
        for i in range(1, 4):
            K_subs[K_build[(a, b, i)]] = omega_syms[(a, b, i)]
    L_conn_alg = sp.expand(L_conn_alg_raw.subs(K_subs))

    L_wK_A = sp.S.Zero
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    eps_abcd = lc4(a, b, c, d)
                    if eps_abcd == 0:
                        continue
                    for e in range(4):
                        for mu in range(4):
                            dea = de_ext_A.get((a, mu), sp.S.Zero)
                            if dea == sp.S.Zero:
                                continue
                            for nu in range(4):
                                ebg_b = e_bg(b, nu)
                                if ebg_b == sp.S.Zero:
                                    continue
                                for rho in range(4):
                                    olc = omega_lc(c, e, rho)
                                    if olc == sp.S.Zero:
                                        continue
                                    for sigma in range(4):
                                        Ksp = K_sp_up(e, d, sigma)
                                        if Ksp == sp.S.Zero:
                                            continue
                                        eps_pos = lc4(mu, nu, rho, sigma)
                                        if eps_pos == 0:
                                            continue
                                        L_wK_A += 2 * eps_abcd * eps_pos * dea * ebg_b * olc * Ksp
    L_wK_A = sp.expand(L_wK_A)

    L_de2_A = sp.S.Zero
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    eps_abcd = lc4(a, b, c, d)
                    if eps_abcd == 0:
                        continue
                    for e in range(4):
                        for mu in range(4):
                            dea = de_ext_A.get((a, mu), sp.S.Zero)
                            if dea == sp.S.Zero:
                                continue
                            for nu in range(4):
                                deb = de_ext_A.get((b, nu), sp.S.Zero)
                                if deb == sp.S.Zero:
                                    continue
                                for rho in range(4):
                                    olc1 = omega_lc(c, e, rho)
                                    if olc1 == sp.S.Zero:
                                        continue
                                    for sigma in range(4):
                                        olc2 = omega_lc(e, d, sigma)
                                        if olc2 == sp.S.Zero:
                                            continue
                                        eps_pos = lc4(mu, nu, rho, sigma)
                                        if eps_pos == 0:
                                            continue
                                        L_de2_A += eps_abcd * eps_pos * dea * deb * olc1 * olc2
    L_de2_A = sp.expand(L_de2_A)

    return {
        "topology": topology,
        "N": N,
        "q": q,
        "n_pert": n_pert,
        "NADM": NADM,
        "s": s,
        "qA_vars": [qA00, qA11, qA01, qA02, qA12],
        "X_vars": [X01, X02, X12],
        "lorentz_pairs": lorentz_pairs,
        "K0_syms": K0_syms,
        "omega_syms": omega_syms,
        "dot_omega_pos_syms": dot_omega_pos_syms,
        "structure": structure,
        "omega_lc_map": omega_lc_map,
        "scalar_curvature_scaled": scalar_curvature_scaled(q, structure),
        "L_kin_A": L_kin_A,
        "L_static_A": sp.expand(L_conn_alg + L_wK_A + L_de2_A),
        "L_total_A": sp.expand(L_kin_A + L_conn_alg + L_wK_A + L_de2_A),
    }


def build_t3_local() -> Dict[str, Any]:
    """Build the T3 (flat Bianchi-I) Branch-A local block.

    Old signature was ``build_t3_local(nil3_mod, engine)`` and patched
    the structure constants to return ``{}``.  Here we call
    ``build_branchA_local_L2_with_lc("T3")`` directly, which dispatches to
    the zero structure constant branch.
    """
    return build_branchA_local_L2_with_lc("T3")


def build_sol3_local() -> Dict[str, Any]:
    """Build the Sol3 Branch-A local block.

    Old signature was ``build_sol3_local(nil3_mod, engine)`` and patched
    the structure constants to return the Sol3 values.
    Here we call ``build_branchA_local_L2_with_lc("Sol3")`` directly.
    """
    return build_branchA_local_L2_with_lc("Sol3")


# ---------------------------------------------------------------------------
# Full local projection and Dirac primary nullspace (helper, kept for
# completeness — not used directly by the paper witnesses but available
# for callers of the library).
# ---------------------------------------------------------------------------

def build_full_local_projection(local: Dict[str, Any]) -> Dict[str, Any]:
    q_cfg: List[sp.Symbol] = []
    velocity_by_q: Dict[sp.Symbol, sp.Symbol] = {}
    sectors: Dict[sp.Symbol, str] = {}

    def add(var: sp.Symbol, sector: str) -> None:
        q_cfg.append(var)
        sectors[var] = sector

    add(local["n_pert"], "multiplier")
    for item in local["NADM"]:
        add(item, "multiplier")
    add(local["s"], "scalar")
    for item in local["qA_vars"]:
        add(item, "spin-2")
    for item in local["X_vars"]:
        add(item, "spin-1")
    for pair in local["lorentz_pairs"]:
        add(local["K0_syms"][pair], "multiplier")
    for a, b in local["lorentz_pairs"]:
        for i in range(1, 4):
            omega = local["omega_syms"][(a, b, i)]
            add(omega, "connection")
            velocity_by_q[omega] = local["dot_omega_pos_syms"][(a, b, i)]

    return {
        "q_cfg": q_cfg,
        "velocity_by_q": velocity_by_q,
        "sectors": sectors,
        "L_proj": local["L_total_A"],
        "L_static": local["L_static_A"],
        "L_kin": local["L_kin_A"],
    }
