"""Finite-amplitude background-adapted Palatini reduction engine.

Background-adapted finite-amplitude Palatini reduction: symbolic series
expansion of the EC density along homogeneous Bianchi paths with connection
perturbations, Routh elimination of connection auxiliaries, and coefficient
extraction at each order in the deformation parameter.

Only scaffolding (dir creation, logging, file I/O, main()) has been stripped;
all pure-compute functions and constants are preserved exactly.
"""

from __future__ import annotations

import time
from collections import defaultdict
from itertools import combinations
from typing import Any, Dict, Tuple

import sympy as sp


ETA = {0: -1, 1: 1, 2: 1, 3: 1}
PAIR_LIST = list(combinations(range(4), 2))


def wedge_one(k1: Tuple[int, ...], k2: Tuple[int, ...]) -> Tuple[int, Tuple[int, ...]]:
    merged = list(k1) + list(k2)
    if len(set(merged)) != len(merged):
        return 0, ()
    inversions = 0
    for i in range(len(merged)):
        for j in range(i + 1, len(merged)):
            if merged[i] > merged[j]:
                inversions += 1
    return (-1 if inversions % 2 else 1), tuple(sorted(merged))


def form_add(*forms: Dict[Tuple[int, ...], sp.Expr]) -> Dict[Tuple[int, ...], sp.Expr]:
    out: Dict[Tuple[int, ...], sp.Expr] = defaultdict(lambda: sp.S.Zero)
    for form in forms:
        for key, value in form.items():
            if value != 0:
                out[key] += value
    return {key: sp.expand(value) for key, value in out.items() if value != 0}


def scalar_mul(value: sp.Expr, form: Dict[Tuple[int, ...], sp.Expr]) -> Dict[Tuple[int, ...], sp.Expr]:
    if value == 0:
        return {}
    return {key: sp.expand(value * coeff) for key, coeff in form.items() if coeff != 0}


def wedge(f: Dict[Tuple[int, ...], sp.Expr], g: Dict[Tuple[int, ...], sp.Expr]) -> Dict[Tuple[int, ...], sp.Expr]:
    out: Dict[Tuple[int, ...], sp.Expr] = defaultdict(lambda: sp.S.Zero)
    for k1, v1 in f.items():
        for k2, v2 in g.items():
            sign, key = wedge_one(k1, k2)
            if sign:
                out[key] += sign * v1 * v2
    return {key: sp.expand(value) for key, value in out.items() if value != 0}


def basis_form(index: int) -> Dict[Tuple[int, ...], sp.Expr]:
    return {(index,): sp.S.One}


def volume_coeff(form: Dict[Tuple[int, ...], sp.Expr]) -> sp.Expr:
    return sp.expand(form.get((0, 1, 2, 3), sp.S.Zero))


def eps4(a: int, b: int, c: int, d: int) -> int:
    return int(sp.LeviCivita(a, b, c, d))


def trunc_exp(coef: sp.Expr, z: sp.Symbol, order: int = 4) -> sp.Expr:
    return sp.expand(sum((coef * z) ** n / sp.factorial(n) for n in range(order + 1)))


def trunc_poly(expr: sp.Expr, z: sp.Symbol, order: int = 4) -> sp.Expr:
    return sp.expand(sp.series(sp.expand(expr), z, 0, order + 1).removeO())


def base_structure(topology: str) -> Dict[Tuple[int, int, int], sp.Expr]:
    if topology == "T3":
        return {}
    if topology == "Nil3":
        return {(3, 1, 2): sp.S.One, (3, 2, 1): -sp.S.One}
    if topology == "Sol3":
        return {
            (2, 1, 2): sp.S.One,
            (2, 2, 1): -sp.S.One,
            (3, 1, 3): -sp.S.One,
            (3, 3, 1): sp.S.One,
        }
    raise ValueError(f"unknown topology: {topology}")


def adapted_structure(topology: str, scales: Dict[int, sp.Expr]) -> Dict[Tuple[int, int, int], sp.Expr]:
    out: Dict[Tuple[int, int, int], sp.Expr] = {}
    for (up, i, j), coeff in base_structure(topology).items():
        out[(up, i, j)] = sp.simplify(coeff * scales[up] / (scales[i] * scales[j]))
    return out


def lc_connection(structure: Dict[Tuple[int, int, int], sp.Expr]) -> Dict[Tuple[int, int, int], sp.Expr]:
    def c(up: int, i: int, j: int) -> sp.Expr:
        return structure.get((up, i, j), sp.S.Zero)

    def gamma(up: int, i: int, j: int) -> sp.Expr:
        return sp.simplify((c(up, i, j) - c(i, j, up) + c(j, up, i)) / 2)

    return {
        (a, b, rho): sp.simplify(-gamma(a, rho, b))
        for a in range(1, 4)
        for b in range(1, 4)
        for rho in range(1, 4)
    }


def d_basis(index: int, structure: Dict[Tuple[int, int, int], sp.Expr]) -> Dict[Tuple[int, ...], sp.Expr]:
    if index == 0:
        return {}
    out: Dict[Tuple[int, ...], sp.Expr] = defaultdict(lambda: sp.S.Zero)
    for (up, i, j), coeff in structure.items():
        if up != index:
            continue
        sign, key = wedge_one((i,), (j,))
        if sign:
            out[key] += -sp.Rational(1, 2) * coeff * sign
    return {key: sp.expand(value) for key, value in out.items() if value != 0}


def d_form(form: Dict[Tuple[int, ...], sp.Expr], structure: Dict[Tuple[int, int, int], sp.Expr]) -> Dict[Tuple[int, ...], sp.Expr]:
    out: Dict[Tuple[int, ...], sp.Expr] = {}
    for key, coeff in form.items():
        if len(key) != 1:
            continue
        out = form_add(out, scalar_mul(coeff, d_basis(key[0], structure)))
    return out


def path_scales(path: str, z: sp.Symbol | None, t_value: sp.Rational | None, route: str = "T") -> Dict[int, sp.Expr]:
    exponents = {
        "nil3_betaminus": {1: sp.Rational(1, 2), 2: sp.Rational(-1, 2), 3: sp.S.Zero},
        # Volume-preserving vertical-vs-horizontal Nil3 candidate; beta_plus = 3 z.
        "nil3_betaplus_u": {1: sp.Rational(-1, 2), 2: sp.Rational(-1, 2), 3: sp.S.One},
        # Sol3 beta_x direction (scalar-axis high-order probes).
        "sol3_betax": {1: sp.Rational(1, 3), 2: sp.Rational(-1, 6), 3: sp.Rational(-1, 6)},
        "t3_shape": {1: sp.Rational(1, 2), 2: sp.Rational(-1, 2), 3: sp.S.Zero},
    }[path]
    scales: Dict[int, sp.Expr] = {}
    for i, coef in exponents.items():
        if route == "F":
            if z is None:
                assert t_value is not None
                scales[i] = sp.exp(coef * t_value)
            else:
                scales[i] = sp.exp(coef * z)
        else:
            if z is None:
                assert t_value is not None
                scales[i] = trunc_exp(coef, sp.Symbol("_zz")).subs(sp.Symbol("_zz"), t_value)
            else:
                scales[i] = trunc_exp(coef, z)
    return scales


def omega_form(
    a: int,
    b: int,
    lc: Dict[Tuple[int, int, int], sp.Expr],
    scales: Dict[int, sp.Expr],
    deltas: Dict[Tuple[int, int, int], sp.Symbol] | None,
) -> Dict[Tuple[int, ...], sp.Expr]:
    if a == b:
        return {}
    sign_ab = sp.S.One if a < b else -sp.S.One
    lo, hi = sorted((a, b))
    out: Dict[Tuple[int, ...], sp.Expr] = {}
    for rho in range(1, 4):
        bg = sp.S.Zero
        if a in (1, 2, 3) and b in (1, 2, 3):
            bg = lc.get((a, b, rho), sp.S.Zero)
        delta = sp.S.Zero if deltas is None else sign_ab * deltas[(lo, hi, rho)]
        coeff = sp.expand((bg + delta) * scales[rho])
        if coeff != 0:
            out[(rho,)] = coeff
    return out


def curvature_ab(
    a: int,
    b: int,
    lc: Dict[Tuple[int, int, int], sp.Expr],
    scales: Dict[int, sp.Expr],
    structure: Dict[Tuple[int, int, int], sp.Expr],
    deltas: Dict[Tuple[int, int, int], sp.Symbol] | None,
) -> Dict[Tuple[int, ...], sp.Expr]:
    result = d_form(omega_form(a, b, lc, scales, deltas), structure)
    for c in range(4):
        result = form_add(result, scalar_mul(ETA[c], wedge(omega_form(a, c, lc, scales, deltas), omega_form(c, b, lc, scales, deltas))))
    return result


def eh_density(topology: str, scales: Dict[int, sp.Expr], deltas: Dict[Tuple[int, int, int], sp.Symbol] | None = None) -> sp.Expr:
    structure = base_structure(topology)
    adapted = adapted_structure(topology, scales)
    lc = lc_connection(adapted)
    e_forms = {
        0: basis_form(0),
        1: scalar_mul(scales[1], basis_form(1)),
        2: scalar_mul(scales[2], basis_form(2)),
        3: scalar_mul(scales[3], basis_form(3)),
    }
    curv = {(a, b): curvature_ab(a, b, lc, scales, structure, deltas) for a, b in PAIR_LIST}
    eh_form: Dict[Tuple[int, ...], sp.Expr] = {}
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    eps = eps4(a, b, c, d)
                    if eps == 0 or c == d:
                        continue
                    lo, hi = sorted((c, d))
                    sign = sp.S.One if c < d else -sp.S.One
                    term = scalar_mul(eps * sign, wedge(wedge(e_forms[a], e_forms[b]), curv[(lo, hi)]))
                    eh_form = form_add(eh_form, term)
    return sp.expand(volume_coeff(eh_form))


def palatini_symbolic_series(topology: str, path: str, route: str = "T") -> Dict[str, Any]:
    z = sp.Symbol("z")
    start = time.perf_counter()
    scales = path_scales(path, z, None, route=route)
    deltas = {(a, b, rho): sp.Symbol(f"D_{a}{b}_{rho}") for a, b in PAIR_LIST for rho in range(1, 4)}
    delta_vars = [deltas[(a, b, rho)] for a, b in PAIR_LIST for rho in range(1, 4)]
    density = sp.expand(eh_density(topology, scales, deltas))
    equations = [sp.diff(density, var) for var in delta_vars]
    mat, vec = sp.linear_eq_to_matrix(equations, delta_vars)
    rank = int(mat.rank())
    aug_rank = int(mat.row_join(vec).rank())
    if rank != len(delta_vars) or aug_rank != rank:
        return {
            "status": "INTRACTABLE_OR_DEGENERATE",
            "rank": rank,
            "aug_rank": aug_rank,
            "unknown_count": len(delta_vars),
            "elapsed_seconds": time.perf_counter() - start,
        }
    sol_vec = mat.LUsolve(vec)
    solution = {var: value for var, value in zip(delta_vars, sol_vec)}
    reduced = sp.cancel(density.subs(solution))
    delta4 = sp.factor(trunc_poly(reduced - reduced.subs(z, 0), z, 4))
    expanded = sp.expand(delta4)
    coeff = {str(order): sp.sstr(sp.expand(expanded.coeff(z, order))) for order in range(5)}
    nonzero_orders = [order for order, value in coeff.items() if order != "0" and sp.sympify(value) != 0]
    return {
        "status": "PASS",
        "route": route,
        "rank": rank,
        "aug_rank": aug_rank,
        "unknown_count": len(delta_vars),
        "reduced_density": sp.sstr(reduced),
        "delta_through_order4": sp.sstr(delta4),
        "coefficients": coeff,
        "nonzero_orders_1_to_4": nonzero_orders,
        "elapsed_seconds": time.perf_counter() - start,
    }
