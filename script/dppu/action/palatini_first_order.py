"""
First-Order Palatini Formulation (class-A minisuperspace)
==========================================================

Minimal first-order Palatini / Einstein-Hilbert density computation for
the four canonical class-A homogeneous three-topologies used in paper07:

    T3, Nil3, Sol3, S3

This module provides:

- ``eh_density``: the scale-adapted EH 4-form coefficient.
- ``palatini_reduced_density``: EH density with auxiliary connection
  fluctuations solved out.
- ``minisuperspace_dewitt_metric``: trace-free DeWitt metric on shape
  space.

**Convention warning**: The Levi-Civita connection computed here via
``_lc_connection`` uses the class-A Koszul formula::

    Gamma^up_{ij} from [e_i, e_j] = +C^k_{ij} e_k  (class-A)

This is distinct from the DPPU engine convention in
``dppu/engine/levi_civita.py`` (CONVENTIONS 3.2)::

    [E_b, E_c] = -C^a_{bc} E_a  (engine)

**Do not merge these two implementations.**  The engine convention
uses the opposite sign for structure constants.  Callers of
``_lc_connection`` must use ``structure_entries_class_a`` entries
(1-based, class-A sign convention).

Author: Muacca
"""

from __future__ import annotations

from collections import defaultdict
from itertools import combinations
from typing import Dict, Sequence, Tuple

import sympy as sp

from ..topology.class_a import structure_entries_class_a
from ..utils.symbolic import canonical_simplify

# ---------------------------------------------------------------------------
# Internal constants
# ---------------------------------------------------------------------------

_PAIR_LIST = list(combinations(range(4), 2))


# ---------------------------------------------------------------------------
# Private exterior algebra helpers
# (kept private to preserve exact expand-based normalization)
# ---------------------------------------------------------------------------

def _wedge_one(
    k1: Tuple[int, ...],
    k2: Tuple[int, ...],
) -> Tuple[int, Tuple[int, ...]]:
    merged = list(k1) + list(k2)
    if len(set(merged)) != len(merged):
        return 0, ()
    inversions = sum(
        1
        for i in range(len(merged))
        for j in range(i + 1, len(merged))
        if merged[i] > merged[j]
    )
    return (-1 if inversions % 2 else 1), tuple(sorted(merged))


def _form_add(
    *forms: Dict[Tuple[int, ...], sp.Expr],
) -> Dict[Tuple[int, ...], sp.Expr]:
    out: Dict[Tuple[int, ...], sp.Expr] = defaultdict(lambda: sp.S.Zero)
    for form in forms:
        for key, value in form.items():
            if value != 0:
                out[key] += value
    return {key: sp.expand(value) for key, value in out.items() if value != 0}


def _scalar_mul(
    value: sp.Expr,
    form: Dict[Tuple[int, ...], sp.Expr],
) -> Dict[Tuple[int, ...], sp.Expr]:
    if value == 0:
        return {}
    return {key: sp.expand(value * coeff) for key, coeff in form.items() if coeff != 0}


def _wedge(
    left: Dict[Tuple[int, ...], sp.Expr],
    right: Dict[Tuple[int, ...], sp.Expr],
) -> Dict[Tuple[int, ...], sp.Expr]:
    out: Dict[Tuple[int, ...], sp.Expr] = defaultdict(lambda: sp.S.Zero)
    for k1, v1 in left.items():
        for k2, v2 in right.items():
            sign, key = _wedge_one(k1, k2)
            if sign:
                out[key] += sign * v1 * v2
    return {key: sp.expand(value) for key, value in out.items() if value != 0}


def _basis_form(index: int) -> Dict[Tuple[int, ...], sp.Expr]:
    return {(index,): sp.S.One}


def _volume_coeff(form: Dict[Tuple[int, ...], sp.Expr]) -> sp.Expr:
    """Extract and expand the coefficient of e^{0123}."""
    return sp.expand(form.get((0, 1, 2, 3), sp.S.Zero))


# ---------------------------------------------------------------------------
# Scale-adapted structure and Levi-Civita connection
# ---------------------------------------------------------------------------

def _adapted_structure(
    topology: str,
    scales: Dict[int, sp.Expr],
) -> Dict[Tuple[int, int, int], sp.Expr]:
    """Scale-adapt class-A structure constants to the current frame.

    Each entry ``C^up_{ij}`` is rescaled by ``scale[up] / (scale[i] * scale[j])``
    to account for the diagonal frame deformation.

    Args:
        topology: Topology label.
        scales:   Dict ``{1: s1, 2: s2, 3: s3}`` of positive scale symbols.

    Returns:
        Rescaled structure-constant dict with 1-based keys.
    """
    out: Dict[Tuple[int, int, int], sp.Expr] = {}
    for (up, i, j), coeff in structure_entries_class_a(topology).items():
        out[(up, i, j)] = canonical_simplify(coeff * scales[up] / (scales[i] * scales[j]))
    return out


def _lc_connection(
    structure: Dict[Tuple[int, int, int], sp.Expr],
) -> Dict[Tuple[int, int, int], sp.Expr]:
    """Compute the Levi-Civita connection components in the class-A frame.

    Uses the Koszul formula for the class-A convention
    ``[e_i, e_j] = +C^k_{ij} e_k``::

        Gamma^a_{b rho} = -(1/2)(C^a_{rho b} - C^rho_{b a} + C^b_{a rho})

    **Sign-convention note**: This formula uses the class-A sign for structure
    constants and is NOT equivalent to the DPPU engine formula in
    ``dppu/engine/levi_civita.py``, which uses ``[E_b, E_c] = -C^a_{bc} E_a``
    (CONVENTIONS 3.2, opposite sign).

    Args:
        structure: Scale-adapted structure constants ``{(up, i, j): coeff}``
                   with 1-based indices.

    Returns:
        Connection dict ``{(a, b, rho): Gamma^a_{b rho}}``, 1-based.
    """
    def c(up: int, i: int, j: int) -> sp.Expr:
        return structure.get((up, i, j), sp.S.Zero)

    def gamma(up: int, i: int, j: int) -> sp.Expr:
        return canonical_simplify((c(up, i, j) - c(i, j, up) + c(j, up, i)) / 2)

    return {
        (a, b, rho): canonical_simplify(-gamma(a, rho, b))
        for a in range(1, 4)
        for b in range(1, 4)
        for rho in range(1, 4)
    }


# ---------------------------------------------------------------------------
# Maurer-Cartan exterior derivative (class-A structure)
# ---------------------------------------------------------------------------

def _d_basis(
    index: int,
    structure: Dict[Tuple[int, int, int], sp.Expr],
) -> Dict[Tuple[int, ...], sp.Expr]:
    if index == 0:
        return {}
    out: Dict[Tuple[int, ...], sp.Expr] = defaultdict(lambda: sp.S.Zero)
    for (up, i, j), coeff in structure.items():
        if up != index:
            continue
        sign, key = _wedge_one((i,), (j,))
        if sign:
            out[key] += -sp.Rational(1, 2) * coeff * sign
    return {key: sp.expand(value) for key, value in out.items() if value != 0}


def _d_form(
    form: Dict[Tuple[int, ...], sp.Expr],
    structure: Dict[Tuple[int, int, int], sp.Expr],
) -> Dict[Tuple[int, ...], sp.Expr]:
    out: Dict[Tuple[int, ...], sp.Expr] = {}
    for key, coeff in form.items():
        if len(key) != 1:
            continue
        out = _form_add(out, _scalar_mul(coeff, _d_basis(key[0], structure)))
    return out


# ---------------------------------------------------------------------------
# Connection 1-form and curvature 2-form
# ---------------------------------------------------------------------------

def _omega_form(
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


def _curvature_ab(
    a: int,
    b: int,
    lc: Dict[Tuple[int, int, int], sp.Expr],
    scales: Dict[int, sp.Expr],
    structure: Dict[Tuple[int, int, int], sp.Expr],
    deltas: Dict[Tuple[int, int, int], sp.Symbol] | None,
) -> Dict[Tuple[int, ...], sp.Expr]:
    result = _d_form(_omega_form(a, b, lc, scales, deltas), structure)
    for c in range(4):
        result = _form_add(
            result,
            _scalar_mul(
                sp.Integer(1 if c == 0 else -1),
                _wedge(
                    _omega_form(a, c, lc, scales, deltas),
                    _omega_form(c, b, lc, scales, deltas),
                ),
            ),
        )
    return result


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def eh_density(
    topology: str,
    scales: Dict[int, sp.Expr],
    deltas: Dict[Tuple[int, int, int], sp.Symbol] | None = None,
) -> sp.Expr:
    """First-order EH density coefficient for the finite Palatini branch witness.

    Computes the coefficient of the 4-volume form in the first-order
    Palatini action, with connection deformations *deltas* (if provided)
    treated as independent variables.

    Args:
        topology: Topology label (``"T3"``, ``"Nil3"``, ``"Sol3"``, ``"S3"``).
        scales:   Dict ``{1: s1, 2: s2, 3: s3}`` of positive scale symbols
                  (time-direction scale is implicit and fixed to 1 via the
                  e^0 basis form).
        deltas:   Optional connection fluctuation symbols
                  ``{(a, b, rho): D_symbol}`` with 1-based frame indices.
                  When ``None`` the background Levi-Civita connection is used.

    Returns:
        SymPy scalar expression: the volume coefficient of the EH 4-form.
    """
    structure = _adapted_structure(topology, scales)
    lc = _lc_connection(structure)
    e_forms = {
        0: _basis_form(0),
        1: _scalar_mul(scales[1], _basis_form(1)),
        2: _scalar_mul(scales[2], _basis_form(2)),
        3: _scalar_mul(scales[3], _basis_form(3)),
    }
    curv = {
        (a, b): _curvature_ab(a, b, lc, scales, structure, deltas)
        for a, b in _PAIR_LIST
    }
    eh_form: Dict[Tuple[int, ...], sp.Expr] = {}
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    eps = int(sp.LeviCivita(a, b, c, d))
                    if eps == 0 or c == d:
                        continue
                    lo, hi = sorted((c, d))
                    sign = sp.S.One if c < d else -sp.S.One
                    term = _scalar_mul(
                        eps * sign,
                        _wedge(_wedge(e_forms[a], e_forms[b]), curv[(lo, hi)]),
                    )
                    eh_form = _form_add(eh_form, term)
    return sp.expand(_volume_coeff(eh_form))


def palatini_reduced_density(
    topology: str,
    vector: Sequence[sp.Expr],
) -> sp.Expr:
    """Solve connection fluctuations and return reduced EH density along a direction.

    Introduces connection fluctuation symbols, solves their equations of
    motion (linear in the fluctuations), substitutes the solution back,
    and returns the on-shell EH density as a function of the scale
    parameters.

    Args:
        topology: Topology label.
        vector:   Direction vector ``[v1, v2, v3]`` in log-scale space.
                  The scale is set to ``s_i = exp(z * v_i)``.

    Returns:
        Reduced EH density as a SymPy expression in ``z``, or ``sp.nan``
        if the fluctuation matrix is rank-deficient.
    """
    z = sp.Symbol("z")
    scales = {i: sp.exp(z * vector[i - 1]) for i in (1, 2, 3)}
    deltas = {
        (a, b, rho): sp.Symbol(f"D_{a}{b}_{rho}")
        for a, b in _PAIR_LIST
        for rho in range(1, 4)
    }
    delta_vars = [deltas[(a, b, rho)] for a, b in _PAIR_LIST for rho in range(1, 4)]
    density = sp.expand(eh_density(topology, scales, deltas))
    equations = [sp.diff(density, var) for var in delta_vars]
    mat, vec = sp.linear_eq_to_matrix(equations, delta_vars)
    if int(mat.rank()) != len(delta_vars) or int(mat.row_join(vec).rank()) != int(mat.rank()):
        return sp.nan
    solution = mat.LUsolve(vec)
    return canonical_simplify(density.subs(dict(zip(delta_vars, solution))))


def minisuperspace_dewitt_metric() -> sp.Matrix:
    """Trace-free DeWitt metric on shape space (homogeneous minisuperspace model).

    Returns the diagonal metric ``G = diag(s2*s3, s1*s3, s1*s2)`` on
    the three-dimensional space of scale factors (s1, s2, s3), which
    corresponds to the trace-free sector of the DeWitt metric for the
    minisuperspace Bianchi class-A model.

    Returns:
        3×3 SymPy diagonal matrix.
    """
    s1, s2, s3 = sp.symbols("s1 s2 s3", positive=True, nonzero=True)
    return sp.diag(s2 * s3, s1 * s3, s1 * s2)
