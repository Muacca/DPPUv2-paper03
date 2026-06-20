"""Dirac constraint analysis utilities.

Provides the fundamental building blocks of the Dirac constraint procedure:
Poisson bracket, coefficient-vector / rank utilities, and linear-independence
test.  These functions operate on plain SymPy expressions; callers that hold
``Constraint`` dataclass objects should extract ``.expr`` before passing to
``linearly_independent``.
"""

from __future__ import annotations

from typing import List, Sequence

import sympy as sp


def poisson_bracket(
    f: sp.Expr,
    g: sp.Expr,
    q_vars: Sequence[sp.Symbol],
    p_vars: Sequence[sp.Symbol],
) -> sp.Expr:
    """Standard Poisson bracket {f, g} = Σ_I (∂f/∂q_I ∂g/∂p_I − ∂f/∂p_I ∂g/∂q_I).

    Parameters
    ----------
    f, g    : SymPy expressions depending on ``q_vars`` and ``p_vars``.
    q_vars  : Ordered sequence of generalised coordinates.
    p_vars  : Corresponding canonical momenta (same length as ``q_vars``).

    Returns
    -------
    sp.Expr
        The Poisson bracket (not simplified; wrap in ``sp.simplify`` if needed).
    """
    result = sp.S.Zero
    for q, p in zip(q_vars, p_vars):
        result += sp.diff(f, q) * sp.diff(g, p) - sp.diff(f, p) * sp.diff(g, q)
    return result


def coeff_vector(expr: sp.Expr, vars_: Sequence[sp.Symbol]) -> List[sp.Expr]:
    """Coefficient vector of *expr* w.r.t. *vars_*, with constant term appended.

    Assumes *expr* is (at most) linear in *vars_*.
    """
    expr = sp.expand(expr)
    subs_zero = {var: 0 for var in vars_}
    coeffs = [sp.simplify(sp.diff(expr, var)) for var in vars_]
    const = sp.simplify(expr.subs(subs_zero))
    return coeffs + [const]


def rank_of_exprs(exprs: Sequence[sp.Expr], vars_: Sequence[sp.Symbol]) -> int:
    """Matrix rank of the coefficient-vector matrix for *exprs* w.r.t. *vars_*."""
    if not exprs:
        return 0
    return sp.Matrix([coeff_vector(expr, vars_) for expr in exprs]).rank()


def linearly_independent(
    expr: sp.Expr,
    existing: Sequence[sp.Expr],
    vars_: Sequence[sp.Symbol],
) -> bool:
    """Return ``True`` if *expr* is linearly independent from *existing*.

    Parameters
    ----------
    expr     : Candidate expression (linear in *vars_*).
    existing : Already-accepted constraint expressions (plain ``sp.Expr``).
    vars_    : Phase-space variables for coefficient extraction.
    """
    expr = sp.simplify(expr)
    if expr == 0:
        return False
    old_rank = rank_of_exprs(list(existing), vars_)
    new_rank = rank_of_exprs(list(existing) + [expr], vars_)
    return new_rank > old_rank
