"""
Effective Potential
===================

Computes V(r) = -S and provides numerical evaluation functions.
"""

from typing import Any, Callable, Dict, List, Optional
from sympy import S, cancel, expand, collect, lambdify


def compute_effective_potential(
    action: Any,
    logger: Optional[Any] = None
) -> Any:
    """
    Compute effective potential V(r) = -Action.

    Args:
        action: Total action S
        logger: Optional logger

    Returns:
        Effective potential V(r)
    """
    if logger:
        logger.info("Computing Effective Potential V(r) = -S...")

    V = cancel(-action)

    if logger:
        logger.info(f"  V(r) = {V}")

    return V


def get_potential_function(
    potential: Any,
    param_symbols: List[Any]
) -> Callable:
    """
    Create lambdified function for fast numerical evaluation.

    Args:
        potential: Symbolic potential V(r)
        param_symbols: List of symbols [r, V, eta, theta, L, kappa]

    Returns:
        Callable that accepts numpy arrays
    """
    return lambdify(param_symbols, potential, modules='numpy')


def decompose_potential(
    potential: Any,
    r: Any
) -> Dict[str, Any]:
    """
    Decompose potential by powers of r.

    Returns dict: {'r^3': coeff, 'r^2': coeff, 'r^1': coeff, ...}
    """
    expanded = expand(potential)
    collected = collect(expanded, r, evaluate=False)

    decomposition = {
        'r^3': S.Zero,
        'r^2': S.Zero,
        'r^1': S.Zero,
        'r^0': S.Zero,
        'total': potential,
    }

    for term, coeff in collected.items():
        if term == r**3:
            decomposition['r^3'] = coeff
        elif term == r**2:
            decomposition['r^2'] = coeff
        elif term == r:
            decomposition['r^1'] = coeff
        elif term == 1 or term == S.One:
            decomposition['r^0'] = coeff

    return decomposition
