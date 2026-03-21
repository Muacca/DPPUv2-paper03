"""
Levi-Civita Symbols
===================

Provides Levi-Civita (totally antisymmetric) symbols for 3D and 4D.

Conventions:
- 3D: ε_{012} = +1 (right-handed)
- 4D: ε_{0123} = +1 (Euclidean signature convention)

Author: Muacca
"""

from typing import Tuple


def epsilon_symbol(i: int, j: int, k: int) -> int:
    """
    3D Levi-Civita symbol for indices 0, 1, 2.

    Definition:
        ε_{ijk} = +1 for even permutations of (0,1,2)
        ε_{ijk} = -1 for odd permutations of (0,1,2)
        ε_{ijk} = 0  if any two indices are equal

    Args:
        i, j, k: Indices in {0, 1, 2}

    Returns:
        +1, -1, or 0

    Examples:
        >>> epsilon_symbol(0, 1, 2)
        1
        >>> epsilon_symbol(2, 1, 0)
        -1
        >>> epsilon_symbol(0, 0, 1)
        0
    """
    if (i, j, k) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
        return 1
    elif (i, j, k) in [(2, 1, 0), (0, 2, 1), (1, 0, 2)]:
        return -1
    else:
        return 0


def levi_civita_4d(mu: int, nu: int, rho: int, sigma: int) -> int:
    """
    4D Levi-Civita symbol for indices 0, 1, 2, 3.

    Definition:
        ε_{μνρσ} = sign of permutation bringing (μ,ν,ρ,σ) to (0,1,2,3)
        ε_{μνρσ} = 0 if any two indices are equal

    Convention:
        ε_{0123} = +1 (Euclidean signature)

    Args:
        mu, nu, rho, sigma: Indices in {0, 1, 2, 3}

    Returns:
        +1, -1, or 0

    Examples:
        >>> levi_civita_4d(0, 1, 2, 3)
        1
        >>> levi_civita_4d(1, 0, 2, 3)
        -1
        >>> levi_civita_4d(0, 0, 2, 3)
        0
    """
    indices = [mu, nu, rho, sigma]

    # Check for repeated indices
    if len(set(indices)) != 4:
        return 0

    # Count inversions (bubble sort)
    inversions = 0
    for i in range(4):
        for j in range(i + 1, 4):
            if indices[i] > indices[j]:
                inversions += 1

    return 1 if inversions % 2 == 0 else -1


def levi_civita_nd(indices: Tuple[int, ...]) -> int:
    """
    N-dimensional Levi-Civita symbol.

    Generalization to arbitrary dimension.

    Args:
        indices: Tuple of n indices, each in {0, 1, ..., n-1}

    Returns:
        +1, -1, or 0

    Examples:
        >>> levi_civita_nd((0, 1))
        1
        >>> levi_civita_nd((0, 1, 2, 3, 4))
        1
    """
    n = len(indices)

    # Check for repeated indices
    if len(set(indices)) != n:
        return 0

    # Count inversions
    inversions = 0
    for i in range(n):
        for j in range(i + 1, n):
            if indices[i] > indices[j]:
                inversions += 1

    return 1 if inversions % 2 == 0 else -1
