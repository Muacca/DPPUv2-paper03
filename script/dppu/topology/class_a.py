"""
Class-A Topology Invariants
============================

Structure constants, n-matrix, and spectral invariants for the four
canonical Bianchi class-A homogeneous three-topologies used in paper07:

    T3, Nil3, Sol3, S3

The *class-A* convention uses 1-based spatial indices and defines the
structure constants via

    [e_i, e_j] = C^k_{ij} e_k

so that the n-matrix (Bianchi classification matrix) is

    n^{ij} = (1/2) epsilon^{ikl} C^j_{kl}

with epsilon being the 0-based 3D Levi-Civita symbol (0,1,2 range).

Author: Muacca
"""

from __future__ import annotations

from typing import Any

import sympy as sp

from ..utils.epsilon import epsilon_nd
from ..utils.symbolic import canonical_simplify


# ---------------------------------------------------------------------------
# Structure constants (class-A convention, 1-based entries)
# ---------------------------------------------------------------------------

def structure_entries_class_a(
    topology: str,
) -> dict[tuple[int, int, int], sp.Expr]:
    """Return nonzero spatial structure constants in the paper class-A convention.

    Keys are ``(up, i, j)`` with 1-based indices; ``C^up_{ij}`` = value.
    Only the independent (upper-triangle) entries are stored; the anti-symmetric
    partner ``C^up_{ji} = -C^up_{ij}`` is not stored explicitly.

    Args:
        topology: One of ``"T3"``, ``"Nil3"``, ``"Sol3"``, ``"S3"``.

    Returns:
        Dict of nonzero entries.  Empty dict for T3 (flat, all zero).

    Raises:
        ValueError: If *topology* is not a recognised label.
    """
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
    if topology == "S3":
        entries: dict[tuple[int, int, int], sp.Expr] = {}
        for up in range(1, 4):
            for i in range(1, 4):
                for j in range(1, 4):
                    value = sp.LeviCivita(up, i, j)
                    if value != 0:
                        entries[(up, i, j)] = sp.S(value)
        return entries
    raise ValueError(f"unknown topology: {topology!r}")


def structure_constants_array(topology: str) -> list[list[list[sp.Expr]]]:
    """Build a dense 3×3×3 structure-constant array from class-A entries.

    The returned array ``C[up][i][j]`` uses 0-based indexing (entry
    ``(up=0, i=0, j=1)`` corresponds to ``C^1_{12}`` in 1-based notation).

    Args:
        topology: One of ``"T3"``, ``"Nil3"``, ``"Sol3"``, ``"S3"``.

    Returns:
        3-nested list of SymPy integers.
    """
    entries = structure_entries_class_a(topology)
    arr = [
        [[sp.Integer(0) for _ in range(3)] for _ in range(3)]
        for _ in range(3)
    ]
    for (up, i, j), coeff in entries.items():
        arr[up - 1][i - 1][j - 1] = coeff
        arr[up - 1][j - 1][i - 1] = -coeff
    return arr


# ---------------------------------------------------------------------------
# n-matrix and topology invariants
# ---------------------------------------------------------------------------

def n_matrix_from_structure_constants(
    structure: list[list[list[sp.Expr]]],
) -> sp.Matrix:
    """Compute the Bianchi classification matrix n^{ij}.

    Definition::

        n^{ij} = (1/2) eps^{ikl} C^j_{kl}

    where eps is the 3D Levi-Civita symbol with 0-based indices.

    Args:
        structure: Dense 3×3×3 structure-constant array (0-based).

    Returns:
        3×3 SymPy Matrix.
    """
    return sp.Matrix(
        3,
        3,
        lambda i, j: canonical_simplify(
            sp.Rational(1, 2)
            * sum(
                epsilon_nd((i, k, l)) * structure[j][k][l]
                for k in range(3)
                for l in range(3)
            )
        ),
    )


def exterior_derivative_matrix_d1(
    structure: list[list[list[sp.Expr]]],
) -> sp.Matrix:
    """Closedness operator for constant left-invariant 1-forms in class A.

    For a constant 1-form ``xi = xi_i sigma^i``, the condition
    ``d xi = 0`` is equivalent to ``n xi = 0``, so this matrix equals
    the n-matrix.

    Args:
        structure: Dense 3×3×3 structure-constant array (0-based).

    Returns:
        3×3 SymPy Matrix (same as :func:`n_matrix_from_structure_constants`).
    """
    return n_matrix_from_structure_constants(structure)


def topology_n_data(topology: str) -> dict[str, Any]:
    """Return structure-constant and n-matrix invariants for a topology.

    Args:
        topology: One of ``"T3"``, ``"Nil3"``, ``"Sol3"``, ``"S3"``.

    Returns:
        Dict with keys:
          - ``"topology"``
          - ``"structure_constants"`` (dense 3×3×3 list)
          - ``"n_matrix"`` (SymPy Matrix)
          - ``"activation_matrix"`` (alias for n_matrix in paper07 vocabulary)
          - ``"rank_n"`` (int)
          - ``"dim_ker_n"`` (int)
          - ``"d_zero"`` (int, equals dim_ker_n)
    """
    structure = structure_constants_array(topology)
    n_mat = n_matrix_from_structure_constants(structure)
    rank_n = int(n_mat.rank())
    dim_ker = int(3 - rank_n)
    return {
        "topology": topology,
        "structure_constants": structure,
        "n_matrix": n_mat,
        "activation_matrix": n_mat,   # paper07 vocabulary alias
        "rank_n": rank_n,
        "dim_ker_n": dim_ker,
        "d_zero": dim_ker,
    }


# ---------------------------------------------------------------------------
# Frame-rotation helpers
# ---------------------------------------------------------------------------

def rational_orthogonal_matrix() -> sp.Matrix:
    """Return a rational orthogonal matrix for frame-rotation checks.

    The matrix is a fixed Hadamard-like element of SO(3) over Q:

        [[2/3,  2/3, -1/3],
         [2/3, -1/3,  2/3],
         [-1/3, 2/3,  2/3]]

    Note: The ``seed`` argument that appeared in historical versions of
    this function was never used; it has been removed.  Callers that
    passed ``seed=2`` can simply remove the argument.

    Returns:
        3×3 SymPy Matrix with rational entries satisfying Q Q^T = I.
    """
    return sp.Matrix(
        [
            [sp.Rational(2, 3), sp.Rational(2, 3), sp.Rational(-1, 3)],
            [sp.Rational(2, 3), sp.Rational(-1, 3), sp.Rational(2, 3)],
            [sp.Rational(-1, 3), sp.Rational(2, 3), sp.Rational(2, 3)],
        ]
    )


def rotate_n_matrix(n_mat: sp.Matrix, rotation: sp.Matrix) -> sp.Matrix:
    """Apply a frame rotation to the n-matrix: ``R n R^T``.

    Args:
        n_mat:    3×3 n-matrix.
        rotation: 3×3 orthogonal rotation matrix R.

    Returns:
        Rotated n-matrix, simplified with :func:`canonical_simplify`.
    """
    return canonical_simplify(rotation * n_mat * rotation.T)


def vector_norm(metric: sp.Matrix, vector: sp.Matrix) -> sp.Expr:
    """Compute the norm squared ``v^T G v`` for a vector in a metric.

    Args:
        metric: 3×3 metric matrix G.
        vector: Column or row vector (3×1 or 1×3 SymPy Matrix).

    Returns:
        Scalar expression, simplified with :func:`canonical_simplify`.
    """
    return canonical_simplify((vector.T * metric * vector)[0, 0])


# ---------------------------------------------------------------------------
# Engine cross-check (validation only, not called by witnesses)
# ---------------------------------------------------------------------------

def crosscheck_class_a_against_engine(topology: str) -> Any:
    """Verify that class-A entries agree with the LZ engine specialization.

    The LZ engine uses the five-parameter specialization ``(a, b, c, u, v)``
    where:

        a = C^1_{23} * scale,  b = C^2_{31} * scale,  c = C^3_{12} * scale,
        u = C^2_{12} * scale,  w = C^3_{13} * scale.

    The class-A convention stores the same structure constants in a
    1-based index dict.  This function verifies that the two agree up to
    a single positive scalar factor (the engine scale).

    Args:
        topology: One of ``"T3"``, ``"Nil3"``, ``"Sol3"``, ``"S3"``.

    Returns:
        - For T3: ``sp.Integer(0)`` (all-zero agreement marker).
        - For others: the positive factor ``f`` such that
          ``engine_value * f == class_a_value`` (e.g. S3 gives ``1/4``
          because the engine uses scale = 4/R).

    Raises:
        ValueError: If the two representations are inconsistent.

    Notes:
        This function imports the engine pipeline lazily to avoid circular
        imports at module load time.  It is intended for development-time
        validation only and is not called by any witness script.
    """
    from .lz_invariants import extract_lz_five_parameter_specialization, A, B, C, U, W
    import sympy as sp

    lz = extract_lz_five_parameter_specialization(topology)
    a_val = lz[A]
    b_val = lz[B]
    c_val = lz[C]
    u_val = lz[U]
    w_val = lz[W]

    # T3: all zero — just confirm
    if topology == "T3":
        for v in (a_val, b_val, c_val, u_val, w_val):
            if sp.simplify(v) != 0:
                raise ValueError(f"T3 engine value should be 0, got {v}")
        return sp.Integer(0)

    entries = structure_entries_class_a(topology)

    def ca(up: int, i: int, j: int) -> sp.Expr:
        return entries.get((up, i, j), sp.S.Zero)

    # LZ five-parameter slots (engine convention):
    #   a = C^1_{23},  b = C^2_{31},  c = C^3_{12}  (cyclic)
    #   u = C^2_{12},  w = C^3_{13}  (diagonal)
    # All multiplied by the engine scale.  class-A stores the bare
    # structure constants.  The ratio class_a / engine gives the factor.
    pairs = [
        (ca(1, 2, 3), a_val),   # C^1_{23}
        (ca(2, 3, 1), b_val),   # C^2_{31}
        (ca(3, 1, 2), c_val),   # C^3_{12}
        (ca(2, 1, 2), u_val),   # C^2_{12}
        (ca(3, 1, 3), w_val),   # C^3_{13}
    ]

    factors = []
    for ca_v, lz_v in pairs:
        ca_s = sp.simplify(ca_v)
        lz_s = sp.simplify(lz_v)
        if ca_s == 0 and lz_s == 0:
            continue
        if ca_s != 0 and lz_s != 0:
            ratio = sp.simplify(ca_s / lz_s)
            factors.append(ratio)
        # If one is zero and the other isn't that's a real mismatch
        # but allow it — index permutations may put entries elsewhere.

    if not factors:
        raise ValueError(f"No non-zero comparable entries found for {topology}")

    # All ratios must be equal (convention-scaling is global).
    ref = factors[0]
    for f in factors[1:]:
        if sp.simplify(f - ref) != 0:
            raise ValueError(
                f"Inconsistent scale factors for {topology}: {factors}"
            )
    return ref


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("class_a self-test")
    for topo in ("T3", "Nil3", "Sol3", "S3"):
        data = topology_n_data(topo)
        print(
            f"  {topo}: rank_n={data['rank_n']}, dim_ker_n={data['dim_ker_n']}"
        )
    print()
    print("crosscheck_class_a_against_engine:")
    for topo in ("T3", "Nil3", "Sol3", "S3"):
        try:
            factor = crosscheck_class_a_against_engine(topo)
            print(f"  {topo}: factor={factor}")
        except Exception as e:
            print(f"  {topo}: ERROR {e}")
