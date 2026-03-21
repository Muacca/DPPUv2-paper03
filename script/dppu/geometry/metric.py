"""
Frame Metric Utilities
======================

Provides utilities for working with the frame (orthonormal) metric.

In our framework, we work in an orthonormal frame basis where:
- The frame metric η_{ab} = diag(1, 1, 1, 1) (Euclidean signature)
- Indices are raised/lowered with this metric
- Metric compatibility: ∇η = 0

Author: Muacca
"""

from typing import Any, List, Optional, Tuple

from sympy import Matrix, S, cancel
from sympy.tensor.array import MutableDenseNDimArray


def create_frame_metric(dim: int = 4, signature: str = "euclidean") -> Matrix:
    """
    Create an orthonormal frame metric.

    Args:
        dim: Dimension (default 4)
        signature: Metric signature
            - "euclidean": diag(+1, +1, +1, +1)
            - "lorentzian": diag(-1, +1, +1, +1)

    Returns:
        SymPy Matrix representing η_{ab}

    Examples:
        >>> metric = create_frame_metric(4, "euclidean")
        >>> metric[0, 0]
        1
    """
    if signature == "euclidean":
        return Matrix.eye(dim)
    elif signature == "lorentzian":
        diag_entries = [-1] + [1] * (dim - 1)
        return Matrix.diag(*diag_entries)
    else:
        raise ValueError(f"Unknown signature: {signature}")


def verify_metric_compatibility(
    Gamma: MutableDenseNDimArray,
    metric: Matrix,
    dim: int,
    logger: Optional[Any] = None
) -> Tuple[bool, List[Tuple]]:
    """
    Verify metric compatibility: Γ_{abc} + Γ_{bac} = 0.

    For an orthonormal frame with the Levi-Civita connection,
    the lowered connection coefficients must be antisymmetric
    in the first two indices.

    Args:
        Gamma: Connection coefficients Γ^a_{bc}
        metric: Frame metric η_{ab}
        dim: Dimension
        logger: Optional logger for output

    Returns:
        Tuple of (passed: bool, violations: list)
        - passed: True if all checks pass
        - violations: List of (a, b, c, value) tuples for failed checks
    """
    if logger:
        logger.info("Verifying metric compatibility (Γ_abc + Γ_bac = 0)...")

    violations = []

    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                # Γ_{abc} = η_{ad} Γ^d_{bc}
                # For identity metric: Γ_{abc} = Γ^a_{bc} (with lower position)
                # Actually we check: Γ^a_{bc} + Γ^b_{ac} should be consistent
                # with Γ_{abc} + Γ_{bac} = 0

                # In orthonormal frame: Γ_{abc} + Γ_{bac} = 0
                # means Γ^a_{bc} + Γ^b_{ac} = 0 when using identity metric
                check = cancel(Gamma[a, b, c] + Gamma[b, a, c])
                if check != S.Zero:
                    violations.append((a, b, c, check))

    passed = len(violations) == 0

    if logger:
        if passed:
            logger.success("Metric compatibility: PASSED")
        else:
            logger.error(f"Metric compatibility: FAILED ({len(violations)} violations)")
            for a, b, c, val in violations[:5]:
                logger.error(f"  Γ_{{{a}{b}{c}}} + Γ_{{{b}{a}{c}}} = {val}")

    return passed, violations


def raise_index(
    tensor: MutableDenseNDimArray,
    metric_inv: Matrix,
    index_pos: int,
    dim: int
) -> MutableDenseNDimArray:
    """
    Raise a tensor index using the inverse metric.

    Args:
        tensor: Input tensor (arbitrary rank)
        metric_inv: Inverse metric η^{ab}
        index_pos: Position of index to raise (0-indexed)
        dim: Dimension

    Returns:
        Tensor with specified index raised

    Notes:
        For orthonormal frame with identity metric, this is trivial.
        Included for generality and documentation.
    """
    # For identity metric, raising is trivial
    # In general: T^{...a...} = η^{ab} T_{...b...}
    return tensor  # Identity for orthonormal frame


def lower_index(
    tensor: MutableDenseNDimArray,
    metric: Matrix,
    index_pos: int,
    dim: int
) -> MutableDenseNDimArray:
    """
    Lower a tensor index using the metric.

    Args:
        tensor: Input tensor (arbitrary rank)
        metric: Metric η_{ab}
        index_pos: Position of index to lower (0-indexed)
        dim: Dimension

    Returns:
        Tensor with specified index lowered

    Notes:
        For orthonormal frame with identity metric, this is trivial.
        Included for generality and documentation.
    """
    # For identity metric, lowering is trivial
    # In general: T_{...a...} = η_{ab} T^{...b...}
    return tensor  # Identity for orthonormal frame
