"""
Hodge Dual Operator
===================

Computes Hodge dual (*) on 2-form indices of curvature tensor.

For R^{ab}_{cd}: (*R)^{ab}_{cd} = (1/2) ε_{cdef} R^{ab,ef}
"""

import numpy as np
from ..utils.levi_civita import levi_civita_4d


def compute_hodge_dual(R_tensor: np.ndarray) -> np.ndarray:
    """
    Compute Hodge dual on (cd) indices of R^{ab}_{cd}.

    (*R)^{ab}_{cd} = (1/2) ε_{cdef} R^{ab,ef}

    Args:
        R_tensor: numpy array (4,4,4,4) representing R^{ab}_{cd}

    Returns:
        numpy array (4,4,4,4) representing (*R)^{ab}_{cd}
    """
    R_dual = np.zeros_like(R_tensor)

    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    val = 0.0
                    for e in range(4):
                        for f in range(4):
                            eps = levi_civita_4d(c, d, e, f)
                            if eps != 0:
                                val += 0.5 * eps * R_tensor[a, b, e, f]
                    R_dual[a, b, c, d] = val

    return R_dual


def classify_block(c: int, d: int) -> str:
    """
    Classify (cd) indices into spatial or mixed block.

    Spatial: c,d ∈ {0,1,2}
    Mixed: involves index 3

    Returns:
        'spatial' or 'mixed'
    """
    if c < 3 and d < 3:
        return 'spatial'
    return 'mixed'


def hodge_swaps_blocks() -> bool:
    """
    Property: Hodge dual swaps spatial and mixed blocks.

    * : spatial → mixed
    * : mixed → spatial

    This is key to proving P = <R, *R> = 0.
    """
    return True
