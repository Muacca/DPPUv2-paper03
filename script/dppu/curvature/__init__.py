"""
DPPU Curvature Layer
====================

Provides curvature-related components for Einstein-Cartan gravity.

Modules:
- riemann: Riemann tensor computation and antisymmetry verification
- ricci: Ricci tensor and scalar
- hodge: Hodge dual operator
- self_duality: Self-duality diagnostics
- pontryagin: Pontryagin inner product diagnostics
"""

from .riemann import RiemannAntisymmetryError, verify_antisymmetry_strict
from .ricci import compute_ricci_scalar, compute_ricci_tensor
from .hodge import compute_hodge_dual
from .self_duality import SDExtensionMixin
from .pontryagin import CurvatureSDDiagnostics

__all__ = [
    'RiemannAntisymmetryError',
    'verify_antisymmetry_strict',
    'compute_ricci_scalar',
    'compute_ricci_tensor',
    'compute_hodge_dual',
    'SDExtensionMixin',
    'CurvatureSDDiagnostics',
]
