"""
DPPU Dirac Constraint Closure Layer
=====================================

Dirac constraint-closure engine for Branch-A first-order Palatini and
Routh-before-Dirac route variant audit.

``constraint_closure`` provides the core DBR engine: Branch-A local L2
builder, harmonic projection, and first/second-class split via Poisson rank.

``route_variants`` provides RBD route variant audit helpers: Routh
elimination (spatial-omega or full-K) followed by scalar projection and
second-order Dirac closure.
"""

from .constraint_closure import (
    build_branchA_local_L2,
    project_lagrangian,
    build_dirac_for_projection,
    build_dirac_for_second_order_lagrangian,
    sample_geometry,
    harmonic_lambda,
    poisson,
    coeff_vector,
    is_independent,
    eps3,
    lc4,
)
from .route_variants import (
    routh_eliminate,
    project_reduced_dirac,
)

__all__ = [
    'build_branchA_local_L2',
    'project_lagrangian',
    'build_dirac_for_projection',
    'build_dirac_for_second_order_lagrangian',
    'sample_geometry',
    'harmonic_lambda',
    'poisson',
    'coeff_vector',
    'is_independent',
    'eps3',
    'lc4',
    'routh_eliminate',
    'project_reduced_dirac',
]
