"""
DPPU Action Layer
=================

Lagrangian, effective potential, EC action builder, and reduced-sector
helpers.  Heavy reduced-sector helpers live in ``dppu.action.reduced_sector``
and are imported explicitly by callers to avoid pipeline import cycles.

``palatini_first_order`` provides the class-A minisuperspace EH density
and Palatini reduction used by paper07 witness scripts.

``finite_reduction`` provides the background-adapted finite-amplitude
Palatini symbolic series engine used by the finite_fate_branch_witness.
"""

from .lagrangian import compute_lagrangian, compute_action
from .potential import compute_effective_potential, get_potential_function, subs_zero_modes
from .ec_action import compute_c2_ec, build_veff_ec, build_veff_ec_func
from .palatini_first_order import (
    eh_density,
    palatini_reduced_density,
    minisuperspace_dewitt_metric,
)
from .finite_reduction import palatini_symbolic_series

__all__ = [
    'compute_lagrangian',
    'compute_action',
    'compute_effective_potential',
    'get_potential_function',
    'subs_zero_modes',
    'compute_c2_ec',
    'build_veff_ec',
    'build_veff_ec_func',
    'eh_density',
    'palatini_reduced_density',
    'minisuperspace_dewitt_metric',
    'palatini_symbolic_series',
]
