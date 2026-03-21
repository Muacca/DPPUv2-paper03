"""
Nil3Geometry — Nil³×S¹ Topology Implementation
================================================

EC+NY engine on Nil³ (Heisenberg nilmanifold) × S¹.

Characteristics of Nil³:
    - Only C²₀₁ is a nonzero structure constant  (1 axis vs. 6 axes for S³)
    - The radial coordinate follows the convention of using uppercase 'R'
      (an alias 'r' is also set for pipeline compatibility)
    - ε (squash) enters C²₀₁ as (1+ε)^{−4/3}
    - s (shear) does not exist
    - Under TWIST, since only dσ² = −σ⁰∧σ¹ is nonzero,
      **only ω₃** contributes to the structure constant C³₀₁ (ω₁, ω₂ are inactive)

Structure constants:
    Spatial block: C²₀₁ = −(1+ε)^{−4/3} / R  (sign: from the exterior derivative of σ²)
    TWIST:         C³₀₁ = −Lω₃ / R²

Volume factor: (2π)⁴ · L · R³
"""

from sympy import pi, S, Rational
from sympy.tensor.array import MutableDenseNDimArray
from typing import Any, Dict

from .base_topology import TopologyEngine


class Nil3Geometry(TopologyEngine):
    """
    Engine dedicated to the Nil³×S¹ topology.

    Simpler than S³ but uses a different symbol naming convention ('R' vs 'r').

    Examples
    --------
    Direct usage::

        from dppu.topology.nil3 import Nil3Geometry
        from dppu.topology.unified import DOFConfig
        from dppu.torsion.mode import Mode

        cfg = DOFConfig.from_engine('TwistNil3S1Engine')
        cfg.torsion_mode = Mode.AX
        engine = Nil3Geometry(cfg)
        engine.run()
    """

    # ── E4.1 hook ─────────────────────────────────────────────────────────────

    def _build_radial_and_deformation_params(self, params: Dict) -> None:
        """
        Add Nil³-specific symbols to params.

            R       : radial coordinate (positive — Nil³ convention)
            r       : alias for R (pipeline compatibility)
            epsilon : squash parameter (Symbol only when enable_squash=True)
            s       : always Zero  (no shear in Nil³)
        """
        from sympy import symbols
        cfg = self.config

        R = symbols('R', positive=True, real=True)
        params['R'] = R
        params['r'] = R   # pipeline-compatibility alias for E4.4 (torsion ansatz), etc.

        params['epsilon'] = (
            symbols('epsilon', real=True) if cfg.enable_squash else S.Zero
        )
        params['s'] = S.Zero   # Nil³: shear not allowed

    # ── E4.2 hook ─────────────────────────────────────────────────────────────

    def _build_structure_constants(
        self, params: Dict, C: MutableDenseNDimArray
    ) -> None:
        """
        Write the Nil³ structure constants into C in-place.

        Spatial block (1 axis only):
            C²₀₁ = −(1+ε)^{−4/3} / R
            C²₁₀ = +(1+ε)^{−4/3} / R  (antisymmetric)

        Additional TWIST terms (only ω₃ is active):
            C³₀₁ = −L ω₃ / R²
            C³₁₀ = +L ω₃ / R²
        """
        from .unified import FiberMode

        R       = params['R']
        L       = params['L']
        epsilon = params['epsilon']

        # Spatial block: C²₀₁
        nil3_factor = (1 + epsilon) ** Rational(-4, 3)
        C[2, 0, 1] = -nil3_factor / R
        C[2, 1, 0] = +nil3_factor / R

        # TWIST: In Nil³, only dσ² = −σ⁰∧σ¹ is nonzero → only ω₃ contributes
        if self.config.fiber_mode == FiberMode.TWIST:
            omega3 = params['omega3']  # iso=True → omega3=omega; iso=False → independent symbol
            twist_factor = L * omega3 / R ** 2
            C[3, 0, 1] = -twist_factor
            C[3, 1, 0] = +twist_factor

    def _compute_volume(self, params: Dict) -> Any:
        """(2π)⁴ · L · R³"""
        return (2 * pi) ** 4 * params['L'] * params['R'] ** 3
