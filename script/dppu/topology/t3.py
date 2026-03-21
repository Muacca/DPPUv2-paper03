"""
T3Geometry — T³×S¹ Topology Implementation
============================================

EC+NY engine on T³ (3-dimensional torus, flat manifold) × S¹.

Characteristics:
    - All structure constants are zero  (C^a_{bc} = 0  for all a,b,c)
    - ε (squash) and s (shear) do not exist  (T³ is not a Lie group)
    - For TWIST, fiber symbols (ω₁,ω₂,ω₃) enter params but have no physical effect
      because dσ^i = 0, so they do not contribute to the structure constants
    - MIXING is specific to S³ and should be restricted at the DOFConfig level

Volume factor: (2π)⁴ · L · r³
"""

from sympy import pi, S
from sympy.tensor.array import MutableDenseNDimArray
from typing import Any, Dict

from .base_topology import TopologyEngine


class T3Geometry(TopologyEngine):
    """
    Engine dedicated to the T³×S¹ topology.

    The simplest implementation: all structure constants are zero.
    ``_build_structure_constants`` does nothing (C remains zero).

    Examples
    --------
    Direct usage::

        from dppu.topology.t3 import T3Geometry
        from dppu.topology.unified import DOFConfig
        from dppu.torsion.mode import Mode

        cfg = DOFConfig.from_engine('VectorTwistT3S1Engine')
        cfg.torsion_mode = Mode.MX
        engine = T3Geometry(cfg)
        engine.run()
    """

    # ── E4.1 hook ─────────────────────────────────────────────────────────────

    def _build_radial_and_deformation_params(self, params: Dict) -> None:
        """
        Add T³-specific symbols to params.

            r       : radial coordinate (positive)
            epsilon : always Zero  (no squash parameter for T³)
            s       : always Zero  (no shear parameter for T³)
        """
        from sympy import symbols
        params['r']       = symbols('r', positive=True, real=True)
        params['epsilon'] = S.Zero   # T³: squash not allowed
        params['s']       = S.Zero   # T³: shear not allowed

    # ── E4.2 hook ─────────────────────────────────────────────────────────────

    def _build_structure_constants(
        self, params: Dict, C: MutableDenseNDimArray
    ) -> None:
        """
        T³ is a flat manifold, so all structure constants are zero.
        Does nothing (C is already zero-initialized).
        """
        # T³: dσ^i = 0  →  C^a_{bc} = 0  for all indices
        pass

    def _compute_volume(self, params: Dict) -> Any:
        """(2π)⁴ · L · r³"""
        return (2 * pi) ** 4 * params['L'] * params['r'] ** 3
