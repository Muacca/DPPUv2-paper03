"""
UnifiedEngine — Configurable Engine for All Geometric Degrees of Freedom
========================================================================

Single engine class that reproduces all 12 specialized topology engines as
special cases, by selecting which geometric DOFs are active at runtime.

DOF Summary
-----------
Topology     : S³×S¹ | T³×S¹ | Nil³×S¹
Squash  (ε)  : axisymmetric squash of M³
Shear   (s)  : 0–1 plane shear (S³ only)
Fiber        : NONE | TWIST (ω₁,ω₂,ω₃) | MIXING (δ)
Torsion      : AX | VT | MX

Equivalence table
-----------------
Legacy engine                       DOFConfig equivalent
---------------------------------   -------------------------------------------------
S3S1Engine                          S3, squash=True,  shear=False, fiber=NONE
T3S1Engine                          T3, squash=False, shear=False, fiber=NONE
Nil3S1Engine                        NIL3, squash=True, shear=False, fiber=NONE
TwistS3S1Engine                     S3, squash=False, shear=False, fiber=TWIST, iso=True
TwistNil3S1Engine                   NIL3, squash=False, shear=False, fiber=TWIST, iso=True
VectorTwistS3S1Engine               S3, squash=False, shear=False, fiber=TWIST, iso=False
VectorTwistT3S1Engine               T3, squash=False, shear=False, fiber=TWIST, iso=False
SquashedVectorTwistS3S1Engine       S3, squash=True,  shear=False, fiber=TWIST, iso=False
FullyAnisotropicS3S1Engine          S3, squash=True,  shear=True,  fiber=TWIST, iso=False
SquashedOffDiagS3S1Engine           S3, squash=True,  shear=False, offdiag=True, fiber=NONE
M1MixingS3S1Engine                  S3, squash=False, shear=False, fiber=MIXING
M1SquashedMixingS3S1Engine          S3, squash=True,  shear=False, fiber=MIXING
FullyAnisotropicMixingS3S1Engine    S3, squash=True,  shear=True,  fiber=MIXING

Usage
-----
    from dppu.topology.unified import UnifiedEngine, DOFConfig, TopologyType, FiberMode
    from dppu.torsion.mode import Mode
    from dppu.torsion.nieh_yan import NyVariant

    # Build a config that replicates TwistS3S1Engine:
    cfg = DOFConfig(
        topology=TopologyType.S3,
        fiber_mode=FiberMode.TWIST,
        isotropic_twist=True,
        torsion_mode=Mode.AX,
    )
    engine = UnifiedEngine(cfg)
    engine.run()

    # Or use a preset:
    cfg = DOFConfig.from_engine('FullyAnisotropicMixingS3S1Engine')
    cfg.torsion_mode = Mode.AX
    engine = UnifiedEngine(cfg)
    engine.run()

Author: Claude (Phase 4 Integration, 2026-03-04)
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, Optional

from sympy import symbols, Matrix, S, pi, Rational, sqrt
from sympy.tensor.array import MutableDenseNDimArray

from ..engine.pipeline import BaseFrameEngine
from ..torsion.mode import Mode
from ..torsion.nieh_yan import NyVariant
from ..utils.levi_civita import epsilon_symbol


# ── Enumerations ──────────────────────────────────────────────────────────────

class TopologyType(Enum):
    """Base topology of the M³×S¹ manifold."""
    S3   = "s3"    # S³×S¹  — SU(2) Lie group × circle
    T3   = "t3"    # T³×S¹  — flat 3-torus × circle
    NIL3 = "nil3"  # Nil³×S¹ — Heisenberg nilmanifold × circle


class FiberMode(Enum):
    """Deformation type applied to the S¹ fiber direction."""
    NONE    = "none"    # ê³ = L dτ  (unmodified)
    TWIST   = "twist"   # ê³ = L(dτ + ω₁σ⁰ + ω₂σ¹ + ω₃σ²)
    MIXING  = "mixing"  # ê², ê³ rotated by mixing angle δ  (S³ only)
    BOTH    = "both"    # TWIST + MIXING simultaneously (S³ only)
                        # Structure constants computed via full frame rotation
                        # (includes O(δω) cross terms automatically)


# ── DOF Configuration ─────────────────────────────────────────────────────────

@dataclass
class DOFConfig:
    """
    Configuration of active geometric degrees of freedom.

    Parameters
    ----------
    topology : TopologyType
        Base manifold topology (S3, T3, NIL3).
    enable_squash : bool
        Activate axisymmetric squash parameter ε.
        S3: breaks SO(3) → U(1).  Nil3: scales C²₀₁ by (1+ε)^{−4/3}.
    enable_shear : bool
        Activate 0–1 plane shear parameter s (S3 only).
        Further breaks U(1) symmetry of squash (f₀ ≠ f₁ ≠ f₂).
    enable_offdiag_shear : bool
        Activate off-diagonal shear parameters q₃, q₄, q₅ (S3 only).
        Corresponds to spin-2 generators T₃ (01-plane), T₄ (02-plane), T₅ (12-plane).
        Frame: F = F_diag(ε,s) × G(q₃,q₄,q₅).  At q₃=q₄=q₅=0: G=I (backward compat).

        WARNING — Out-of-Memory risk: Running the full pipeline (engine.run())
        with enable_offdiag_shear=True causes OOM even on machines with 64 GB RAM.
        Steps E4.3a and E4.7 (Riemann tensor) produce 5-variable rational products
        of ~3M SymPy ops that exhaust heap memory before completion.
        For off-diagonal Hessian computations, use the single-variable Method 2
        approach (z = exp(q/√2) rational substitution) or a numerical pipeline.
    fiber_mode : FiberMode
        Type of S¹ fiber deformation (NONE | TWIST | MIXING).
    isotropic_twist : bool
        Used only when fiber_mode=TWIST.
        True  → single symbol ω shared across all three twist components.
        False → independent symbols ω₁, ω₂, ω₃.
    torsion_mode : Mode
        AX (η≠0, V=0) | VT (η=0, V≠0) | MX (both).
    ny_variant : NyVariant
        Nieh-Yan term variant: TT | REE | FULL.

    Class Methods
    -------------
    from_engine(name) → DOFConfig
        Return a preset that reproduces the named legacy engine geometry
        (torsion_mode / ny_variant are left at defaults and should be set
        by the caller).
    """
    topology:           TopologyType = TopologyType.S3
    enable_squash:      bool         = False
    enable_shear:       bool         = False
    enable_offdiag_shear: bool       = False
    fiber_mode:         FiberMode    = FiberMode.NONE
    isotropic_twist:    bool         = True
    torsion_mode:       Mode         = Mode.AX
    ny_variant:         NyVariant    = NyVariant.FULL
    enable_velocity:    bool         = False  # Velocity mode for G metric computation (BOTH mode only)
    skip_antisymmetry_check: bool    = False  # Skip Riemann antisymmetry verification (speedup Option B)
                                              # Enable only when correctness is guaranteed by construction,
                                              # e.g., off-diagonal computations using rationalized variables (z_i).

    # ------------------------------------------------------------------
    @classmethod
    def from_engine(cls, engine_name: str) -> "DOFConfig":
        """
        Return the DOFConfig preset that matches a legacy specialized engine.

        Geometry DOFs are set to match the legacy engine exactly.
        The caller must set `torsion_mode` and `ny_variant` separately.

        Parameters
        ----------
        engine_name : str
            Name of the legacy engine class.  Accepted values:
            'S3S1Engine', 'T3S1Engine', 'Nil3S1Engine',
            'TwistS3S1Engine', 'TwistNil3S1Engine',
            'VectorTwistS3S1Engine', 'VectorTwistT3S1Engine',
            'SquashedVectorTwistS3S1Engine',
            'FullyAnisotropicS3S1Engine',
            'M1MixingS3S1Engine', 'M1SquashedMixingS3S1Engine',
            'FullyAnisotropicMixingS3S1Engine'

        Returns
        -------
        DOFConfig
        """
        _PRESETS: Dict[str, "DOFConfig"] = {
            'S3S1Engine': cls(
                topology=TopologyType.S3,
                enable_squash=True, enable_shear=False,
                fiber_mode=FiberMode.NONE),
            'T3S1Engine': cls(
                topology=TopologyType.T3,
                enable_squash=False, enable_shear=False,
                fiber_mode=FiberMode.NONE),
            'Nil3S1Engine': cls(
                topology=TopologyType.NIL3,
                enable_squash=True, enable_shear=False,
                fiber_mode=FiberMode.NONE),
            'TwistS3S1Engine': cls(
                topology=TopologyType.S3,
                enable_squash=False, enable_shear=False,
                fiber_mode=FiberMode.TWIST, isotropic_twist=True),
            'TwistNil3S1Engine': cls(
                topology=TopologyType.NIL3,
                enable_squash=False, enable_shear=False,
                fiber_mode=FiberMode.TWIST, isotropic_twist=True),
            'VectorTwistS3S1Engine': cls(
                topology=TopologyType.S3,
                enable_squash=False, enable_shear=False,
                fiber_mode=FiberMode.TWIST, isotropic_twist=False),
            'VectorTwistT3S1Engine': cls(
                topology=TopologyType.T3,
                enable_squash=False, enable_shear=False,
                fiber_mode=FiberMode.TWIST, isotropic_twist=False),
            'SquashedVectorTwistS3S1Engine': cls(
                topology=TopologyType.S3,
                enable_squash=True, enable_shear=False,
                fiber_mode=FiberMode.TWIST, isotropic_twist=False),
            'FullyAnisotropicS3S1Engine': cls(
                topology=TopologyType.S3,
                enable_squash=True, enable_shear=True,
                fiber_mode=FiberMode.TWIST, isotropic_twist=False),
            # WARNING: engine.run() with this config causes Out-of-Memory even on
            # 64 GB machines.  Structure constants (E4.2) and the LC connection
            # (E4.3) complete successfully, but the Riemann tensor steps (E4.3a,
            # E4.7) generate ~3M-op intermediate polynomials and exhaust heap.
            # Use the single-variable Method 2 (z = exp(q/√2)) or a numerical
            # finite-difference pipeline for off-diagonal Hessian calculations.
            'SquashedOffDiagS3S1Engine': cls(
                topology=TopologyType.S3,
                enable_squash=True, enable_shear=False,
                enable_offdiag_shear=True,
                fiber_mode=FiberMode.NONE),
            'M1MixingS3S1Engine': cls(
                topology=TopologyType.S3,
                enable_squash=False, enable_shear=False,
                fiber_mode=FiberMode.MIXING),
            'M1SquashedMixingS3S1Engine': cls(
                topology=TopologyType.S3,
                enable_squash=True, enable_shear=False,
                fiber_mode=FiberMode.MIXING),
            'FullyAnisotropicMixingS3S1Engine': cls(
                topology=TopologyType.S3,
                enable_squash=True, enable_shear=True,
                fiber_mode=FiberMode.MIXING),
            'FullyAnisotropicBothS3S1Engine': cls(
                topology=TopologyType.S3,
                enable_squash=True, enable_shear=True,
                fiber_mode=FiberMode.BOTH, isotropic_twist=False),
        }
        if engine_name not in _PRESETS:
            raise ValueError(
                f"Unknown engine name: {engine_name!r}.\n"
                f"Available: {sorted(_PRESETS.keys())}"
            )
        return _PRESETS[engine_name]


# ── Unified Engine — backward-compatible factory ──────────────────────────────

class UnifiedEngine:
    """
    Backward-compatible factory.

    Accepts a DOFConfig and returns an instance of the appropriate
    topology implementation class (S3Geometry / T3Geometry / Nil3Geometry).

    Interface is identical to the old monolithic class::

        engine = UnifiedEngine(cfg, logger)
        engine.run()
        engine.data['potential']     # works as before
        engine.get_riemann_lambdified()
        engine.get_free_params()
        engine.get_config_summary()

    Parameters
    ----------
    config : DOFConfig
        Geometry and torsion configuration (topology / DOFs / mode).
    logger : optional
        Logger instance; NullLogger is used when omitted.
    checkpoint_mgr : optional
        Checkpoint manager for step-level checkpointing.
    """

    def __new__(
        cls,
        config: "DOFConfig",
        logger=None,
        checkpoint_mgr=None,
    ):
        from .s3   import S3Geometry
        from .t3   import T3Geometry
        from .nil3 import Nil3Geometry

        if not isinstance(config, DOFConfig):
            raise TypeError(
                f"config must be DOFConfig, got {type(config)}"
            )

        _MAP = {
            TopologyType.S3:   S3Geometry,
            TopologyType.T3:   T3Geometry,
            TopologyType.NIL3: Nil3Geometry,
        }
        impl_cls = _MAP.get(config.topology)
        if impl_cls is None:
            raise ValueError(
                f"Unknown topology: {config.topology!r}"
            )

        # Delegate to the concrete topology class.
        # Its __init__ (TopologyEngine.__init__) is called normally.
        return impl_cls(config, logger, checkpoint_mgr)

