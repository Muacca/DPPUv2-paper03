"""
DPPU Topology Layer
===================

Topology-specific implementations for M³×S¹ manifolds.

Unified Engine (recommended)
-----------------------------
- UnifiedEngine : Single configurable engine covering all 12 legacy cases.
- DOFConfig     : Dataclass specifying which geometric DOFs are active.
- TopologyType  : Enum — S3 | T3 | NIL3
- FiberMode     : Enum — NONE | TWIST | MIXING

  Quick start::

      from dppu.topology import UnifiedEngine, DOFConfig, TopologyType, FiberMode
      cfg = DOFConfig.from_engine('TwistS3S1Engine')
      cfg.torsion_mode = Mode.AX
      engine = UnifiedEngine(cfg)
      engine.run()

"""

# ── Unified engine & topology implementations ─────────────────────────────────
from .unified       import UnifiedEngine, DOFConfig, TopologyType, FiberMode
from .base_topology import TopologyEngine
from .s3            import S3Geometry
from .t3            import T3Geometry
from .nil3          import Nil3Geometry

# ── Legacy engines (preserved for backward compatibility) ─────────────────────
__all__ = [
    # Unified factory & topology implementations
    'UnifiedEngine', 'DOFConfig', 'TopologyType', 'FiberMode',
    'TopologyEngine', 'S3Geometry', 'T3Geometry', 'Nil3Geometry',
    # Legacy
    'S3S1Engine', 'T3S1Engine', 'Nil3S1Engine',
    'TwistS3S1Engine', 'TwistNil3S1Engine',
    'VectorTwistS3S1Engine', 'VectorTwistT3S1Engine',
    'SquashedVectorTwistS3S1Engine',
    'FullyAnisotropicS3S1Engine',
    'M1MixingS3S1Engine',
    'M1SquashedMixingS3S1Engine',
    'FullyAnisotropicMixingS3S1Engine',
]
