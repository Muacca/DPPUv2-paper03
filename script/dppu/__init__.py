"""
DPPU: Differential Geometry Package for Physics Universe
=========================================================

A modular Python package for Einstein-Cartan gravity with Nieh-Yan term,
designed for publication-quality computations.

Package Structure:
- geometry/    : Metric, volume forms, structure constants
- connection/  : Levi-Civita, contortion, EC connection
- curvature/   : Riemann, Ricci, Hodge dual, self-duality
- torsion/     : Mode enum, ansatz, Nieh-Yan density
- action/      : Lagrangian, effective potential, stability
- topology/    : S3xS1, T3xS1, Nil3xS1 implementations
- engine/      : Computation pipeline, logging, checkpoints
- scanning/    : Parameter scans, Phase 1 integration
- utils/       : Shared utilities

Author: Muacca
Version: 2.0
Date: 2026-02
"""

__version__ = "2.0.0"
__author__ = "Muacca"

# Lazy imports for main components
def __getattr__(name):
    """Lazy loading of submodules."""
    if name == "Mode":
        from .torsion.mode import Mode
        return Mode
    elif name == "NyVariant":
        from .torsion.nieh_yan import NyVariant
        return NyVariant
    elif name == "BaseFrameEngine":
        from .engine.pipeline import BaseFrameEngine
        return BaseFrameEngine
    elif name == "S3S1Engine":
        from .topology.s3s1 import S3S1Engine
        return S3S1Engine
    elif name == "T3S1Engine":
        from .topology.t3s1 import T3S1Engine
        return T3S1Engine
    elif name == "Nil3S1Engine":
        from .topology.nil3s1 import Nil3S1Engine
        return Nil3S1Engine
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
