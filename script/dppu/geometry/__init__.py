"""
DPPU Geometry Layer
===================

Provides fundamental geometric structures.

Modules:
- metric: Frame metric and related utilities
- volume: Volume forms and integration measures
- structure_constants: Lie algebra structure constants
"""

from .metric import create_frame_metric, verify_metric_compatibility

__all__ = [
    'create_frame_metric',
    'verify_metric_compatibility',
]
