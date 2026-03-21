"""
DPPU Engine Layer
=================

Provides computation pipeline infrastructure.

Modules:
- pipeline: BaseFrameEngine and step orchestration
- logger: ComputationLogger and NullLogger
- checkpoint: CheckpointManager for save/restore
"""

from .logger import ComputationLogger, NullLogger
from .checkpoint import CheckpointManager

__all__ = [
    'ComputationLogger',
    'NullLogger',
    'CheckpointManager',
]
