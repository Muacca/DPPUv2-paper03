# Engine Layer

⇒ [日本語](README_ja.md)

Module group providing computation pipeline foundations.

## Overview

Provides the 15-step computation pipeline, logging, and checkpoint functionality.

## Modules

### pipeline.py

Base class for the computation pipeline.

**BaseFrameEngine:**

Abstract base class for all topology engines. Executes 15 computation steps sequentially:

1. **E4.1**: Parameter setup
2. **E4.2**: Metric and frame definition
3. **E4.3**: Christoffel symbols (general Koszul formula)
4. **E4.4**: Torsion tensor construction
5. **E4.5**: Contortion tensor
6. **E4.6**: Einstein-Cartan connection
7. **E4.7**: Riemann tensor (with strict verification)
8. **E4.8**: Ricci scalar
9. **E4.9**: Torsion scalar
10. **E4.10**: Nieh-Yan density (all variants)
11. **E4.11**: Lagrangian construction
12. **E4.12**: Angular integration
13. **E4.13**: Effective potential extraction
14. **E4.14**: Stability analysis
15. **E4.15**: Summary

**Abstract methods (implemented by subclasses):**

```python
@abstractmethod
def _define_structure_constants(self):
    """Define structure constants C^a_{bc}"""

@abstractmethod
def _define_metric(self):
    """Define frame metric η_{ab}"""

@abstractmethod
def _get_volume_factor(self):
    """Return volume factor"""
```

### logger.py

Computation progress logging.

**ComputationLogger:**

Log output to file and console:

```python
from dppu.engine import ComputationLogger

logger = ComputationLogger('dppu_run.log')
logger.step(1, "E4.1", "Parameter setup")
logger.info("Setting r = 1.0")
logger.success("Step completed")
```

**NullLogger:**

Disable log output (for testing):

```python
from dppu.engine import NullLogger

logger = NullLogger()  # No output
```

### checkpoint.py

Save and restore computation state.

**CheckpointManager:**

Save state at each step via pickle:

```python
from dppu.engine import CheckpointManager

ckpt = CheckpointManager('checkpoints/', enabled=True)

# Save
ckpt.save('E4.7', {'R': riemann_tensor})

# Restore
data = ckpt.load('E4.7')
```

**Enabling checkpoints:**

```bash
python scripts/pipeline/run_s3s1.py --mode MX --ny-variant FULL \
    --checkpoint-dir ./checkpoints
```

## Usage

```python
from dppu.topology import S3S1Engine
from dppu.torsion import Mode, NyVariant
from dppu.engine import ComputationLogger, CheckpointManager

# Set up logger and checkpoint
logger = ComputationLogger('run.log')
ckpt = CheckpointManager('ckpt/', enabled=True)

# Run engine
engine = S3S1Engine(Mode.MX, NyVariant.FULL, logger, ckpt)
engine.run()
```

## Pipeline Customization

Execute specific steps only:

```python
engine = S3S1Engine(Mode.MX, NyVariant.FULL)

# Execute steps 1-7 only (up to curvature computation)
engine.run(stop_after='E4.7')

# Resume from step 7
engine.run(resume_from='E4.7')
```

## Dependencies

- SymPy (symbolic computation)
- pickle (checkpoints)

## Related Modules

- [topology](../topology/README.md): Topology-specific engines
