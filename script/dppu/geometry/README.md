# Geometry Layer

⇒ [日本語](README_ja.md)

Module group providing geometric foundations.

## Overview

Defines fundamental structures of differential geometry: metric tensor, volume forms, etc.

## Modules

### metric.py

Definition and manipulation of metric tensors.

**Key classes/functions:**

- `create_frame_metric(topology)`: Generate frame metric for a given topology
- `get_volume_factor(topology)`: Compute volume factor

**Metric form:**

```
η_ab = diag(-1, +1, +1, +1)  (Lorentzian signature)
```

After compactification, this effectively becomes Euclidean signature.

## Topology-specific Metrics

| Topology | Metric Structure | Volume |
|----------|-----------------|--------|
| S³×S¹ | bi-invariant on SU(2) | 2π²Lr³ |
| T³×S¹ | flat | (2π)⁴LR₁R₂R₃ |
| Nil³×S¹ | left-invariant on Heisenberg | (2π)⁴LR³ |

## Usage

```python
from dppu.geometry import create_frame_metric

# Frame metric for S³
metric = create_frame_metric('S3')
```

## Dependencies

- SymPy (symbolic computation)
- NumPy (numerical evaluation)

## Related Modules

- [connection](../connection/README.md): Compute connections from metric
- [topology](../topology/README.md): Topology-specific metric definitions
