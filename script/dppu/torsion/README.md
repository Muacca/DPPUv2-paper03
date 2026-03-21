# Torsion Layer

⇒ [日本語](README_ja.md)

Module group defining torsion structures.

## Overview

Provides torsion ansatz and Nieh-Yan density computation.

## Modules

### mode.py

Torsion mode definitions.

**Mode Enum:**

| Mode | Description | Parameters |
|------|-------------|-----------|
| `AX` | Axial torsion only | η |
| `VT` | Vector trace only | V |
| `MX` | Mixed mode | η, V |

```python
from dppu.torsion import Mode

mode = Mode.MX  # Mixed mode
```

### nieh_yan.py

Nieh-Yan variant definitions and NY density computation.

**NyVariant Enum:**

| Variant | Definition | Physical Meaning |
|---------|-----------|-----------------|
| `TT` | T^a ∧ T_a | Torsion-torsion term |
| `REE` | e^a ∧ e^b ∧ R_ab | Curvature-derived term |
| `FULL` | TT - REE | Full Nieh-Yan density |

**Key functions:**

- `compute_nieh_yan_TT(T)`: Compute TT term
- `compute_nieh_yan_REE(R)`: Compute REE term
- `compute_nieh_yan_full(T, R)`: Full NY density

### ansatz.py

Torsion tensor ansatz construction.

**T1: Totally antisymmetric (axial)**

```
T^{ijk} = (2η/r) ε^{ijk}  (i,j,k ∈ {0,1,2})
```

Exists only in spatial directions, defining the axial vector S^μ:
```
S^μ = (1/6) ε^{μνρσ} T_{νρσ}
```

**T2: Vector trace**

```
T^{abc} = (1/3)(η^{ac}T^b - η^{ab}T^c)
T^μ = (0, 0, 0, V)
```

Exists only in the S¹ direction.

**Key functions:**

- `build_torsion_tensor(mode, eta, V, r)`: Construct T^{abc} from ansatz

### scalar.py

Torsion scalar computation.

**Torsion scalar T:**

```
T = T_{abc} T^{abc}
```

## Usage

```python
from dppu.torsion import Mode, NyVariant, build_torsion_tensor

# Construct torsion tensor
T = build_torsion_tensor(Mode.MX, eta=-1.0, V=2.0, r=1.0)

# Nieh-Yan density
from dppu.torsion import compute_nieh_yan_full
N = compute_nieh_yan_full(T, R)
```

## Physical Meaning

- **Axial torsion (AX)**: Related to spin-spin interactions
- **Vector torsion (VT)**: Related to mass generation mechanisms
- **Mixed mode (MX)**: General setting including both effects

## Dependencies

- [utils](../utils/README.md): Levi-Civita symbol

## Related Modules

- [connection](../connection/README.md): Compute contortion from torsion
- [curvature](../curvature/README.md): REE part of Nieh-Yan term
