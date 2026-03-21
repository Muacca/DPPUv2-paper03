# Action Layer

⇒ [日本語](README_ja.md)

Module group for action and stability analysis.

## Overview

Provides Lagrangian construction, effective potential, and stability classification.

## Modules

### lagrangian.py

Einstein-Cartan + Nieh-Yan Lagrangian construction.

**Action:**

```
S = ∫ d⁴x √|g| L
L = R/(2κ²) + θ_NY × N
```

**Key functions:**

- `build_lagrangian(R, N, kappa, theta_NY)`: Construct Lagrangian
- `integrate_angular(L, topology)`: Integrate over angular directions

### potential.py

Effective potential computation.

**Effective potential:**

After angular integration, extract the r-dependent effective potential V(r):

```
S = ∫ dr × V_eff(r)
```

**Key functions:**

- `extract_effective_potential(S_integrated)`: Extract effective potential
- `find_extrema(V, r_range)`: Search for extrema
- `compute_barrier_height(V, r_min)`: Compute barrier height

### stability.py

Stability analysis and classification.

**Stability types:**

| Type | Condition | Physical Interpretation |
|------|-----------|------------------------|
| type-I | V''(r*) > 0, V(r*) > 0 | Metastable (with barrier) |
| type-II | V''(r*) > 0, V(r*) < 0 | True minimum (below vacuum) |
| type-III | No minimum or V''(r*) ≤ 0 | Unstable |

**Key functions:**

- `classify_stability(V, r_star)`: Classify stability
- `StabilityResult`: Data class for classification results

```python
@dataclass
class StabilityResult:
    stability_type: str  # 'type-I', 'type-II', 'type-III'
    r_star: float        # Position of stable point
    V_star: float        # Potential value at stable point
    V_second: float      # Second derivative
    barrier_height: float  # Barrier height (type-I only)
```

## Usage

```python
from dppu.action import build_lagrangian, classify_stability

# Construct Lagrangian
L = build_lagrangian(R, N, kappa=1.0, theta_NY=0.5)

# Classify stability
result = classify_stability(V_eff, r_star=1.5)
print(f"Type: {result.stability_type}")
print(f"r* = {result.r_star:.3f}")
```

## Phase 1 Integration

Phase 1 scan results are saved as CSV and used for SD analysis in Phase 2:

```csv
topology,ny_variant,V,eta,theta,r0,delta_V,stability_type
S3,FULL,2.0,-1.0,0.0,1.234,0.056,type-I
```

## Dependencies

- [curvature](../curvature/README.md): Ricci scalar
- [torsion](../torsion/README.md): Nieh-Yan density

## Related Modules

- [scanning](../scanning/README.md): Parameter space scanning
