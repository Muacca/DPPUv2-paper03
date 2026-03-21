# Topology Layer

⇒ [日本語](README_ja.md)

Module group providing topology-specific computation engines.

## Overview

The topology layer implements the geometry of $M^3 \times S^1$ manifolds ($M^3$ = S³, T³, or Nil³) within the 15-step EC+NY computation pipeline. Paper 03 consolidates all topology engines into a unified configurable interface (`UnifiedEngine` + `DOFConfig`), replacing the individual `S3S1Engine`, `T3S1Engine`, `Nil3S1Engine` classes from paper 02.

## Class Hierarchy

```
BaseFrameEngine          (dppu/engine/pipeline.py)
    └── TopologyEngine   (base_topology.py)   ← abstract
            ├── S3Geometry    (s3.py)
            ├── T3Geometry    (t3.py)
            └── Nil3Geometry  (nil3.py)
```

`UnifiedEngine` (in `unified.py`) is a factory that selects the appropriate subclass based on a `DOFConfig` instance.

---

## Modules

### unified.py — UnifiedEngine and DOFConfig

Entry point for all paper 03 computations.

**Key classes:**

- `UnifiedEngine`: Factory creating the appropriate topology engine from a `DOFConfig`
- `DOFConfig`: Configuration dataclass specifying active geometric DOFs
- `TopologyType`: Enum — `S3` / `T3` / `NIL3`
- `FiberMode`: Enum — `NONE` / `TWIST` / `MIXING` / `BOTH`

**DOFConfig parameters:**

| Parameter | Type | Meaning |
|---|---|---|
| `topology` | `TopologyType` | Base topology |
| `enable_squash` | `bool` | Activate ε deformation |
| `enable_shear` | `bool` | Activate s (T₂) shear |
| `fiber_mode` | `FiberMode` | Fiber DOF selection |
| `isotropic_twist` | `bool` | Use single ω for all three directions |
| `torsion_mode` | `Mode` | `AX` / `VT` / `MX` |
| `ny_variant` | `NyVariant` | Nieh-Yan variant |
| `enable_velocity` | `bool` | Activate velocity symbols for G metric |

**Usage:**

```python
from dppu.topology.unified import UnifiedEngine, DOFConfig, TopologyType, FiberMode
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant

cfg = DOFConfig(
    topology=TopologyType.S3,
    enable_squash=False,
    fiber_mode=FiberMode.BOTH,
    isotropic_twist=False,
    torsion_mode=Mode.MX,
    ny_variant=NyVariant.FULL,
)
engine = UnifiedEngine(cfg)
engine.run()

Veff = engine.data['potential']
fp   = engine.get_free_params()
```

**Preset names (backward compatibility):**

```python
cfg = DOFConfig.from_engine('S3S1Engine')    # S³, squash=True, fiber=NONE
cfg = DOFConfig.from_engine('T3S1Engine')    # T³, fiber=NONE
cfg = DOFConfig.from_engine('Nil3S1Engine')  # Nil³, squash=True, fiber=NONE
```

See the docstring of `unified.py` for the complete equivalence table (13 legacy engines).

---

### base_topology.py — TopologyEngine Abstract Base

Abstract base class implementing the Template Method pattern for steps E4.1/E4.2.

**Abstract methods (subclasses must implement):**

- `_build_radial_and_deformation_params(params)` — define r/R, ε, s symbols
- `_build_structure_constants(params, C)` — fill structure constant array C
- `_compute_volume(params)` — return volume factor expression

**Common concrete helpers:**

- `_build_fiber_params(params)` — add ω / δ / cd / sd symbols (TWIST/MIXING/BOTH modes)
- `_add_s3_twist_C(C, ...)` — embed twist terms $C^3_{jk}$
- `_apply_s3_mixing_rotation_to_C(C, params)` — apply mixing rotation $M(\delta_0,\delta_1,\delta_2)$
- `get_free_params()` — return dict of active SymPy Symbols
- `get_riemann_lambdified()` — return lambdified $R_{abcd}$ for numerical scanning

---

### s3.py — S³×S¹ (S3Geometry)

**Mathematical structure:**

- Lie group: SU(2)
- Structure constants: $C^i_{jk} = (4/r)\,\lambda_i\,\varepsilon_{ijk}$ (scaling by ε, s)
- Metric: bi-invariant
- Background curvature: $R_{\rm LC} = 24/r^2$

**Supported DOFs:** ε (squash), s (shear), q₃,q₄,q₅ (off-diagonal shear), TWIST (ω₁,ω₂,ω₃), MIXING (δ₀,δ₁,δ₂)

**Volume:** $V = 2\pi^2 L r^3$

**Usage:**

```python
from dppu.topology.s3 import S3Geometry
from dppu.topology.unified import DOFConfig, TopologyType

cfg = DOFConfig.from_engine('S3S1Engine')
engine = S3Geometry(cfg)
engine.run()
```

---

### t3.py — T³×S¹ (T3Geometry)

**Mathematical structure:**

- Lie group: U(1)³ (Abelian)
- Structure constants: all zero
- Metric: flat
- Background curvature: $R_{\rm LC} = 0$

**Volume:** $V = (2\pi)^4 L r^3$ (isotropic scaling $R_1 = R_2 = R_3 = r$)

---

### nil3.py — Nil³×S¹ (Nil3Geometry)

**Mathematical structure:**

- Lie group: Heisenberg group
- Structure constants: $[E_0, E_1] = (1/R)E_2$
- Metric: left-invariant (**NOT bi-invariant**)
- Background curvature: $R_{\rm LC} = -1/(2R^2)$

**Important:** Nil³ does not admit a bi-invariant metric, so the general Koszul formula is used (see CONVENTIONS §5).

**Volume:** $V = (2\pi)^4 L R^3$

---

## Topology Comparison

| Property | S³×S¹ | T³×S¹ | Nil³×S¹ |
|----------|-------|-------|---------|
| Structure constants | $\varepsilon_{ijk}$ | 0 | $[E_0,E_1]=E_2$ |
| Background curvature | +24/r² | 0 | −1/(2R²) |
| bi-invariant | Yes | Yes | **No** |
| Koszul formula | Simplified | Trivial | General |
| Supported squash | Yes | No | Yes |
| Supported fiber | TWIST/MIXING/BOTH | TWIST | TWIST |

## Dependencies

- [engine](../engine/README.md): `BaseFrameEngine`
- [geometry](../geometry/README.md): Metric definitions
- [connection](../connection/README.md): Connection computation
- [curvature](../curvature/README.md): Curvature computation
- [torsion](../torsion/README.md): Torsion construction
- [action](../action/README.md): Action and stability
