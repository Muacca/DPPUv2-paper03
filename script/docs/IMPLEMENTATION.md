# Numerical implementation guide

[日本語](IMPLEMENTATION_ja.md) / [Conventions](CONVENTIONS.md) / [Reproduction](../README.md)

## Paper calculation

Run Python from `script`, as for the command-line entry point. `paper_results`
uses NumPy arrays and closed formulas.
`cartan` independently assembles C, geometric T, the LC and EC connections,
and the constant-component curvature before contraction. `potential` evaluates
the closed result (M3).

```python
import numpy as np
from paper_results.model import cartan, potential
from paper_results.circle import joint_bulk

geometry = cartan(r=3., L=1., eta=1., V=1.)
U = potential(3., eta=1., V=1.)
assert np.isclose(-geometry["volume"] * geometry["R_scalar"] / 2, U)
assert geometry["torsion_residual"] < 1e-12
bulk = joint_bulk(r=3., M=np.eye(3), M_prime=np.zeros((3, 3)))
assert np.isclose(bulk["bulk"], U)
```

The arrays are floating point. Numerical residuals should be compared with a
tolerance; the manuscript supplies the exact identities. `cartan` is static.
For periodic profiles use `circle.joint_bulk`, which implements (M6) after
integrating total derivatives. Integrate its `bulk` value over y in [0,1].
Pointwise values of this bulk function need not equal the unintegrated local
curvature density.

| Function | Return value |
|---|---|
| `cartan` | Dictionary including coframe, C, Gamma, T, R, R_LC, R_scalar, TT, REE, NY, P, B, volume and torsion_residual |
| `scalar_data` | Gradient and 4 by 4 Hessian in (r,L,eta,V) order |
| `shape_minimum` | Nonnegative doublet amplitude and minimum shape function |
| `twist_bulk` | Static twist contribution, derivative contribution and 3 by 3 velocity Hessian |
| `scalar_bulk` | Round scalar bulk density |
| `joint_bulk` | static, velocity_quadratic, velocity_linear, bulk, gamma, lapse, shift |
| `diagonal_cross_coefficients` | G_rw, G_sw, G_rs, G_rL, G_sL |
| `coupled_periodic_family` | Frequency series, static/derivative/total arrays, C_A and round reference value |

The model file checks positive scales in its principal action functions.
Use real finite arrays and determinant-one shapes; `joint_bulk` additionally
checks array dimensions and the tangent condition for M_prime.
When using matrix exponentials, the shear basis is symmetric and tracefree.

## Data and drawing

`data.recompute(output)` writes new fixed-grid datasets.
`data.load_data(source)` checks the schema and saved data hashes, returning
CSV arrays and scalar JSON. `figures.render(source,output)` draws ten pairs
from these arrays. The manifest stores parameters, grids and hashes of
model.py, circle.py and the data. Loading checks the data hashes; it does not
automatically compare the current source files with the recorded source hashes.
Explicit recomputation creates a new source/data association.
