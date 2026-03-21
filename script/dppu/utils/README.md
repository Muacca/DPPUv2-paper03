# Utils Layer

⇒ [日本語](README_ja.md)

Module group providing shared utilities.

## Overview

Common functionality used across multiple modules: Levi-Civita symbol, symbolic computation helpers, and visualization utilities.

## Modules

### levi_civita.py

Levi-Civita symbol (totally antisymmetric tensor).

**Key functions:**

- `epsilon_symbol(i, j, k)`: 3D Levi-Civita symbol
- `levi_civita_4d(a, b, c, d)`: 4D Levi-Civita symbol

**3D:**

```python
from dppu.utils import epsilon_symbol

eps = epsilon_symbol(0, 1, 2)  # = +1
eps = epsilon_symbol(1, 0, 2)  # = -1
eps = epsilon_symbol(0, 0, 2)  # = 0
```

**4D:**

```python
from dppu.utils import levi_civita_4d

eps = levi_civita_4d(0, 1, 2, 3)  # = +1
eps = levi_civita_4d(1, 0, 2, 3)  # = -1
```

**Properties:**

- Totally antisymmetric: sign flips under any 2-index swap
- Zero for repeated indices
- Normalization: ε₀₁₂₃ = +1 (or ε₀₁₂ = +1)

### symbolic.py

Helper functions for symbolic computation.

**Key functions:**

- `prove_zero(expr)`: Prove an expression is zero
- `find_nonzero_witness(expr, params)`: Search for non-zero counterexamples
- `generate_test_points(param_ranges, n_points)`: Generate test points

**Zero proof:**

```python
from dppu.utils import prove_zero
import sympy as sp

x = sp.Symbol('x')
expr = x**2 - x**2

is_zero, proof_type = prove_zero(expr)
# is_zero=True, proof_type='symbolic'
```

**Counterexample search:**

```python
from dppu.utils import find_nonzero_witness

expr = some_complex_expression
witness = find_nonzero_witness(expr, {'r': (0.1, 5.0), 'eta': (-5, 5)})

if witness:
    print(f"Non-zero at: {witness['params']}")
    print(f"Value: {witness['value']}")
else:
    print("No non-zero witness found")
```

### visualization.py

Standardized plot styles and helper functions for paper figures.

**Key functions:**

- `set_style()`: Set global matplotlib style for high-quality figures (seaborn whitegrid, serif font, 300 dpi)
- `save_plot(filename, output_dir="output")`: Save plot with standard settings

**Usage:**

```python
from dppu.utils.visualization import set_style, save_plot
import matplotlib.pyplot as plt

set_style()
plt.plot(x, y)
plt.xlabel(r'$\varepsilon$')
plt.ylabel(r'$\mu^2$')
save_plot("spin2_spectrum.pdf", output_dir="output")
```

## Usage Examples

### Levi-Civita in Torsion Construction

```python
from dppu.utils import epsilon_symbol

# Axial torsion T^{ijk} = (2η/r) ε^{ijk}
def build_axial_torsion(eta, r):
    T = {}
    for i in range(3):
        for j in range(3):
            for k in range(3):
                T[(i,j,k)] = (2*eta/r) * epsilon_symbol(i, j, k)
    return T
```

### 4D Levi-Civita in Hodge Dual

```python
from dppu.utils import levi_civita_4d

# (*R)^{ab}_{cd} = (1/2) ε_{cdef} R^{ab,ef}
def hodge_dual(R):
    R_star = np.zeros((4,4,4,4))
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    for e in range(4):
                        for f in range(4):
                            R_star[a,b,c,d] += 0.5 * levi_civita_4d(c,d,e,f) * R[a,b,e,f]
    return R_star
```

## Mathematical Background

### Levi-Civita Symbol Definition

**3D:**
```
ε_{ijk} = { +1  (i,j,k) is an even permutation of (0,1,2)
          { -1  (i,j,k) is an odd permutation of (0,1,2)
          {  0  any index is repeated
```

**4D:**
```
ε_{abcd} = { +1  (a,b,c,d) is an even permutation of (0,1,2,3)
           { -1  (a,b,c,d) is an odd permutation of (0,1,2,3)
           {  0  any index is repeated
```

### Identities

```
ε_{ijk} ε_{ilm} = δ_{jl}δ_{km} - δ_{jm}δ_{kl}
ε_{abcd} ε_{abef} = 2(δ_{ce}δ_{df} - δ_{cf}δ_{de})
```

## Dependencies

- SymPy (symbolic computation)
- NumPy (numerical computation)
- matplotlib + seaborn (visualization)
- mpmath (high-precision arithmetic, optional)

## Related Modules

- [curvature](../curvature/README.md): Hodge dual computation
- [torsion](../torsion/README.md): Torsion construction
