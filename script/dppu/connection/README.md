# Connection Layer

⇒ [日本語](README_ja.md)

Module group responsible for connection computation.

## Overview

Computes Levi-Civita connection, contortion, and Einstein-Cartan connection.

## Modules

### levi_civita.py

Computes the torsion-free Levi-Civita connection.

**Key functions:**

- `compute_christoffel_symbols(structure_constants, metric)`: Computation via general Koszul formula

**Koszul formula:**

```
Γ^a_{bc} = (1/2)(C^a_{bc} + η^{ad}η_{be}C^e_{dc} - η^{ad}η_{ce}C^e_{bd})
```

For bi-invariant metrics, this simplifies to:
```
Γ^a_{bc} = (1/2)C^a_{bc}
```

### contortion.py

Computes the contortion tensor arising from torsion.

**Key functions:**

- `compute_contortion(torsion, metric)`: Compute contortion tensor

**Definition:**

```
K_{abc} = (1/2)(T_{abc} + T_{bca} - T_{cab})
```

### ec_connection.py

Integrates the Einstein-Cartan connection (LC + contortion).

**Key functions:**

- `compute_ec_connection(christoffel, contortion)`: Compute EC connection

**Definition:**

```
Γ^{EC,a}_{bc} = Γ^{LC,a}_{bc} + K^a_{bc}
```

## Usage

```python
from dppu.connection import compute_christoffel_symbols, compute_contortion

# Levi-Civita connection
christoffel = compute_christoffel_symbols(C, eta)

# Contortion
K = compute_contortion(T, eta)

# EC connection
Gamma_EC = christoffel + K
```

## Metric Compatibility

The EC connection satisfies metric compatibility:
```
∇_a η_{bc} = 0
```

This is verified automatically.

## Dependencies

- [geometry](../geometry/README.md): Metric tensor
- [torsion](../torsion/README.md): Torsion tensor

## Related Modules

- [curvature](../curvature/README.md): Compute curvature from connection
