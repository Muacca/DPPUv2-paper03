"""
DPPU Connection Layer
=====================

Provides connection-related computations for Einstein-Cartan gravity.

Modules:
- levi_civita: Levi-Civita (torsion-free) connection via Koszul formula
- contortion: Contortion tensor from torsion
- ec_connection: Einstein-Cartan connection = LC + contortion

Connection Types:
- Levi-Civita: Γ^a_{bc} (torsion-free, metric-compatible)
- Einstein-Cartan: ω^a_{bc} = Γ^a_{bc} + K^a_{bc}

Where:
- K^a_{bc} is the contortion tensor
- T^a_{bc} is the torsion tensor
- K^a_{bc} = (1/2)(T^a_{bc} + T^b_{ca} - T^c_{ab})
"""

from .levi_civita import compute_christoffel_frame
from .contortion import compute_contortion
from .ec_connection import compute_ec_connection

__all__ = [
    'compute_christoffel_frame',
    'compute_contortion',
    'compute_ec_connection',
]
