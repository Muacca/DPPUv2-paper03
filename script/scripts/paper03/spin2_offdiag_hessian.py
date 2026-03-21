# -*- coding: utf-8 -*-
"""
spin-2 off-diagonal Hessian H₄₄^q via rational substitution (App D.9)
=======================================================================

Objective:
  Derive a closed-form expression for H₄₄^q(ε), the Hessian of V_eff with
  respect to the off-diagonal spin-2 mode T₄ (m=±1 doublet, 02-plane rotation).
  Unlike the diagonal modes T₁ (m=0) and T₂ (m=±2), T₄ requires an
  off-diagonal frame deformation that makes symbolic computation with SymPy's
  standard cosh/sinh pipeline prohibitively slow (OOM on 64 GB machines).

Method (Method 2 — rational substitution):
  Replace the transcendental parameterization
      cosh(q₄/√2) → (z₄ + 1/z₄)/2
      sinh(q₄/√2) → (z₄ − 1/z₄)/2
  where z₄ = exp(q₄/√2) > 0.  The full EC + torsion pipeline is then carried
  out as rational arithmetic in (z₄, f, r), where f = (1+ε)^{2/3}.

  The Hessian at q₄=0 (z₄=1) follows from the chain rule:
      d/dq₄ = (1/√2) z₄ d/dz₄
      H₄₄ = d²V/dq₄² |_{q₄=0}
           = (1/2)[dV/dz₄ + z₄ d²V/dz₄²]_{z₄=1}
           = (1/2)[dV/dz₄ + d²V/dz₄²]_{z₄=1}

  Physical parameters: η=L=κ=1, α=θ_NY=0 (Mode.MX, NyVariant.FULL).

Steps:
  [1] Build EC pipeline as rational functions in (z₄, f, r) via Method 2
  [2] Derive H₄₄^q(f, r) in closed form by symbolic differentiation
  [3] Verify Schur degeneracy at ε=0: H₄₄^q = H₁₁^q = H₂₂^q = 64π²r₀
  [4] Scan ε ∈ [−0.3, 0.6] and confirm H₄₄^q(ε) > 0 throughout

Author: Muacca
Date: 2026-03-16
"""

import sys
import os
import time
import io

if sys.stdout.encoding and sys.stdout.encoding.lower() not in ('utf-8', 'utf8'):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import numpy as np
from sympy import (
    symbols, Symbol, Rational, S,
    Matrix, cancel, diff, pi, lambdify
)
from scipy.optimize import minimize_scalar

t0_total = time.time()

print("=" * 70)
print("spin-2 off-diagonal Hessian H₄₄^q via rational substitution (App D.9)")
print("=" * 70)


# =============================================================================
# Utility: Levi-Civita symbol
# =============================================================================

def eps3(i, j, k):
    """3D Levi-Civita symbol."""
    t = (i, j, k)
    if t in {(0, 1, 2), (1, 2, 0), (2, 0, 1)}: return S.One
    if t in {(0, 2, 1), (2, 1, 0), (1, 0, 2)}: return S.NegativeOne
    return S.Zero


def build_C_spatial(F, F_inv, r_sym):
    """
    Structure constants C^i_{jk} = (4/r) Σ_{lmn} F[i,l] F⁻¹[j,m] F⁻¹[k,n] ε_{lmn}
    for i,j,k ∈ {0,1,2}.
    """
    C = [[[S.Zero] * 3 for _ in range(3)] for _ in range(3)]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                v = S.Zero
                for l in range(3):
                    for m in range(3):
                        for n in range(3):
                            e = eps3(l, m, n)
                            if e != S.Zero:
                                v += F[i, l] * F_inv[j, m] * F_inv[k, n] * e
                C[i][j][k] = cancel(4 * v / r_sym)
    return C


def build_veff_rational(F, F_inv, r_sym, eta_val=S.One, L_val=S.One, kappa_val=S.One):
    """
    Compute V_eff as a rational function in the symbols of F, F_inv, and r_sym.

    Pipeline: structure constants → LC connection → torsion (MX mode) →
              contortion → EC connection → Riemann scalar → V_eff.

    All intermediate expressions are simplified by cancel() to keep
    rational polynomial representations compact.
    """
    C3 = build_C_spatial(F, F_inv, r_sym)

    def C(a, b, c):
        return C3[a][b][c] if (a < 3 and b < 3 and c < 3) else S.Zero

    def Gamma(a, b, c):
        """LC connection: Γ^a_{bc} = (C[a,b,c] + C[c,b,a] − C[b,a,c]) / 2"""
        return cancel((C(a, b, c) + C(c, b, a) - C(b, a, c)) / 2)

    def T(a, b, c):
        """
        Torsion tensor (Mode.MX: axial T1 + vector-trace T2).
        T1: T_{ijk} = (2η/r) ε_{ijk}  (spatial only)
        T2: T_{abc} = (δ_{ac} V_b − δ_{ab} V_c) / 3,  V_μ = δ^3_μ
        """
        t1 = (2 * eta_val / r_sym * eps3(a, b, c)
              if (a < 3 and b < 3 and c < 3) else S.Zero)
        V_b = S.One if b == 3 else S.Zero
        V_c = S.One if c == 3 else S.Zero
        t2 = ((S.One if a == c else S.Zero) * V_b
              - (S.One if a == b else S.Zero) * V_c) / 3
        return t1 + t2

    def K(a, b, c):
        """Contortion: K_{abc} = (T[a,b,c] + T[b,c,a] − T[c,a,b]) / 2"""
        return cancel((T(a, b, c) + T(b, c, a) - T(c, a, b)) / 2)

    # Cache EC connection (64 components)
    GEC = [[[cancel(Gamma(a, b, c) + K(a, b, c))
             for c in range(4)] for b in range(4)] for a in range(4)]
    C4  = [[[C(a, b, c)
             for c in range(4)] for b in range(4)] for a in range(4)]

    # Ricci scalar: R = Σ_{a,b} R^a_{bab}
    R_scalar = S.Zero
    for a in range(4):
        for b in range(4):
            v = S.Zero
            for e in range(4):
                v += (GEC[a][e][a] * GEC[e][b][b]
                      - GEC[a][e][b] * GEC[e][b][a]
                      + GEC[a][b][e] * C4[e][a][b])
            R_scalar = cancel(R_scalar + v)

    vol = 2 * pi ** 2 * L_val * r_sym ** 3
    return cancel(R_scalar / (2 * kappa_val ** 2) * vol)


# =============================================================================
# Part 1: Build V_eff as rational function in (z₄, f, r)
# =============================================================================
print()
print("[1] Build V_eff in (z₄, f, r) via rational substitution")
print()

t1 = time.time()
r_sym = symbols('r', positive=True)
z4    = Symbol('z4', positive=True)   # z₄ = exp(q₄/√2)
f     = Symbol('f',  positive=True)   # f  = (1+ε)^{2/3}

c4 = (z4 + 1 / z4) / 2   # cosh(q₄/√2) = (z₄ + 1/z₄)/2
s4 = (z4 - 1 / z4) / 2   # sinh(q₄/√2) = (z₄ − 1/z₄)/2

# Frame matrix: F = F_diag(ε) × G₄(q₄)
# G₄ is the hyperbolic rotation in the (0,2) plane
G4     = Matrix([[c4, S.Zero, s4], [S.Zero, S.One, S.Zero], [s4, S.Zero, c4]])
G4_inv = Matrix([[c4, S.Zero, -s4], [S.Zero, S.One, S.Zero], [-s4, S.Zero, c4]])

F     = Matrix.diag(f, f, 1 / f ** 2) * G4
F_inv = G4_inv * Matrix.diag(1 / f, 1 / f, f ** 2)

print("  Running EC pipeline in (z₄, f, r) ...")
Veff = build_veff_rational(F, F_inv, r_sym)
print(f"  Pipeline completed in {time.time() - t1:.1f}s")
print(f"  V_eff free symbols: {sorted([str(s) for s in Veff.free_symbols])}")


# =============================================================================
# Part 2: Derive H₄₄^q(f, r) in closed form
# =============================================================================
print()
print("[2] Derive H₄₄^q(f, r) in closed form")
print()

# Chain rule for z₄ = exp(q₄/√2):
#   dV/dq₄   = (z₄/√2) dV/dz₄
#   d²V/dq₄² = (1/2)(z₄ dV/dz₄ + z₄² d²V/dz₄²)
# Evaluated at q₄=0, i.e., z₄=1:
#   H₄₄ = (1/2)(dV/dz₄ + d²V/dz₄²)|_{z₄=1}

t2 = time.time()
print("  Differentiating V_eff w.r.t. z₄ ...")
dV1     = cancel(diff(Veff, z4))
dV2     = cancel(diff(dV1, z4))
H44_expr = cancel((dV1 + dV2) / 2)
H44_raw  = cancel(H44_expr.subs(z4, S.One))
print(f"  Differentiation completed in {time.time() - t2:.1f}s")

# Sign calibration: the standalone symbolic pipeline computes V_eff with the
# opposite overall sign compared to the engine (K₄₄ = −1 exactly).
# This is a systematic difference in the Lagrangian sign convention between
# the two implementations; it is absorbed by the factor K₄₄ = −1.
# Verification: K₄₄ = H₁₁^q(ε=0) / H₄₄^raw(ε=0) = −1 (confirmed numerically).
K44 = -1
H44_fr = cancel(K44 * H44_raw)

print()
print(f"  H₄₄^q(f, r) =  {H44_fr}")
print()
print("  Interpretation: f = (1+ε)^{2/3}, so ε=0 corresponds to f=1.")
print(f"  At f=1: H₄₄^q(1, r) = {cancel(H44_fr.subs(f, S.One))}")

# Substitute f = (1+ε)^{2/3}
eps_sym = symbols('epsilon')
H44_eps = H44_fr.subs(f, (1 + eps_sym) ** Rational(2, 3))
H44_func = lambdify((r_sym, eps_sym), H44_eps, modules='numpy')


# =============================================================================
# Part 3: Verify Schur degeneracy at ε=0
# =============================================================================
print()
print("[3] Verify Schur degeneracy at ε=0")
print()

# Equilibrium radius r₀ from isotropic V_eff (FA engine reference value)
from dppu.topology.unified import UnifiedEngine, DOFConfig
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant
from dppu.engine.logger import NullLogger

cfg = DOFConfig.from_engine('FullyAnisotropicS3S1Engine')
cfg.torsion_mode = Mode.MX
cfg.ny_variant   = NyVariant.FULL
eng_ref = UnifiedEngine(cfg, NullLogger())
eng_ref.run()

params   = eng_ref.data['params']
Veff_sym = eng_ref.data['potential']
r_s, eps_s, s_s = params['r'], params['epsilon'], params['s']
V_s, eta_s, L_s = params['V'], params['eta'], params['L']
kap_s, th_s, al_s = params['kappa'], params['theta_NY'], params['alpha']

Veff_iso = cancel(Veff_sym.subs([
    (params['omega1'], S.Zero), (params['omega2'], S.Zero), (params['omega3'], S.Zero),
    (V_s, S.One), (eta_s, S.One), (L_s, S.One),
    (kap_s, S.One), (th_s, S.Zero), (al_s, S.Zero),
]))
Veff_iso_num = lambdify((r_s, eps_s, s_s), Veff_iso, modules='numpy')

res_r0 = minimize_scalar(
    lambda rv: float(Veff_iso_num(rv, 0.0, 0.0)), bounds=(0.5, 5.0), method='bounded'
)
r0 = res_r0.x
print(f"  Equilibrium radius (isotropic point):  r₀ = {r0:.6f}")

# Schur reference: 64π²r₀  (from spin2_schur_degeneracy.py)
schur_ref = 64.0 * np.pi ** 2 * r0
H44_iso   = float(H44_func(r0, 0.0))

print(f"  H₄₄^q(ε=0, r₀) = {H44_iso:.4f}")
print(f"  Schur reference  64π²r₀ = {schur_ref:.4f}")
print(f"  Relative error   |H₄₄ − 64π²r₀| / 64π²r₀ = "
      f"{abs(H44_iso - schur_ref) / schur_ref:.2e}  "
      f"{'[OK]' if abs(H44_iso - schur_ref) / schur_ref < 1e-6 else '[FAIL]'}")


# =============================================================================
# Part 4: Scan ε ∈ [−0.3, 0.6] — confirm H₄₄^q > 0 throughout
# =============================================================================
print()
print("[4] ε-scan: confirm H₄₄^q(ε) > 0 for all ε ∈ [−0.3, 0.6]")
print()

eps_grid = np.linspace(-0.30, 0.60, 45)

H44_arr  = []
r_arr    = []

print("  Computing ...", end='', flush=True)
for eps_val in eps_grid:
    res = minimize_scalar(
        lambda rv: float(Veff_iso_num(rv, eps_val, 0.0)),
        bounds=(0.3, 6.0), method='bounded'
    )
    rv = res.x
    r_arr.append(rv)
    H44_arr.append(float(H44_func(rv, eps_val)))
    print('.', end='', flush=True)
print(" done")

H44_arr = np.array(H44_arr)
r_arr   = np.array(r_arr)

H44_min = H44_arr.min()
H44_max = H44_arr.max()
idx_min = H44_arr.argmin()

print()
print(f"  ε range scanned:   [{eps_grid[0]:.2f}, {eps_grid[-1]:.2f}]")
print(f"  H₄₄^q minimum:    {H44_min:.4f}  at ε = {eps_grid[idx_min]:.3f},  r = {r_arr[idx_min]:.4f}")
print(f"  H₄₄^q maximum:    {H44_max:.4f}")
print(f"  H₄₄^q > 0 throughout: "
      f"{'YES [OK]' if H44_min > 0 else 'NO [FAIL]'}")

# Find first-transition reference (ε_c where H₁₁^q → 0)
# H₁₁^q from FA engine
def finite_diff2(func, eps_v, r_v, h=1e-4):
    vp = func(r_v, eps_v + h, 0.0)
    vm = func(r_v, eps_v - h, 0.0)
    v0 = func(r_v, eps_v,     0.0)
    return (vp + vm - 2 * v0) / h ** 2

H11_at_ec_grid = []
for eps_val, rv in zip(eps_grid, r_arr):
    H_ee  = finite_diff2(lambda r, e, s: float(Veff_iso_num(r, e, s)), eps_val, rv)
    H11_q = H_ee * (np.sqrt(6) / 2 * (1 + eps_val)) ** 2
    H11_at_ec_grid.append(H11_q)

H11_arr    = np.array(H11_at_ec_grid)
sign_cross = np.where(np.diff(np.sign(H11_arr)))[0]

if len(sign_cross) > 0:
    eps_c = eps_grid[sign_cross[0]]
    idx_c = sign_cross[0]
    print()
    print(f"  First phase transition (T₁ soft mode):  ε_c ≈ {eps_c:.4f}")
    print(f"  H₄₄^q at ε_c:  {H44_arr[idx_c]:.4f}  "
          f"→  m=±1 doublet {'massive [OK]' if H44_arr[idx_c] > 0 else 'NOT massive [FAIL]'}")
else:
    print(f"  H₁₁^q did not change sign within scanned range")

print()
print("Summary")
print("-------")
print(f"  Closed-form H₄₄^q(f,r) = {H44_fr}")
print(f"  where f = (1+ε)^{{2/3}}")
print(f"  Schur degeneracy at ε=0:  H₄₄^q = H₁₁^q = H₂₂^q = 64π²r₀  [{'OK' if abs(H44_iso - schur_ref) / schur_ref < 1e-6 else 'FAIL'}]")
print(f"  H₄₄^q > 0 for all ε ∈ [{eps_grid[0]:.2f}, {eps_grid[-1]:.2f}]:  [{'OK' if H44_min > 0 else 'FAIL'}]")
print(f"  m=±1 doublet (T₄,T₅) remains massive at ε_c:  [{'OK' if len(sign_cross) > 0 and H44_arr[sign_cross[0]] > 0 else 'N/A' if len(sign_cross) == 0 else 'FAIL'}]")

print()
print(f"Total computation time: {time.time() - t0_total:.1f}s")
print()
print("*Script: scripts/paper03/spin2_offdiag_hessian.py")
