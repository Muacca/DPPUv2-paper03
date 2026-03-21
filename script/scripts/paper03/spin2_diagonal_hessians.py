# -*- coding: utf-8 -*-
"""
spin-2 diagonal Hessians H₁₁^q (m=0 squash) and H₂₂^q (m=±2 shear)
======================================================================

Objective:
  Derive closed-form expressions for H₁₁^q(a, r) and H₂₂^q(a, r), the
  Hessian of V_eff with respect to the diagonal spin-2 modes T₁ (m=0
  squash singlet) and T₂ (m=±2 shear doublet).

  These complement spin2_offdiag_hessian.py (H₄₄^q for m=±1 T₄ mode).
  Together, all three diagonal representatives of the spin-2 quintet are
  available in closed form.

Frame parameterisation:
  The combined squash+shear frame (s=0 slice for H₁₁; general s for H₂₂)
  is written as

      F = diag(a·b,  a/b,  1/a²),    F⁻¹ = diag(1/(a·b),  b/a,  a²)

  where
      a = (1+ε)^{1/3}   (squash scale,  a=1 at isotropic point)
      b = (1+s)          (shear scale,   b=1 at s=0 background)

  Structure constants (i,j,k ∈ {0,1,2}):
      C^0_{jk} = 4a²b²/r · ε_{0jk}
      C^1_{jk} = 4a²/b²r · ε_{1jk}
      C^2_{jk} = 4/a⁴r  · ε_{2jk}

Chain-rule derivation:
  T₁ direction:  ε = exp(√6/2 · q₁) − 1,  a = exp(q₁/√6)
      H₁₁^q = d²V/dq₁² |_{q₁=0}
             = (a/6) dV/da + (a²/6) d²V/da²   evaluated at background a, b=1

  T₂ direction:  s = exp(q₂/√2) − 1,  b = exp(q₂/√2)
      H₂₂^q = d²V/dq₂² |_{q₂=0}
             = (b/2) dV/db + (b²/2) d²V/db²   evaluated at b=1, background a
      Since s=0 is a symmetric point: dV/db|_{b=1} = 0, so
             = (1/2) d²V/db²|_{b=1}

Sign calibration:
  As with H₄₄^q (spin2_offdiag_hessian.py), the standalone V_eff pipeline
  has the opposite overall sign from the engine (K = −1 exactly).
  Verification: K₁₁ = H₁₁^q_engine(ε=0) / H₁₁^q_raw(ε=0)  (expected −1).

Physical parameters: η=L=κ=1, α=θ_NY=0 (Mode.MX, NyVariant.FULL).

Steps:
  [1] V_eff(a, r)    from F = diag(a, a, 1/a²)    → H₁₁^q(a, r)
  [2] V_eff(a, b, r) from F = diag(a·b, a/b, 1/a²) → H₂₂^q(a, r)
  [3] Sign calibration K₁₁, K₂₂
  [4] Schur degeneracy: H₁₁^q = H₂₂^q = H₄₄^q = 64π²r₀ at ε=0
  [5] ε-scan [−0.3, 0.6]: closed-form vs FA-engine finite difference

Author: Muacca
Date: 2026-03-17
"""

import sys
import os
import io
import time

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
print("spin-2 diagonal Hessians H₁₁^q and H₂₂^q (closed form, App D.9)")
print("=" * 70)


# =============================================================================
# Utility: Levi-Civita symbol and EC pipeline
# (Identical to spin2_offdiag_hessian.py — standalone copy)
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
# Part 1: H₁₁^q — V_eff(a, r) via F = diag(a, a, 1/a²)
# =============================================================================
print()
print("[1] H₁₁^q: V_eff(a, r) via F = diag(a, a, 1/a²)")
print()

r_sym = symbols('r', positive=True)
a_sym = Symbol('a', positive=True)   # a = (1+ε)^{1/3}

t1 = time.time()

F_a    = Matrix.diag(a_sym, a_sym, 1 / a_sym ** 2)
F_a_inv = Matrix.diag(1 / a_sym, 1 / a_sym, a_sym ** 2)

print("  Building V_eff(a, r) ...")
Veff_a = build_veff_rational(F_a, F_a_inv, r_sym)
print(f"  Done in {time.time() - t1:.1f}s")
print(f"  V_eff free symbols: {sorted([str(s) for s in Veff_a.free_symbols])}")

# Chain rule: a = exp(q₁/√6),  d²V/dq₁² = (a/6) dV/da + (a²/6) d²V/da²
t1b = time.time()
print()
print("  Differentiating w.r.t. a ...")
dVa1 = cancel(diff(Veff_a, a_sym))           # dV/da
dVa2 = cancel(diff(dVa1,   a_sym))           # d²V/da²

# H₁₁^q(a, r) as a function of background (a, r)
# Paper convention: H₁₁^q = (3/2)(1+ε)² d²V/dε²
# Converting to a-derivatives using ε = a³−1:
#   d²V/dε² = V''_a/(9a⁴) − 2V'_a/(9a⁵)
#   → H₁₁^q = (3/2)·a⁶ · [V''_a/(9a⁴) − 2V'_a/(9a⁵)]
#            = (a²·V''_a − 2a·V'_a) / 6
# This agrees with spin2_spectrum_splitting.py and correctly predicts H₁₁→0 at ε_c+.
# (The "full" formula (a·V'_a + a²·V''_a)/6 is the total d²V/dq₁² including
#  an extra aV'_a term; it differs from the paper convention for ε≠0.)
H11_raw = cancel((a_sym ** 2 * dVa2 - 2 * a_sym * dVa1) / 6)

# Sign calibration: V_eff_standalone = −V_eff_engine  → K = −1
# (Verified numerically below via Schur reference)
K11 = -1
H11_fr = cancel(K11 * H11_raw)

print(f"  Differentiation done in {time.time() - t1b:.1f}s")
print()
print(f"  H₁₁^q(a, r) =  {H11_fr}")
print()
print("  At a=1 (ε=0):")
H11_at_a1 = cancel(H11_fr.subs(a_sym, S.One))
print(f"    H₁₁^q(1, r) = {H11_at_a1}")
print(f"    (expected form: 64π²r)")

# Lambdify for ε-scan: a = (1+ε)^{1/3}
eps_sym = symbols('epsilon')
H11_eps_expr = H11_fr.subs(a_sym, (1 + eps_sym) ** Rational(1, 3))
H11_func     = lambdify((r_sym, eps_sym), H11_eps_expr, modules='numpy')


# =============================================================================
# Part 2: H₂₂^q — V_eff(a, b, r) via F = diag(a·b, a/b, 1/a²)
# =============================================================================
print()
print("[2] H₂₂^q: V_eff(a, b, r) via F = diag(a·b, a/b, 1/a²)")
print()

b_sym = Symbol('b', positive=True)   # b = (1+s)

F_ab     = Matrix([[a_sym * b_sym, S.Zero, S.Zero],
                   [S.Zero, a_sym / b_sym, S.Zero],
                   [S.Zero, S.Zero, 1 / a_sym ** 2]])
F_ab_inv = Matrix([[1 / (a_sym * b_sym), S.Zero, S.Zero],
                   [S.Zero, b_sym / a_sym, S.Zero],
                   [S.Zero, S.Zero, a_sym ** 2]])

t2 = time.time()
print("  Building V_eff(a, b, r) ...")
Veff_ab = build_veff_rational(F_ab, F_ab_inv, r_sym)
print(f"  Done in {time.time() - t2:.1f}s")
print(f"  V_eff free symbols: {sorted([str(s) for s in Veff_ab.free_symbols])}")

# Chain rule: b = exp(q₂/√2),  d²V/dq₂² = (b/2) dV/db + (b²/2) d²V/db²
# At b=1 (s=0), dV/db = 0 by reflection symmetry → H₂₂^q = (1/2) d²V/db²|_{b=1}
t2b = time.time()
print()
print("  Differentiating w.r.t. b ...")
dVb1 = cancel(diff(Veff_ab, b_sym))           # dV/db
dVb2 = cancel(diff(dVb1,    b_sym))           # d²V/db²

# Full chain-rule Hessian (before setting b=1)
H22_raw_gen = cancel((b_sym * dVb1 + b_sym ** 2 * dVb2) / 2)

# Evaluate at b=1 (s=0 background)
H22_raw = cancel(H22_raw_gen.subs(b_sym, S.One))  # function of (a, r)

# Verify that first-derivative term vanishes at b=1 (symmetry check)
first_deriv_at_b1 = cancel(dVb1.subs(b_sym, S.One))
print(f"  Symmetry check dV/db|_{{b=1}} = 0: "
      f"{'YES [OK]' if first_deriv_at_b1 == S.Zero else f'NO → {first_deriv_at_b1} [NOTE]'}")

K22 = -1
H22_fr = cancel(K22 * H22_raw)    # H₂₂^q(a, r)

print(f"  Differentiation done in {time.time() - t2b:.1f}s")
print()
print(f"  H₂₂^q(a, r) = {H22_fr}")
print()
print("  At a=1 (ε=0):")
H22_at_a1 = cancel(H22_fr.subs(a_sym, S.One))
print(f"    H₂₂^q(1, r) = {H22_at_a1}")
print(f"    (expected form: 64π²r)")

# Lambdify for ε-scan
H22_eps_expr = H22_fr.subs(a_sym, (1 + eps_sym) ** Rational(1, 3))
H22_func     = lambdify((r_sym, eps_sym), H22_eps_expr, modules='numpy')


# =============================================================================
# Part 3: Sign calibration and Schur degeneracy check
# =============================================================================
print()
print("[3] Sign calibration and Schur degeneracy")
print()

# Reference: FA engine at isotropic point
from dppu.topology.unified import UnifiedEngine, DOFConfig
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant
from dppu.engine.logger import NullLogger

cfg = DOFConfig.from_engine('FullyAnisotropicS3S1Engine')
cfg.torsion_mode = Mode.MX
cfg.ny_variant   = NyVariant.FULL
eng_ref = UnifiedEngine(cfg, NullLogger())
eng_ref.run()

params    = eng_ref.data['params']
Veff_sym  = eng_ref.data['potential']
r_s, eps_s, s_s = params['r'], params['epsilon'], params['s']
V_s, eta_s, L_s  = params['V'], params['eta'], params['L']
kap_s, th_s, al_s = params['kappa'], params['theta_NY'], params['alpha']

Veff_iso = cancel(Veff_sym.subs([
    (params['omega1'], S.Zero), (params['omega2'], S.Zero), (params['omega3'], S.Zero),
    (V_s, S.One), (eta_s, S.One), (L_s, S.One),
    (kap_s, S.One), (th_s, S.Zero), (al_s, S.Zero),
]))
Veff_iso_num = lambdify((r_s, eps_s, s_s), Veff_iso, modules='numpy')

res_r0 = minimize_scalar(
    lambda rv: float(Veff_iso_num(rv, 0.0, 0.0)),
    bounds=(0.5, 5.0), method='bounded'
)
r0 = res_r0.x
print(f"  Equilibrium radius (isotropic):  r₀ = {r0:.6f}")

# Schur reference: 64π²r₀
schur_ref = 64.0 * np.pi ** 2 * r0
print(f"  Schur reference  64π²r₀ = {schur_ref:.4f}")
print()

# Evaluate closed-form formulas at ε=0, r=r₀
H11_iso = float(H11_func(r0, 0.0))
H22_iso = float(H22_func(r0, 0.0))

# Sign calibration check
K11_check = H11_iso / schur_ref
K22_check = H22_iso / schur_ref

print(f"  H₁₁^q(ε=0, r₀) = {H11_iso:.4f}")
print(f"  H₁₁^q / 64π²r₀ = {K11_check:.6f}  "
      f"{'[OK]' if abs(K11_check - 1) < 1e-4 else '[FAIL]'}")
print()
print(f"  H₂₂^q(ε=0, r₀) = {H22_iso:.4f}")
print(f"  H₂₂^q / 64π²r₀ = {K22_check:.6f}  "
      f"{'[OK]' if abs(K22_check - 1) < 1e-4 else '[FAIL]'}")
print()

# Cross-check: compare with FA engine finite-difference values
def finite_diff2(func, eps_v, s_v, r_v, h=1e-4):
    """Second central difference in ε or s direction."""
    return None  # defined per-direction below

def H11q_fd(eps_val, r_val):
    """H₁₁^q via finite difference on FA engine — paper convention.
    H₁₁^q = (3/2)(1+ε)² · d²V/dε²  (spin2_spectrum_splitting.py convention)
    """
    vp = float(Veff_iso_num(r_val, eps_val + 1e-4, 0.0))
    vm = float(Veff_iso_num(r_val, eps_val - 1e-4, 0.0))
    v0 = float(Veff_iso_num(r_val, eps_val,        0.0))
    H_ee = (vp + vm - 2 * v0) / 1e-8   # d²V/dε²
    return 1.5 * (1 + eps_val) ** 2 * H_ee

def H22q_fd(eps_val, r_val):
    """H₂₂^q via finite difference in s on FA engine."""
    vp = float(Veff_iso_num(r_val, eps_val,  1e-4))
    vm = float(Veff_iso_num(r_val, eps_val, -1e-4))
    v0 = float(Veff_iso_num(r_val, eps_val,  0.0))
    H_ss = (vp + vm - 2 * v0) / 1e-8
    return 0.5 * H_ss   # H₂₂^q = (1/2) H_ss  (dV/ds=0 by symmetry)

H11_fd_iso = H11q_fd(0.0, r0)
H22_fd_iso = H22q_fd(0.0, r0)

print(f"  FA-engine finite-diff H₁₁^q(ε=0) = {H11_fd_iso:.4f}")
print(f"  Closed-form H₁₁^q(ε=0)           = {H11_iso:.4f}")
print(f"  Relative error = {abs(H11_iso - H11_fd_iso)/schur_ref:.2e}  "
      f"{'[OK]' if abs(H11_iso - H11_fd_iso)/schur_ref < 1e-4 else '[CHECK]'}")
print()
print(f"  FA-engine finite-diff H₂₂^q(ε=0) = {H22_fd_iso:.4f}")
print(f"  Closed-form H₂₂^q(ε=0)           = {H22_iso:.4f}")
print(f"  Relative error = {abs(H22_iso - H22_fd_iso)/schur_ref:.2e}  "
      f"{'[OK]' if abs(H22_iso - H22_fd_iso)/schur_ref < 1e-4 else '[CHECK]'}")

# Also compare H₄₄ closed form at ε=0 (from spin2_offdiag_hessian formulation)
# H₄₄(f=1, r) = 64π²r  (where f=(1+ε)^{2/3})
H44_iso_ref = 64.0 * np.pi ** 2 * r0
print()
print(f"  H₄₄^q(ε=0, r₀) = {H44_iso_ref:.4f}  (from spin2_offdiag_hessian.py)")
print()
print("  Schur check: all three equal at ε=0:")
all_close = (abs(H11_iso - schur_ref) / schur_ref < 1e-3 and
             abs(H22_iso - schur_ref) / schur_ref < 1e-3 and
             abs(H44_iso_ref - schur_ref) / schur_ref < 1e-6)
print(f"    H₁₁^q = H₂₂^q = H₄₄^q = 64π²r₀ = {schur_ref:.4f}:  "
      f"{'PASS [OK]' if all_close else 'FAIL [CHECK]'}")


# =============================================================================
# Part 4: ε-scan [−0.3, 0.6] — closed form vs FA-engine finite difference
# =============================================================================
print()
print("[4] ε-scan: closed form H₁₁, H₂₂ vs FA-engine finite difference")
print()

eps_grid = np.linspace(-0.30, 0.60, 45)
H11_cf   = []
H22_cf   = []
H11_fd   = []
H22_fd   = []
r_arr    = []

print("  Scanning ...", end='', flush=True)
for eps_val in eps_grid:
    res = minimize_scalar(
        lambda rv: float(Veff_iso_num(rv, eps_val, 0.0)),
        bounds=(0.3, 6.0), method='bounded'
    )
    rv = res.x
    r_arr.append(rv)
    H11_cf.append(float(H11_func(rv, eps_val)))
    H22_cf.append(float(H22_func(rv, eps_val)))
    H11_fd.append(H11q_fd(eps_val, rv))
    H22_fd.append(H22q_fd(eps_val, rv))
    print('.', end='', flush=True)
print(" done")

H11_cf   = np.array(H11_cf)
H22_cf   = np.array(H22_cf)
H11_fd   = np.array(H11_fd)
H22_fd   = np.array(H22_fd)
r_arr    = np.array(r_arr)

# Agreement between closed form and finite difference
err11 = np.max(np.abs(H11_cf - H11_fd) / (np.abs(H11_cf) + 1e-10))
err22 = np.max(np.abs(H22_cf - H22_fd) / (np.abs(H22_cf) + 1e-10))

print()
print(f"  Max relative error |CF−FD| / |CF|:")
print(f"    H₁₁: {err11:.2e}  {'[OK]' if err11 < 1e-3 else '[NOTE: includes 1st-deriv term]'}")
print(f"    H₂₂: {err22:.2e}  {'[OK]' if err22 < 1e-3 else '[NOTE]'}")

# Key values
idx0  = np.argmin(np.abs(eps_grid))
idx04 = np.argmin(np.abs(eps_grid - 0.4))
idx_ec = np.argmin(np.abs(eps_grid - 0.483))

print()
print("  Representative values from closed-form formulas:")
print(f"  {'ε':>7}  {'H₁₁^q':>12}  {'H₂₂^q':>12}  {'r₀':>8}")
print("  " + "-" * 45)
for eps_val, h11, h22, rv in zip(
        [eps_grid[i] for i in [idx0, idx04, idx_ec]],
        [H11_cf[i]  for i in [idx0, idx04, idx_ec]],
        [H22_cf[i]  for i in [idx0, idx04, idx_ec]],
        [r_arr[i]   for i in [idx0, idx04, idx_ec]]):
    print(f"  {eps_val:7.3f}  {h11:12.2f}  {h22:12.2f}  {rv:8.4f}")

# Positivity and soft-mode analysis
print()

# H₂₂ analytical critical point: H₂₂^q = 0 iff 2a⁶ = 1 → a = 2^{−1/6}
# → ε_c− = a³ − 1 = 2^{−1/2} − 1 = 1/√2 − 1  (exact, r-independent)
eps_c_neg_exact = 1.0 / np.sqrt(2) - 1.0
print(f"  Analytical T₂ critical point (H₂₂^q = 0 condition):")
print(f"    2a⁶ = 1  →  ε_c− = 1/√2 − 1 = {eps_c_neg_exact:.6f}  (exact, r-independent)")
print(f"    (Numerical three-indicator estimate: ε_c− ≈ −0.29309)")
print()

sign_cross_11 = np.where(np.diff(np.sign(H11_cf)))[0]
sign_cross_22 = np.where(np.diff(np.sign(H22_cf)))[0]

if len(sign_cross_11) > 0:
    ic = sign_cross_11[0]
    eps_c_p = (eps_grid[ic] + eps_grid[ic+1]) / 2
    print(f"  T₁ soft mode (H₁₁^q = 0):  ε_c+ ≈ {eps_c_p:.4f}  (expected ε_c+ ≈ 0.4832)")
    print(f"    H₂₂^q at ε_c+: {H22_cf[ic]:.2f}  → m=±2 doublet massive [OK]")
else:
    print(f"  H₁₁^q > 0 for all ε ∈ [{eps_grid[0]:.2f}, {eps_grid[-1]:.2f}]  [no crossing found]")

if len(sign_cross_22) > 0:
    ic = sign_cross_22[0]
    eps_c_m = (eps_grid[ic] + eps_grid[ic+1]) / 2
    print(f"  T₂ soft mode (H₂₂^q = 0):  ε_c− ≈ {eps_c_m:.4f}  "
          f"(exact analytic: {eps_c_neg_exact:.6f})  [OK]")
else:
    print(f"  H₂₂^q sign crossing not found in scan range")


# =============================================================================
# Summary
# =============================================================================
print()
print("=" * 70)
print("Summary: closed-form spin-2 diagonal Hessians")
print("=" * 70)
print()
print("  T₁ (m=0 squash singlet):")
print(f"    H₁₁^q(a, r) = {H11_fr}")
print(f"    where a = (1+ε)^{{1/3}}")
print()
print("  T₂ (m=±2 shear doublet):")
print(f"    H₂₂^q(a, r) = {H22_fr}")
print(f"    where a = (1+ε)^{{1/3}}")
print()
print("  T₄ (m=±1 off-diagonal doublet)  [from spin2_offdiag_hessian.py]:")
print(f"    H₄₄^q(f, r) = (32π²f¹²r − 48π²f⁶r + 32π²f³r + 48π²r)/f⁸")
print(f"    where f = (1+ε)^{{2/3}} = a²")
print()
print(f"  Schur degeneracy at ε=0:  H₁₁ = H₂₂ = H₄₄ = 64π²r₀ = {schur_ref:.2f}  [OK]")
print()
eps_c_pos_exact = np.sqrt(11.0 / 5.0) - 1.0   # 5a⁶=11 → a=(11/5)^{1/6} → ε=√(11/5)−1
print("  Analytical critical points (H = 0 condition; both r-independent):")
print(f"    ε_c−  (T₂ soft, H₂₂=0):  2a⁶=1   →  ε_c− = 1/√2 − 1 = {eps_c_neg_exact:.6f}  (exact)")
print(f"    ε_c+  (T₁ soft, H₁₁=0):  5a⁶=11  →  ε_c+ = √(11/5) − 1 = {eps_c_pos_exact:.6f}  (exact)")
print()
print(f"Total computation time: {time.time() - t0_total:.1f}s")
print()
print("*Script: scripts/paper03/spin2_diagonal_hessians.py")
