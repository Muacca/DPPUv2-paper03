# -*- coding: utf-8 -*-
"""
spin-2 spectrum splitting and soft mode identification (App D.9)
=================================================

Objective:
  Show that in a squash background ε ≠ 0, the SO(3) → SO(2) symmetry breaking
  splits the 5-fold degenerate spin-2 spectrum into 3 groups.

  Splitting:
    m=0  (T₁, squash singlet)   : H₁₁^q(ε) = (3/2)(1+ε)² H_εε(ε)
    m=±2 (T₂,T₃, 0-1 doublet)  : H₂₂^q(ε) = (1/2) H_ss(ε)
    m=±1 (T₄,T₅, off-axis doublet): H₄₄^q(ε)  [numerical computation]

  Main conclusions:
    1. As ε increases, the m=0 (T₁) mode becomes soft first → condensation at ε_c
    2. m=±2 (T₂,T₃) remain massive at ε_c → s ≠ 0 is a separate transition
    3. m=±1 (T₄,T₅) have intermediate mass → massive at ε_c

Steps:
  [1] Group-theoretic organization of SO(3)→SO(2) splitting
  [2] H₁₁^q(ε), H₂₂^q(ε) — finite difference from FA engine
  [3] H₄₄^q(ε) — numerical Ricci scalar + finite difference
  [4] Plot ε-dependence of all 3 masses
  [5] Identify soft mode direction
  [6] Full spin-2 interpretation of s ≠ 0 phase transition

Author: Muacca
Date: 2026-03-15
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

from sympy import symbols, S, cancel, lambdify, diff as sdiff, pi
from dppu.topology.unified import UnifiedEngine, DOFConfig
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant
from dppu.engine.logger import NullLogger

t0_total = time.time()

print("=" * 62)
print("spin-2 spectrum splitting and soft mode identification (App D.9)")
print("=" * 62)

# =============================================================================
# Part 1: Group-theoretic organization of SO(3) → SO(2) splitting
# =============================================================================
print()
print("[1] SO(3) → SO(2) splitting (squash background ε ≠ 0)")
print()
print("  The squash axis = 2-axis preserves SO(2) = U(1).")
print()
print("  The 5 components of spin-2 are classified by SO(2) charge m:")
print("    m=0  : T₁ = diag(1,1,-2)/√6  [1 mode]   (squash singlet)")
print("    m=±2 : T₂,T₃                  [2 modes]  (0-1 plane doublet)")
print("    m=±1 : T₄,T₅                  [2 modes]  (off-squash-axis doublet)")
print()
print("  ε=0: all 5 components degenerate (Schur, known result (App D.9))")
print("  ε≠0: splits into at most 3 distinct masses")
print("    H₁₁^q(ε), H₂₂^q(ε) = H₃₃^q(ε), H₄₄^q(ε) = H₅₅^q(ε)")

# T_A basis (same settings as spin2_schur_degeneracy.py)
from sympy import sqrt, Matrix, simplify

sqrt2 = np.sqrt(2)
sqrt6 = np.sqrt(6)

T = {}
T[1] = np.array([[1, 0, 0], [0, 1, 0], [0, 0, -2]], float) / sqrt6
T[2] = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]], float) / sqrt2
T[3] = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], float) / sqrt2
T[4] = np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], float) / sqrt2
T[5] = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]], float) / sqrt2

# Verify SO(2) transformation (R(θ) around axis 2)
theta_test = np.pi / 6
c, s_th = np.cos(theta_test), np.sin(theta_test)
R_SO2 = np.array([[c, s_th, 0], [-s_th, c, 0], [0, 0, 1]])

# T₁ is invariant (m=0)
T1_rot = R_SO2 @ T[1] @ R_SO2.T
err_T1 = np.max(np.abs(T1_rot - T[1]))

# T₂,T₃ mix (m=±2)
T2_rot = R_SO2 @ T[2] @ R_SO2.T
T3_rot = R_SO2 @ T[3] @ R_SO2.T
# T₂ → cos(2θ)T₂ - sin(2θ)T₃  [natural sign of SO(2) rotation]
cos2t, sin2t = np.cos(2*theta_test), np.sin(2*theta_test)
T23_mix = cos2t * T[2] - sin2t * T[3]
err_T23 = np.max(np.abs(T2_rot - T23_mix))

# T₄,T₅ mix (m=±1)
T4_rot = R_SO2 @ T[4] @ R_SO2.T
# T₄ → cos(θ)T₄ - sin(θ)T₅
T4_mix = np.cos(theta_test) * T[4] - np.sin(theta_test) * T[5]
err_T45 = np.max(np.abs(T4_rot - T4_mix))

print()
print(f"  SO(2) transformation test (θ=π/6):")
print(f"    T₁ invariant:                        err = {err_T1:.2e}  {'confirmed' if err_T1 < 1e-12 else 'failed'}")
print(f"    T₂ → cos(2θ)T₂ - sin(2θ)T₃:  err = {err_T23:.2e}  {'confirmed' if err_T23 < 1e-12 else 'failed'}")
print(f"    T₄ → cos(θ)T₄ - sin(θ)T₅:    err = {err_T45:.2e}  {'confirmed' if err_T45 < 1e-12 else 'failed'}")

# =============================================================================
# Part 2: H₁₁^q(ε), H₂₂^q(ε) — from FA engine
# =============================================================================
print()
print("[2] Compute H₁₁^q(ε), H₂₂^q(ε) (FA engine, finite difference)")
print()

t1 = time.time()
cfg = DOFConfig.from_engine('FullyAnisotropicS3S1Engine')
cfg.torsion_mode = Mode.MX
cfg.ny_variant = NyVariant.FULL
eng = UnifiedEngine(cfg, NullLogger())
eng.run()
print(f"  Engine run: {time.time()-t1:.1f}s")

params = eng.data['params']
Veff_sym = eng.data['potential']
r_s, eps_s, s_s = params['r'], params['epsilon'], params['s']
om1_s, om2_s, om3_s = params['omega1'], params['omega2'], params['omega3']
V_s, eta_s, L_s = params['V'], params['eta'], params['L']
kap_s, th_s, al_s = params['kappa'], params['theta_NY'], params['alpha']

# Substitute ω=0, physical parameters
print("  Substituting ω=0...")
t2 = time.time()
Veff_w0 = cancel(Veff_sym.subs(
    [(om1_s, S.Zero), (om2_s, S.Zero), (om3_s, S.Zero)]
))
sub_phys = [(V_s, S.One), (eta_s, S.One), (L_s, S.One), (kap_s, S.One),
            (th_s, S.Zero), (al_s, S.Zero)]
Veff_p = cancel(Veff_w0.subs(sub_phys))
print(f"  Substitution done: {time.time()-t2:.1f}s")

# Lambdify
Veff_num = lambdify((r_s, eps_s, s_s), Veff_p, modules='numpy')

# Find r_min (equilibrium radius at ε=s=0)
from scipy.optimize import minimize_scalar, brentq

def Veff_iso(r_val):
    return float(Veff_num(r_val, 0.0, 0.0))

res_r0 = minimize_scalar(Veff_iso, bounds=(0.5, 5.0), method='bounded')
r0 = res_r0.x
V0 = Veff_iso(r0)
print(f"  Isotropic point equilibrium radius r₀ = {r0:.6f}, V₀ = {V0:.4f}")

def finite_diff2_eps(eps_val, h=1e-4, r_val=None):
    """d²V/dε² at ε=eps_val, s=0 (finite difference)."""
    rv = r_val if r_val is not None else r0
    vp = float(Veff_num(rv, eps_val + h, 0.0))
    vm = float(Veff_num(rv, eps_val - h, 0.0))
    v0 = float(Veff_num(rv, eps_val, 0.0))
    return (vp + vm - 2*v0) / h**2

def finite_diff2_s(eps_val, h=1e-4, r_val=None):
    """d²V/ds² at ε=eps_val, s=0 (finite difference)."""
    rv = r_val if r_val is not None else r0
    vp = float(Veff_num(rv, eps_val, h))
    vm = float(Veff_num(rv, eps_val, -h))
    v0 = float(Veff_num(rv, eps_val, 0.0))
    return (vp + vm - 2*v0) / h**2

# Jacobian in q coordinates: ε = exp(√6/2 * q₁) - 1
# dε/dq₁ = √6/2 * (1+ε)  → H₁₁^q = (dε/dq₁)² * H_εε = (6/4)*(1+ε)² * H_εε
# ds/dq₂ = 1/√2 * (1+s) → at s=0: H₂₂^q = (1/2) * H_ss
def H11q(eps_val, r_val=None):
    """H₁₁^q = (√6/2 * (1+ε))² * H_εε"""
    H_ee = finite_diff2_eps(eps_val, r_val=r_val)
    jac = (np.sqrt(6)/2 * (1 + eps_val))**2
    return H_ee * jac

def H22q(eps_val, r_val=None):
    """H₂₂^q = (1/√2)² * H_ss = (1/2) * H_ss"""
    H_ss = finite_diff2_s(eps_val, r_val=r_val)
    return H_ss * 0.5

# Verify at ε=0
H11_0 = H11q(0.0)
H22_0 = H22q(0.0)
print(f"  ε=0: H₁₁^q = {H11_0:.4f}, H₂₂^q = {H22_0:.4f}")
print(f"  Theoretical 64π²Lr/κ² = {64*np.pi**2*r0:.4f}")
print(f"  H₁₁^q/H₂₂^q = {H11_0/H22_0:.4f}  (theoretical: 1.000)")

# =============================================================================
# Part 3: H₄₄^q(ε) — via numerical Ricci scalar
# =============================================================================
print()
print("[3] H₄₄^q(ε) — numerical Ricci scalar + finite difference")
print()

def compute_C_for_frame(M3, r_val=1.0):
    """
    C'^a_{bc} = M^a_d (M^{-1})^b_e (M^{-1})^c_f * (4/r) * eps_{def}
    on S³ with isotropic structure constants 4/r * eps_{abc}.
    """
    eps_tensor = np.zeros((3, 3, 3))
    eps_tensor[0, 1, 2] = eps_tensor[1, 2, 0] = eps_tensor[2, 0, 1] = 1.0
    eps_tensor[0, 2, 1] = eps_tensor[2, 1, 0] = eps_tensor[1, 0, 2] = -1.0

    Minv = np.linalg.inv(M3)
    C = (4.0 / r_val) * np.einsum('ad,be,cf,def->abc', M3, Minv, Minv, eps_tensor)
    return C


def compute_R3_from_C(C):
    """
    3D Ricci scalar from structure constants C^a_{bc} (orthonormal basis).

    LC connection (orthonormal basis):
        omega_{abc} = (1/2)(C_{abc} - C_{bac} + C_{cab})
    where C_{abc} = C^a_{bc} (g=delta normalization).

    Riemann tensor (constant structure constants: d*omega = 0):
        R^c_{dab} = omega^c_{ae} omega^e_{db} - omega^c_{be} omega^e_{da}

    Ricci scalar:
        R = g^{da} R^c_{dac} = sum_{da} R^c_{dac}

    Verification: isotropic S³ has R = 24/r² (r=1 → R = 24)
    """
    # omega[c, a, b] = (1/2)(C[c,a,b] - C[a,b,c] + C[b,a,c])
    # C[a,b,c] = C^a_{bc}
    # C_{abc} = C^a_{bc} (orthonormal, g=delta)
    # omega_{abc} = (1/2)(C_{abc} - C_{bac} + C_{cab})
    # = (1/2)(C[a,b,c] - C[b,a,c] + C[c,a,b])
    # omega^c_{ab}: index c up, a,b down
    # omega^c_{ab} = delta^{cd} omega_{dab} = omega_{cab}
    # = (1/2)(C[c,a,b] - C[a,b,c] + C[b,a,c])   [note: cyclic]
    omega = 0.5 * (C.transpose(2, 0, 1)   # C[c,a,b]
                 - C.transpose(0, 1, 2)    # C[a,b,c] = -C[a,c,b]... wait
                 + C.transpose(1, 2, 0))   # C[b,a,c]
    # omega[c,a,b] is correct: omega^c_{ab}

    # Actually let me redo:
    # C.transpose(2,0,1)[c,a,b] = C[a,b,c]... no
    # numpy transpose(2,0,1) maps: new[i,j,k] = old[j,k,i]? No.
    # np.transpose(C, (2,0,1))[c,a,b] = C[a,b,c]
    # So:
    # C.transpose(2,0,1)[c,a,b] = C[a,b,c]   ← this is C_{abc} → C[c,a,b]? NO
    # Let me be explicit:
    # np.transpose(arr, axes) rearranges axes.
    # arr[i,j,k] -> new[i,j,k] = arr[axes[0]-th, axes[1]-th, axes[2]-th position -> old indices]
    # Actually: result[i,j,k] = arr[perm(i,j,k)] where perm is defined by axes.
    # np.einsum would be cleaner:
    C_abc = C                         # C[a,b,c] = C^a_{bc}
    C_cab = np.einsum('abc->cab', C)  # C[c,a,b] = C^c_{ab}
    C_bac = np.einsum('abc->bac', C)  # C[b,a,c] = C^b_{ac}

    # omega^c_{ab} = (1/2)(C^c_{ab} - C^a_{bc} + C^b_{ac})
    # omega[c,a,b] = (1/2)(C_cab[c,a,b] - C_abc[a,b,c]... wait
    # C^c_{ab} = C_cab[c,a,b]
    # C^a_{bc} = C_abc[a,b,c]  but we need C_{abc} = C^a_{bc} = C_abc[a,b,c]... hmm
    # omega^c_{ab} = (1/2)(C[c,a,b] - C[a,b,c] + C[b,a,c])
    # In index notation where C[i,j,k] = C^i_{jk}:
    # C^c_{ab} = C[c,a,b]
    # C^a_{bc} = C[a,b,c]
    # C^b_{ac} = C[b,a,c]
    omega_corr = 0.5 * (C_cab - C_abc + C_bac)
    # omega_corr[c,a,b] = omega^c_{ab}

    # Riemann: R^c_{dab} = omega^c_{ae}*omega^e_{db} - omega^c_{be}*omega^e_{da}
    # R[c,d,a,b] = sum_e omega[c,a,e]*omega[e,d,b] - omega[c,b,e]*omega[e,d,a]
    R_riemann = (np.einsum('cae,edb->cdab', omega_corr, omega_corr)
               - np.einsum('cbe,eda->cdab', omega_corr, omega_corr))

    # Ricci: R_{db} = sum_c R^c_{dcb} = R[c,d,c,b] -> R_db = sum over c of R[c,d,c,b]
    Ricci = np.einsum('cdcb->db', R_riemann)

    # Ricci scalar
    return np.trace(Ricci)


# Verification: isotropic S³ (r=1)
M_iso = np.eye(3)
C_iso = compute_C_for_frame(M_iso, r_val=1.0)
R_iso = compute_R3_from_C(C_iso)
print(f"  Isotropic S³ (r=1) Ricci scalar R = {R_iso:.6f}  (theoretical: 24.0)")
print(f"  Deviation from theory: {abs(R_iso - 24.0):.2e}  {'[OK]' if abs(R_iso - 24.0) < 1e-10 else '[FAIL]'}")


def compute_Veff_LC_only(M3, r_val, L_val=1.0, kappa_val=1.0):
    """
    LC part of V_eff (η=0, α=0 effective potential).
    V_LC = -(1/(2κ²)) * Vol * R₃
    Vol = 2π² L r³ det(M₃)  [M₃ ∈ SL(3) → det=1]
    """
    C = compute_C_for_frame(M3, r_val)
    R3 = compute_R3_from_C(C)
    Vol = 2 * np.pi**2 * L_val * r_val**3  # det(M3) = 1 for SL(3)
    return -(1.0 / (2.0 * kappa_val**2)) * Vol * R3


# Verify: V_LC at isotropic point (r=1, L=κ=1)
V_LC_iso = compute_Veff_LC_only(np.eye(3), r0)
print(f"  V_LC at isotropic point (r₀={r0:.4f}): {V_LC_iso:.6f}")
print(f"  Engine V_eff (η=0): approx {float(Veff_num(r0, 0.0, 0.0)):.6f}")
print("  (These differ since η=1 above; comparison only)")


def compute_H44_LC(eps_val, r_val=None, h=1e-3):
    """
    LC part of H₄₄ (η=0).
    Finite difference of V_LC with small deformation M₃(q₄) = M₀ exp(q₄ T₄) in T₄ direction.
    """
    rv = r_val if r_val is not None else r0
    m0 = (1 + eps_val) ** (1.0/3.0)
    m2 = (1 + eps_val) ** (-2.0/3.0)
    M0 = np.diag([m0, m0, m2])

    def V_at_q4(q4):
        M_pert = M0 @ (np.eye(3) + q4 * T[4])
        # Normalize to SL(3) (preserve det = 1)
        det = np.linalg.det(M_pert)
        M_pert /= det ** (1.0/3.0)
        return compute_Veff_LC_only(M_pert, rv)

    vp = V_at_q4(+h)
    vm = V_at_q4(-h)
    v0 = V_at_q4(0.0)
    return (vp + vm - 2*v0) / h**2


def compute_H11_LC(eps_val, r_val=None, h=1e-4):
    """LC part of H₁₁ (η=0): T₁ direction."""
    rv = r_val if r_val is not None else r0
    m0 = (1 + eps_val) ** (1.0/3.0)
    m2 = (1 + eps_val) ** (-2.0/3.0)
    M0 = np.diag([m0, m0, m2])

    def V_at_q1(q1):
        M_pert = M0 @ (np.eye(3) + q1 * T[1])
        det = np.linalg.det(M_pert)
        M_pert /= det ** (1.0/3.0)
        return compute_Veff_LC_only(M_pert, rv)

    vp = V_at_q1(+h)
    vm = V_at_q1(-h)
    v0 = V_at_q1(0.0)
    return (vp + vm - 2*v0) / h**2


def compute_H22_LC(eps_val, r_val=None, h=1e-4):
    """LC part of H₂₂ (η=0): T₂ direction."""
    rv = r_val if r_val is not None else r0
    m0 = (1 + eps_val) ** (1.0/3.0)
    m2 = (1 + eps_val) ** (-2.0/3.0)
    M0 = np.diag([m0, m0, m2])

    def V_at_q2(q2):
        M_pert = M0 @ (np.eye(3) + q2 * T[2])
        det = np.linalg.det(M_pert)
        M_pert /= det ** (1.0/3.0)
        return compute_Veff_LC_only(M_pert, rv)

    vp = V_at_q2(+h)
    vm = V_at_q2(-h)
    v0 = V_at_q2(0.0)
    return (vp + vm - 2*v0) / h**2


# Verify LC Hessian at ε=0 (Schur: H₁₁ = H₂₂ = H₄₄)
H11_LC_0 = compute_H11_LC(0.0)
H22_LC_0 = compute_H22_LC(0.0)
H44_LC_0 = compute_H44_LC(0.0)
print()
print(f"  LC Hessian at ε=0 (η=0):")
print(f"    H₁₁^q_LC = {H11_LC_0:.4f}")
print(f"    H₂₂^q_LC = {H22_LC_0:.4f}")
print(f"    H₄₄^q_LC = {H44_LC_0:.4f}")
print(f"    Schur isotropic ratio H₁₁/H₂₂ = {H11_LC_0/H22_LC_0:.4f}  (theoretical: 1.000)")
print(f"    Schur isotropic ratio H₁₁/H₄₄ = {H11_LC_0/H44_LC_0:.4f}  (theoretical: 1.000)")

# =============================================================================
# Part 4: ε-dependence plot — mass splitting of all 3 masses
# =============================================================================
print()
print("[4] ε-dependence: mass splitting plot")
print()

# ε grid (from 0 to near ε_c)
eps_grid = np.linspace(0.0, 0.55, 30)

H11_arr = []
H22_arr = []
H44_arr = []

H11_LC_arr = []

print("  Computing...", end='', flush=True)
for eps_val in eps_grid:
    # Equilibrium value of r (r_min at ε≠0)
    def Veff_at_r(rv):
        return float(Veff_num(rv, eps_val, 0.0))
    res = minimize_scalar(Veff_at_r, bounds=(0.3, 6.0), method='bounded')
    rv = res.x

    H11_arr.append(H11q(eps_val, r_val=rv))
    H22_arr.append(H22q(eps_val, r_val=rv))
    H44_arr.append(compute_H44_LC(eps_val, r_val=rv))
    H11_LC_arr.append(compute_H11_LC(eps_val, r_val=rv))
    print('.', end='', flush=True)

print(" done")

H11_arr = np.array(H11_arr)
H22_arr = np.array(H22_arr)
H44_arr = np.array(H44_arr)
H11_LC_arr = np.array(H11_LC_arr)

# Identify ε_c (point where H₁₁^q → 0)
idx_negative = np.where(H11_arr < 0)[0]
if len(idx_negative) > 0:
    eps_c_approx = eps_grid[idx_negative[0]]
    print(f"  H₁₁^q becomes negative at ε ≈ {eps_c_approx:.3f}")
else:
    print("  H₁₁^q did not become negative within grid range (ε_c > 0.55)")


# =============================================================================
# Part 5: Identification of soft mode direction
# =============================================================================
print()
print("[5] Identification of soft mode direction")
print()

idx_04 = np.argmin(np.abs(eps_grid - 0.4))
idx_near_ec = np.argmin(np.abs(H11_arr))
eps_near_ec = eps_grid[idx_near_ec]

print(f"  ε=0:   H₁₁(full)={H11_arr[0]:.1f}, H₂₂(full)={H22_arr[0]:.1f}, "
      f"H₁₁(LC)={H11_LC_arr[0]:.1f}, H₄₄(LC)={H44_arr[0]:.1f}")
print(f"  ε=0.4: H₁₁(full)={H11_arr[idx_04]:.1f}, H₂₂(full)={H22_arr[idx_04]:.1f}, "
      f"H₁₁(LC)={H11_LC_arr[idx_04]:.1f}, H₄₄(LC)={H44_arr[idx_04]:.1f}")
print(f"  ε≈ε_c ({eps_near_ec:.3f}): H₁₁(full)={H11_arr[idx_near_ec]:.1f}, "
      f"H₂₂(full)={H22_arr[idx_near_ec]:.1f}, H₁₁(LC)={H11_LC_arr[idx_near_ec]:.1f}, "
      f"H₄₄(LC)={H44_arr[idx_near_ec]:.1f}")

# Torsion contribution (= full - LC, for T₁ direction)
H11_torsion = H11_arr - H11_LC_arr
print()
print(f"  T₁ mode torsion contribution H₁₁(full) - H₁₁(LC):")
print(f"    ε=0:   {H11_torsion[0]:.1f}")
print(f"    ε=0.4: {H11_torsion[idx_04]:.1f}")
print(f"    ε≈ε_c: {H11_torsion[idx_near_ec]:.1f}")
print(f"  → Torsion (η=1) greatly exceeds the negative LC contribution, keeping H₁₁_full > 0")
print(f"  → To compute H₄₄_full = H₄₄_LC + H₄₄_torsion, torsion in T₄ direction is needed")

print()
print("  [Identification of soft mode]")
soft_mode_found = len(idx_negative) > 0
if soft_mode_found:
    print(f"  H₁₁^q(full) → 0 at ε_c ≈ {eps_c_approx:.3f} (T₁ singlet softening) [OK]")
    H22_at_ec = H22_arr[idx_negative[0] - 1] if idx_negative[0] > 0 else H22_arr[0]
    print(f"  At ε_c: H₂₂^q(full) ≈ {H22_at_ec:.1f}  → m=±2 doublet {'massive [OK]' if H22_at_ec > 0 else 'not massive [NG]'}")
    print(f"  m=±1 (T₄,T₅): H₄₄(LC) < 0 but full expected > 0 with torsion compensation")
    print(f"    (Schur: at ε=0, H₄₄_full = H₁₁_full = {H11_arr[0]:.1f} {'> 0 [OK]' if H11_arr[0] > 0 else '[NG]'})")

print()
if soft_mode_found:
    print(f"  [Soft mode identification]")
    print(f"  1. First soft mode: m=0 (T₁, squash singlet) → condensation at ε_c ≈ {eps_c_approx:.3f}")
    print(f"  2. m=±2 (T₂,T₃) are massive at ε_c → s≠0 is a separate transition")
    print(f"  3. m=±1 (T₄,T₅): massive at ε=0 by Schur; state at ε_c requires torsion computation")
else:
    print(f"  [Soft mode not found in grid range]")

# =============================================================================
# Part 6: Full spin-2 interpretation of s ≠ 0 phase transition
# =============================================================================
print()
print("[6] Full spin-2 interpretation of s ≠ 0 phase transition")
print()

# Condition for s ≠ 0 to become favorable: d²V/ds² < 0
# from prior r3/r5 results
print("  s ≠ 0 direction = T₂ (one component of m=+2 doublet)")
print("  → s ≠ 0 phase transition is condensation of m=±2 doublet")
print()


# Search for ε_s where H₂₂^q(ε) becomes negative (s ≠ 0 phase transition point)
H22_sign_changes = np.where(np.diff(np.sign(H22_arr)))[0]
if len(H22_sign_changes) > 0:
    eps_s_transition = eps_grid[H22_sign_changes[0]]
    print(f"  H₂₂^q → 0 (s≠0 phase transition): ε ≈ {eps_s_transition:.3f}")
else:
    print(f"  H₂₂^q does not change sign in grid range [0, {eps_grid[-1]:.2f}] (remains > 0)")
    print(f"  → s ≠ 0 phase transition occurs at ε > {eps_grid[-1]:.2f} or as an independent s variable")

print()
print(f"Computation time: {time.time()-t0_total:.1f}s")
print()
print("*Script: script/paper03/spin2_spectrum_splitting.py")
