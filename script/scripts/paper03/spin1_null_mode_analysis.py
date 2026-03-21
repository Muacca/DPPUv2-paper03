"""
Null mode interpretation and physical subspace analysis (App E.6)
==========================================================

Objective:
  (1) Identify the null eigenvector direction of G
        - Extract null eigenvector and express in (omega, delta) basis
        - Compute delta_e^a for null direction, determine gauge / auxiliary / redundancy
  (2) Project onto physical subspace
        - Project H onto positive eigenvalue directions of G
        - Solve projected generalized eigenvalue problem to find lambda_phys
  (3) Dimensional interpretation of lambda_phys ≈ 4.19
        - Compare with known mass^2 values (mu^2_omega, mu^2_delta, spin-2 isotropic point value)
        - Verify basis rescaling invariance

Input data (precomputed):
  - G matrix: velocity-mode engine result (spin1_kinetic_metric.py)
  - H matrix: Slice Engine numerical Hessian result (spin1_kinetic_metric.py)
  - r0 = 3.0, L = 1.0

Author: Muacca
Date: 2026-03-15
"""

import sys, os
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import time
import numpy as np

_t_start = time.time()

# ============================================================
# § 0: Load precomputed values
# ============================================================

print("=" * 65)
print("Null mode interpretation / physical subspace analysis (App E.6)")
print("=" * 65)

# --- Physical parameters (isotropic point) ---
r0 = 3.0
L  = 1.0

# --- G (6×6) field-space metric ---
# Source: spin1_kinetic_metric.py (velocity-mode engine)
# Schur structure confirmed: each block is proportional to delta_kl
# Values below are block-diagonal (same 2×2 structure for each axis k=0,1,2)
G_oo_val = 236.870435   # G^{omega omega}  (same for each axis)
G_dd_val = 532.958162   # G^{delta delta}  (same for each axis)
G_od_val = -355.305547  # G^{omega delta}  (same for each axis)

# Build 6×6 G explicitly (2×2 block for each axis k=0,1,2)
# Basis order: (v_omega1, v_omega2, v_omega3, v_delta0, v_delta1, v_delta2)
G_6x6 = np.zeros((6, 6))
for k in range(3):
    G_6x6[k,   k  ] = G_oo_val
    G_6x6[k+3, k+3] = G_dd_val
    G_6x6[k,   k+3] = G_od_val
    G_6x6[k+3, k  ] = G_od_val

# --- H (6×6) Hessian (spin-1 sector, Slice Engine reproduced values) ---
# Source: spin1_kinetic_metric.py (Slice Engine verification)
mu2_oo = 52.6379    # mu^2_omega
mu2_dd = 3789.9295  # mu^2_delta
m_od   = -631.6544  # omega-delta bilinear mixing coefficient

# 6×6 H (Schur structure → 2×2 block per axis)
H_6x6 = np.zeros((6, 6))
for k in range(3):
    H_6x6[k,   k  ] = mu2_oo
    H_6x6[k+3, k+3] = mu2_dd
    H_6x6[k,   k+3] = m_od
    H_6x6[k+3, k  ] = m_od

# ============================================================
# § 1: Eigendecomposition of G / null mode extraction
# ============================================================

print("\n" + "=" * 65)
print("§1. Eigendecomposition of G and null mode extraction")
print("=" * 65)

G_eigs, G_evecs = np.linalg.eigh(G_6x6)  # ascending order

print("\n[Eigenvalues of G (6×6)]")
for i, lam in enumerate(G_eigs):
    print(f"  g_{i} = {lam:12.4f}")

# Null mode determination (numerical threshold)
g_tol = 1e-6 * max(abs(G_eigs))
null_indices = [i for i, g in enumerate(G_eigs) if abs(g) < g_tol]
phys_indices = [i for i, g in enumerate(G_eigs) if abs(g) >= g_tol]

print(f"\n  Zero eigenvalue indices (|g| < {g_tol:.2e}): {null_indices}")
print(f"  Physical eigenvalue indices: {phys_indices}")

# ---- Display null eigenvectors ----
labels = ['omega1', 'omega2', 'omega3', 'delta0', 'delta1', 'delta2']
print("\n[Zero eigenvalue directions of G (null modes)]")
null_vecs = []
for idx in null_indices:
    v = G_evecs[:, idx]
    print(f"\n  g_{idx} = {G_eigs[idx]:.2e}, eigenvector:")
    for j, lbl in enumerate(labels):
        if abs(v[j]) > 1e-10:
            print(f"    {lbl:>10s}: {v[j]:+.6f}")
    null_vecs.append(v)

# ---- Canonical per-axis null basis verification ----
print("\n[Canonical null basis (per-axis): n_k = (v_omega_k=1, v_delta_k=2L/r0) for k=0,1,2]")
print("  (eigh returns an arbitrary orthogonal basis within the null space; the following is the physically natural basis)")
ratio = 2.0 * L / r0
for k in range(3):
    n_can = np.zeros(6)
    n_can[k]   = 1.0       # v_omega_{k}
    n_can[k+3] = ratio     # v_delta_{k} = (2L/r0) v_omega_{k}
    # Numerically verify G n_can = 0
    Gn = G_6x6 @ n_can
    residual = np.linalg.norm(Gn)
    print(f"  n_{k} = (omega{k+1}=1, delta{k}={ratio:.4f}): |G n_{k}| = {residual:.2e}")
print(f"  (For each k: C^3_{{3k}} = 2/r0 - {ratio:.4f}/L = 0 holds)")


# ============================================================
# § 1b: Physical content of null mode (compute delta_e^a)
# ============================================================

print("\n" + "=" * 65)
print("§1b. Physical content of null mode (compute delta e^a)")
print("=" * 65)

# Null direction for each axis k is (v_omega_k, v_delta_k) ∝ (3, 2) ≡ (1, 2L/r0)
# Verify: 2*v_omega/r0 - v_delta/L = 2*(1)/3 - (2/3)/1 = 2/3 - 2/3 = 0 [OK]
null_vw = 1.0
null_vd = 2.0 * L / r0

print(f"\n  r0 = {r0}, L = {L}")
print(f"  Null direction: (v_omega_k, v_delta_k) = ({null_vw:.4f}, {null_vd:.4f})")
print(f"  C^3_{{3k}} velocity term = 2*v_omega/r0 - v_delta/L"
      f" = {2*null_vw/r0:.6f} - {null_vd/L:.6f} = {2*null_vw/r0 - null_vd/L:.2e} [OK]")

print("""
  [Analytic computation of delta e^a (for axis k=0)]
  Coframe at isotropic point (omega=delta=epsilon=s=0, r=r0, L=L):
    e^0 = (r0/2) sigma^0
    e^1 = (r0/2) sigma^1
    e^2 = (r0/2) sigma^2
    e^3 = L dtau

  Variation along null direction: (delta_omega1, delta_delta0) = (1, 2L/r0) = (1, 2/3)""")

# Effect of twist delta_omega1 = 1:
#   e^3 → L(dtau + omega1*sigma^0) → delta e^3 = L * sigma^0 = L * (2/r0) e^0
delta_e3_from_twist = L * (2.0/r0)  # coefficient (in e^0 direction)
print(f"\n  Change from twist (delta_omega1 = 1):")
print(f"    delta e^3 = L * sigma^0 = (2L/r0) e^0 = {delta_e3_from_twist:+.4f} e^0")
print(f"    delta e^0 = 0")

# Effect of mixing delta_delta0 = 2L/r0 (first-order expansion at isotropic point):
#   R_{03}(theta_0) at theta_0 ~ delta_0:
#     e^0 → e^0 + delta_0 * e^3  (e^3 component added to e^0)
#     e^3 → e^3 - delta_0 * e^0  (e^0 component added to e^3)
delta_delta0 = null_vd  # = 2L/r0
delta_e3_from_mixing = -delta_delta0  # coefficient (in e^0 direction)
delta_e0_from_mixing = +delta_delta0  # coefficient (in e^3 direction)
print(f"\n  Change from mixing (delta_delta0 = {delta_delta0:.4f}):")
print(f"    delta e^0 = +delta_delta0 * e^3 = {delta_e0_from_mixing:+.4f} e^3")
print(f"    delta e^3 = -delta_delta0 * e^0 = {delta_e3_from_mixing:+.4f} e^0")

# Total
total_delta_e3 = delta_e3_from_twist + delta_e3_from_mixing
total_delta_e0 = delta_e0_from_mixing
print(f"\n  [Total delta e^a (null direction)]")
print(f"    delta e^0 = {total_delta_e0:+.4f} e^3 = {total_delta_e0*L:+.4f} dtau component")
print(f"    delta e^3 = {total_delta_e3:+.6f} e^0  (= twist + mixing = {delta_e3_from_twist:+.4f} + {delta_e3_from_mixing:+.4f})")
print(f"    delta e^1 = delta e^2 = 0")

# Compute metric change
# delta g_munu = (delta e^a_mu) e^a_nu + e^a_mu (delta e^a_nu)
# delta e^0 = (2L/r0) dtau component → contributes to delta g_{tau, sigma^0 direction}
delta_g_tau_s0 = total_delta_e0 * L + (L * (2.0/r0)) * 0.0  # e^3_tau = L, e^0_{sigma^0} = r0/2
# More precisely:
# (delta e^0)_tau = delta_delta0 * L  (e^3 = L dtau, so (e^3)_tau = L)
# (e^0)_{sigma^0} = r0/2
# delta g_{tau, sigma^0} = (delta e^0)_tau * (e^0)_{sigma^0} + 0
#                        = delta_delta0 * L * (r0/2)
delta_g_cross = delta_delta0 * L * (r0/2.0)

print(f"""
  [Physical metric change delta g_{{mu nu}}]
    delta g_{{tau, sigma^0}} = (delta e^0)_tau * (e^0)_{{sigma^0}}
                             = {delta_delta0:.4f} * {L:.1f} * {r0/2:.4f}
                             = {delta_g_cross:.4f}   (non-zero!)
    delta g_{{sigma^0, tau}} = {delta_g_cross:.4f}  (symmetric)
    Others: = 0

  delta e^3 = 0 → g_{{{3},{3}}} component does not change

  [Determination]
    delta g_{{mu nu}} != 0  →  null mode is not a pure gauge transformation
    However delta e^3 = 0  →  cannot be called gauge in terms of torsion/twist
    → Case B candidate: non-dynamical auxiliary (non-propagating direction with zero kinetic term)""")


# ============================================================
# § 2: Projection onto physical subspace
# ============================================================

print("=" * 65)
print("§2. Projection onto physical subspace and lambda_phys")
print("=" * 65)

# Extract positive eigenvalue directions of G
U_plus = G_evecs[:, phys_indices]   # shape: (6, 3)
G_plus = np.diag(G_eigs[phys_indices])

print(f"\n  Physical subspace dimension: {len(phys_indices)}")
print(f"  G eigenvalues of physical directions: {G_eigs[phys_indices]}")

print(f"\n[Eigenvectors of physical directions (positive eigenvalue directions of G)]")
for ii, idx in enumerate(phys_indices):
    v = G_evecs[:, idx]
    print(f"  u_{ii} (g={G_eigs[idx]:.2f}): ", end="")
    parts = [f"{v[j]:+.4f}*{labels[j]}" for j in range(6) if abs(v[j]) > 1e-10]
    print(", ".join(parts))

# Project H onto physical subspace
H_phys = U_plus.T @ H_6x6 @ U_plus   # (3, 3)
G_phys = U_plus.T @ G_6x6 @ U_plus   # (3, 3), should be diag(g+)

print(f"\n[H_phys (3×3, projected Hessian)]")
np.set_printoptions(precision=4, suppress=True)
print(H_phys)

print(f"\n[G_phys (3×3, projected metric, should be diagonal)]")
print(G_phys)

# Generalized eigenvalue problem: H_phys v = lambda G_phys v
# If G_phys is positive definite, can be solved with scipy.linalg.eigh
from scipy.linalg import eigh as scipy_eigh

G_phys_eigs = np.linalg.eigvalsh(G_phys)
print(f"\n  G_phys eigenvalues: {G_phys_eigs}")
G_phys_posdef = all(e > 0 for e in G_phys_eigs)
print(f"  G_phys positive definite: {'YES' if G_phys_posdef else 'NO'}")

if G_phys_posdef:
    lam_phys, evec_phys = scipy_eigh(H_phys, G_phys)
    print(f"\n[Generalized eigenvalues lambda_phys = H_phys v / G_phys v]")
    for i, lam in enumerate(lam_phys):
        print(f"  lambda_{i} = {lam:10.4f}  {'(< 0: unstable candidate)' if lam < 0 else '(> 0: stable direction)'}")
    lam_min = lam_phys[0]
else:
    print("  Cannot solve generalized eigenvalue problem since G_phys is not positive definite.")
    lam_min = float('nan')

# Manual verification: physical direction u = (-2, 3) for 1-axis (k=0)
print("\n[Manual verification: physical direction u = (-2, 3) for 1-axis (k=0)]")
# 2×2 block (omega1, delta0 subspace)
H_2 = np.array([[mu2_oo, m_od], [m_od, mu2_dd]])
G_2 = np.array([[G_oo_val, G_od_val], [G_od_val, G_dd_val]])
u   = np.array([-2.0, 3.0])
u_norm = u / np.linalg.norm(u)

uHu = u @ H_2 @ u
uGu = u @ G_2 @ u
lam_manual = uHu / uGu if abs(uGu) > 1e-12 else float('nan')
print(f"  u = (-2, 3),  u^T H u = {uHu:.4f},  u^T G u = {uGu:.4f}")
print(f"  lambda_phys (manual) = {uHu:.4f} / {uGu:.4f} = {lam_manual:.6f}")

# ============================================================
# § 3: Dimensional interpretation of lambda_phys ≈ 4.19
# ============================================================

print("\n" + "=" * 65)
print("§3. Dimensional and normalization interpretation of lambda_phys")
print("=" * 65)


# Comparison with known mass^2 values
print("  [Known mass^2 values (same conventions)]")
mu2_spin2_iso = 48 * np.pi**2   # spin-2 isotropic point mass (known result (App D.9))
print(f"  spin-2 isotropic point mu^2    = 48 pi^2 ~= {mu2_spin2_iso:.4f}")
print(f"  spin-1 omega  mu^2   = {mu2_oo:.4f}")
print(f"  spin-1 delta  mu^2   = {mu2_dd:.4f}")
print(f"  lambda_phys          = {lam_manual:.4f}")
print()
print(f"  Ratio lambda_phys / mu^2_spin2 = {lam_manual / mu2_spin2_iso:.6f}")
print(f"  Ratio lambda_phys / mu^2_omega = {lam_manual / mu2_oo:.6f}")

# Verify basis rescaling invariance
print("\n  [Basis rescaling invariance]")
print("  For omega_k → c * omega_k, delta_k → c * delta_k:")
print("  H_ij → H_ij / c^2, G_ij → G_ij / c^2")
print("  → lambda_phys = (u^T H u / c^2) / (u^T G u / c^2) = u^T H u / u^T G u")
print("  → lambda_phys is invariant under isotropic rescaling")

# Anisotropic scaling: omega → a * omega, delta → b * delta
print("\n  Anisotropic scaling: omega → a*omega, delta → b*delta:")
a, b = 2.0, 0.5   # example
H_scaled = H_2.copy()
H_scaled[0, 0] *= 1.0/(a*a)  # d^2/d(a*omega)^2 = (1/a^2) d^2/domega^2
H_scaled[1, 1] *= 1.0/(b*b)
H_scaled[0, 1] *= 1.0/(a*b)
H_scaled[1, 0] *= 1.0/(a*b)
G_scaled = G_2.copy()
G_scaled[0, 0] *= 1.0/(a*a)
G_scaled[1, 1] *= 1.0/(b*b)
G_scaled[0, 1] *= 1.0/(a*b)
G_scaled[1, 0] *= 1.0/(a*b)

if G_scaled[0,0]*G_scaled[1,1] - G_scaled[0,1]**2 > 1e-12:
    lams_scaled, _ = scipy_eigh(H_scaled, G_scaled)
    print(f"  a={a}, b={b}: lambda_phys = {lams_scaled[0]:.4f}, {lams_scaled[1]:.4f}")
    print(f"  → lambda_phys changes under anisotropic scaling (differs from single-variable eigenvalue)")
else:
    # singular G after scaling → use manual calculation
    u_scaled = np.array([-2.0*a, 3.0*b])  # u transformed according to scaling
    lam_sc = (u_scaled @ H_2 @ u_scaled) / (u_scaled @ G_2 @ u_scaled)
    print(f"  a={a}, b={b}: lambda_phys (manual) = {lam_sc:.4f}")
    print(f"  → lambda_phys changes under anisotropic scaling ({lam_manual:.4f} → {lam_sc:.4f})")

# ============================================================
# § 4: Interpretation of bilinear mixing in physical basis
# ============================================================

print("\n" + "=" * 65)
print("§4. Interpretation of omega-delta mixing in physical basis")
print("=" * 65)

# physical direction: u_phys = (-2, 3) / ||(-2, 3)||
# null direction:     u_null = (1, 2/3) / ||(1, 2/3)||
norm_phys = np.sqrt(4 + 9)
norm_null = np.sqrt(1 + (2.0/3)**2)
u_phys_vec = np.array([-2.0, 3.0]) / norm_phys
u_null_vec = np.array([1.0, 2.0/3]) / norm_null

print(f"\n  Physical basis X (propagating direction):  ({-2/norm_phys:+.4f}, {3/norm_phys:+.4f})")
print(f"    omega component = {-2/norm_phys:+.4f}, delta component = {3/norm_phys:+.4f}")
print(f"    |delta/omega| = {3.0/2:.4f}  (large delta component)")
print()
print(f"  Null basis Y (non-propagating direction): ({1/norm_null:+.4f}, {(2.0/3)/norm_null:+.4f})")
print(f"    omega component = {1/norm_null:+.4f}, delta component = {(2.0/3)/norm_null:+.4f}")
print(f"    |delta/omega| = {(2.0/3):.4f}  (omega dominant)")


# ============================================================
# § 5: Summary of determinations
# ============================================================

print("=" * 65)
print("§5. Summary of determinations")
print("=" * 65)

print(f"""
  [1] Null mode of G
  Physical subspace dimension: {len(phys_indices)}
  Null eigenvector: Y_k ∝ (omega_k, delta_k) = (1, 2L/r0) = (1, {null_vd:.4f})
    → C^3_{{3k}} = 2*v_omega_k/r0 - v_delta_k/L = {2*null_vw/r0 - null_vd/L:.2e}  (= 0: kinetic term zero)
    → delta g_{{mu nu}} != 0  (physically distinct configuration)

  [2] Physical subspace
  Physical basis: X_k = (-2 omega_k + 3 delta_k)/sqrt(13)  (k=0,1,2)
  lambda_phys = {lam_manual:.4f}  {'> 0 [OK]' if lam_manual > 0 else '< 0 [NG]'}  (3-fold degenerate)

  [3] Mass^2 comparison
  mu^2_omega = {mu2_oo:.4f},  mu^2_delta = {mu2_dd:.4f},  lambda_phys = {lam_manual:.4f}
""")

print("=" * 65)
print("Done")
print("=" * 65)

print(f"\nComputation time: {time.time() - _t_start:.1f}s")
print()
print("*Script: scripts/paper03/spin1_null_mode_analysis.py")
