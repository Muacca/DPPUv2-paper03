"""
Confirmation of Case C (variable redundancy) exclusion (App E.6)
=========================================================

Objective:
  Confirm that the null mode (delta_omega_k=1, delta_delta_k=2L/r0) is NOT

    Case C: implementation degeneracy (redundancy in variable definition)

  For Case C to hold:
    (omega_k, delta_k) = (0, 0) and (omega_k + eps, delta_k + 2L*eps/r0)
    must give the same physical coframe e^a.

  This script confirms:
    (A) Explicit coframe matrix: null direction gives delta g_{mu nu} != 0 → physically distinct
    (B) Curvature of H in null direction: n^T H n > 0 → V_eff has non-zero curvature → constraint equation exists
    (C) Confirmation of Case B (non-dynamical):
          G_null = 0, H_null > 0 → "infinite mass" constraint direction
    (D) Numerical: verify divergence of lambda_eff = H_null / G_null

Author: Muacca
Date: 2026-03-15
"""

import sys, os
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from dppu.utils.tee_logger import setup_log, teardown_log
setup_log(__file__, log_dir='output')

_t_start = time.time()

# ============================================================
# § 0: Parameters and matrices
# ============================================================

print("=" * 65)
print("Case C (variable redundancy) exclusion confirmation (App E.6)")
print("=" * 65)

r0 = 3.0
L  = 1.0

# G_red (2×2, common per-axis block for each axis)
G_oo = 236.870435
G_dd = 532.958162
G_od = -355.305547
G_red = np.array([[G_oo, G_od],
                  [G_od, G_dd]])

# H_red (2×2, common per-axis block for each axis)
mu2_oo = 52.6379
mu2_dd = 3789.9295
m_od   = -631.6544
H_red = np.array([[mu2_oo, m_od],
                  [m_od,   mu2_dd]])

# null / physical directions
null_vw = 1.0
null_vd = 2.0 * L / r0          # = 2/3
n_vec = np.array([null_vw, null_vd])       # null direction
u_vec = np.array([-2.0,    3.0])           # physical direction (orthogonal to n)

# ============================================================
# Explicit coframe matrix
# ============================================================

print("\n" + "=" * 65)
print("Coframe matrix: null direction != identity transformation")
print("=" * 65)

e_iso = np.array([
    [0,    r0/2, 0,    0   ],   # e^0
    [0,    0,    r0/2, 0   ],   # e^1
    [0,    0,    0,    r0/2],   # e^2
    [L,    0,    0,    0   ],   # e^3
])

print("  Isotropic point coframe matrix e^a_mu:")
mu_labels = ['tau', 'sigma0', 'sigma1', 'sigma2']
print("       " + "  ".join(f"{l:>8s}" for l in mu_labels))
for a, lbl in enumerate(['e^0', 'e^1', 'e^2', 'e^3']):
    row = "  " + f"{lbl}  " + "  ".join(f"{e_iso[a, mu]:8.4f}" for mu in range(4))
    print(row)

# After deformation in null direction (small finite deformation eps = 0.001)
eps = 0.001
d_omega1 = eps * null_vw      # = 0.001
d_delta0 = eps * null_vd      # = 0.001 * 2/3

# first-order twist: delta e^3 += L * d_omega1 * sigma^0  →  e^3[sigma^0] += L*d_omega1
# first-order mixing: delta e^0 += d_delta0 * e^3_iso     →  e^0[tau] += d_delta0 * L
#                     delta e^3 += -d_delta0 * e^0_iso    →  e^3[sigma^0] -= d_delta0 * r0/2
e_null = e_iso.copy().astype(float)
e_null[3, 1] += L * d_omega1             # e^3[sigma^0] += L*d_omega1 (twist)
e_null[0, 0] += d_delta0 * L            # e^0[tau]     += d_delta0 * L (mixing)
e_null[3, 1] -= d_delta0 * (r0 / 2)     # e^3[sigma^0] -= d_delta0 * r0/2 (mixing)

delta_e = e_null - e_iso

print(f"\n  After deformation in null direction (eps={eps}): (delta_omega1={d_omega1:.4f}, delta_delta0={d_delta0:.4f})")
print("  delta e^a_mu (changes):")
print("       " + "  ".join(f"{l:>8s}" for l in mu_labels))
for a, lbl in enumerate(['e^0', 'e^1', 'e^2', 'e^3']):
    row = "  " + f"{lbl}  " + "  ".join(f"{delta_e[a, mu]:8.4f}" for mu in range(4))
    print(row)

# Metric change
g_iso  = e_iso.T  @ e_iso    # g_iso = sum_a e^a_mu e^a_nu
g_null = e_null.T @ e_null
delta_g = g_null - g_iso

print(f"\n  delta g_{{mu nu}} = g_null - g_iso (numerical, eps={eps}):")
print("       " + "  ".join(f"{l:>8s}" for l in mu_labels))
for mu, lbl in enumerate(mu_labels):
    row = "  " + f"{lbl}  " + "  ".join(f"{delta_g[mu, nu]:8.4f}" for nu in range(4))
    print(row)

delta_g_max = np.max(np.abs(delta_g))
print(f"\n  max |delta g| = {delta_g_max:.6f}  {'(non-zero: Case C excluded [OK])' if delta_g_max > 1e-10 else '(zero: Case C holds?)'}")
print(f"  Main non-zero component: g_{{tau, sigma^0}} = {delta_g[0, 1]:.6f}, g_{{sigma^0, tau}} = {delta_g[1, 0]:.6f}")

# Verify against analytic value
# Analytic: delta g_{tau, sigma^0} = delta_delta0 * L * (r0/2) = (2L/r0 * eps) * L * (r0/2)
#       = eps * L^2 = 0.001 (L=1)
analytic_delta_g = eps * L**2
print(f"\n  Comparison with analytic value:")
print(f"    Analytic: delta g_{{tau,sigma^0}} = eps * L^2 = {analytic_delta_g:.6f}")
print(f"    Numerical: delta g_{{tau,sigma^0}} = {delta_g[0, 1]:.6f}")
print(f"    Match: {'[OK]' if abs(delta_g[0, 1] - analytic_delta_g) < 1e-8 else '[NG]'}")

# ============================================================
# Curvature of H in null direction (2nd condition for Case C exclusion)
# ============================================================

print("\n" + "=" * 65)
print("Curvature of H in null direction: compute n^T H n")
print("=" * 65)

G_null = n_vec @ G_red @ n_vec
H_null = n_vec @ H_red @ n_vec

print(f"  Null direction n = ({null_vw}, {null_vd:.4f})")
print(f"  G_null = n^T G n = {G_null:.6f}  {'(= 0: kinetic term zero [OK])' if abs(G_null) < 1e-3 else '(non-zero)'}")
print(f"  H_null = n^T H n = {H_null:.6f}  {'(> 0: V_eff has positive curvature [Case B])' if H_null > 0 else '(= 0: flat → Case C holds?)'}")

print()
if abs(G_null) < 1e-3 and H_null > 0:
    print("  Determination: G_null = 0, H_null > 0")
    print("    → Null direction generates the constraint equation (d/dq_null)V = 0")
    print("    → 'Infinite effective mass' non-propagating direction → Case B (non-dynamical auxiliary) [OK]")
    print("    → Case C (variable redundancy) does NOT hold")
elif abs(G_null) < 1e-3 and abs(H_null) < 1e-3:
    print("  Determination: G_null = 0, H_null = 0")
    print("    → V_eff is also flat → Case C (variable redundancy) possible!")
else:
    print(f"  Determination: unexpected result. G_null={G_null:.6f}, H_null={H_null:.6f}")

# ============================================================
# Euler-Lagrange equation for null direction (constraint equation)
# ============================================================

print("\n" + "=" * 65)
print("Euler-Lagrange equation for null direction (constraint equation)")
print("=" * 65)

print(f"  H_null = {H_null:.4f}  {'> 0: constraint equation forcing q_n = 0' if H_null > 0 else '= 0: degree of freedom (Case C)'}")

# ============================================================
# Numerical divergence of effective eigenvalue lambda_eff = H/G
# ============================================================

print("\n" + "=" * 65)
print("Divergence of effective eigenvalue lambda_eff = (n^T H n) / (n^T G n)")
print("=" * 65)

# Normalize n and u
G_XX = u_vec @ G_red @ u_vec
H_XX = u_vec @ H_red @ u_vec
H_cross = n_vec @ H_red @ u_vec
lambda_proj = H_XX / G_XX

print(f"  G_null  = n^T G n = {G_null:.6f}")
print(f"  H_null  = n^T H n = {H_null:.6f}")
print(f"  G_XX    = u^T G u = {G_XX:.4f}")
print(f"  H_XX    = u^T H u = {H_XX:.4f}")
print(f"  H_cross = n^T H u = {H_cross:.4f}")
print(f"  lambda_proj = H_XX / G_XX = {lambda_proj:.6f}  (pre-Schur projection quotient)")

print("\n  lambda_eff when varying t/s ratio (mixed direction):")
print(f"  {'t/s':>8s}  {'G_v':>14s}  {'H_v':>14s}  {'lambda_eff':>14s}")
for ts in [0.0, 0.1, 1.0, 10.0, 100.0]:
    # G, H for direction v = ts*n/|n| + u/|u|
    # Compute without normalization
    v_test = ts * n_vec + u_vec
    Gv = v_test @ G_red @ v_test
    Hv = v_test @ H_red @ v_test
    lam = Hv / Gv if abs(Gv) > 1e-12 else float('inf')
    print(f"  {ts:>8.1f}  {Gv:>14.4f}  {Hv:>14.4f}  {lam:>14.4f}")

if abs(G_null) < 1e-3 and H_null > 0:
    print("  lambda_eff -> inf as null direction dominates  [Case C excluded: OK]")
else:
    print("  lambda_eff does not diverge as expected  [NG]")
print()

# Condition flags for conclusion
cond_A = delta_g_max > 1e-10
cond_B = H_null > 0
cond_C = abs(G_null) < 1e-3 and H_null > 0
case_c_excluded = cond_A and cond_B and cond_C

# ============================================================
# Critical difference between Case C and Case B
# ============================================================

print("=" * 65)
print("Critical differences between Case C and Case B")
print("=" * 65)

print(f"""
  [Necessary conditions for Case C: variable redundancy]
    1. delta g_{{mu nu}} = 0  (deformation in null direction is physically void)
    2. H_null = 0           (V_eff is also flat)
    3. lambda_eff finite   (eigenvalue is finite)

  [Results of this computation]
    1. max|delta g| = {delta_g_max:.6f} {'!= 0' if cond_A else '= 0 [unexpected]'}   → Case C necessary condition 1 {'not satisfied [OK]' if cond_A else 'SATISFIED [NG]'}
    2. H_null = {H_null:.4f} {'> 0' if cond_B else '<= 0 [unexpected]'}           → Case C necessary condition 2 {'not satisfied [OK]' if cond_B else 'SATISFIED [NG]'}
    3. lambda_eff -> {'inf (t/s >> 1)' if cond_C else 'finite [unexpected]'}       → Case C necessary condition 3 {'not satisfied [OK]' if cond_C else 'SATISFIED [NG]'}

  [Conclusion]
    Case C (variable redundancy) {'is excluded.' if case_c_excluded else 'is NOT excluded. [NG]'}

  [Confirmation of Case B (re-verification)]
    G_null = {G_null:.2e} ~= 0   (kinetic term is zero)
    H_null = {H_null:.4f} {'> 0   (potential generates constraint)' if cond_B else '<= 0   (no constraint) [NG]'}
    delta g_{{mu nu}} {'!= 0' if cond_A else '= 0 [NG]'}       ({'physically distinct configurations' if cond_A else 'physically equivalent — Case C possible'})

    → Null mode is a "{'non-propagating constrained direction (non-dynamical constraint)' if case_c_excluded else '?'}"
    → {'Case B: non-dynamical auxiliary is the most appropriate classification' if case_c_excluded else 'Case B cannot be confirmed [NG]'}

""")

# ============================================================
# Summary
# ============================================================

print("=" * 65)
print("Summary: Step 1 candidate completed")
print("=" * 65)

print(f"""
  [Confirmed items: exclusion of Case C (variable redundancy)]

    (A) delta g_{{mu nu}} != 0:
        max|delta g| = {delta_g_max:.6e}  (non-zero)
        Main component: delta g_{{tau,sigma^0}} = {delta_g[0,1]:.6f}
        Analytic value: eps * L^2 = {eps * L**2:.6f}  {'[match OK]' if abs(delta_g[0, 1] - analytic_delta_g) < 1e-8 else '[match NG]'}
        → Null direction gives physically distinct coframes

    (B) H_null = n^T H n = {H_null:.4f} {'> 0' if cond_B else '<= 0 [NG]'}:
        → V_eff has positive curvature in null direction → constraint equation exists
        → Excludes Case C (requires H_null=0)

    (C) lambda_eff = H/G -> inf (when null direction dominates):
        → Infinite effective mass → non-propagating direction
        → Excludes Case C (requires finite lambda_eff)

  [Final determination]
    Case C (variable redundancy): {'excluded [OK]' if case_c_excluded else 'NOT excluded [NG]'}
    Case B (non-dynamical auxiliary): {'adopted [determined]' if case_c_excluded else 'inconclusive [NG]'}

""")

print("=" * 65)
print("Done: " + ("Case C exclusion confirmed" if case_c_excluded else "Case C exclusion FAILED [NG]"))
print("=" * 65)

print(f"\nComputation time: {time.time() - _t_start:.1f}s")
print()
print("*Script: scripts/paper03/spin1_auxiliary_exclusion.py")
teardown_log()
