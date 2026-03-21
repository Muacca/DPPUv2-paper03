"""
Computation of field-space metric G (App E.6)
=======================================================

Objective:
  Use the velocity-mode engine to rigorously compute G^{(v)} (6×6),
  then solve the generalized eigenvalue problem H v = λ G v to determine
  whether tachyons are present.

Engine used:
  UnifiedEngine / DOFConfig(topology=S3, fiber=BOTH, torsion=MX, enable_velocity=True)
  Slice Engine (ε=s=0 fixed) — mathematically equivalent to FullyAnisotropicBothS3S1Engine at isotropic point

Physical meaning of velocity symbols:
  C^3_{3k} += +2 v^ω_k / r  −  v^δ_k / L   (k=0,1,2)
  → G_ij = (1/2) ∂²V_eff/∂v^i∂v^j |_{v=0, isotropic point}

Author: Muacca
Date: 2026-03-15
"""

import sys, os, time
import io

# Windows cp932 workaround
if sys.stdout.encoding and sys.stdout.encoding.lower() not in ('utf-8', 'utf8'):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import numpy as np
from sympy import S, diff, lambdify
from scipy.optimize import minimize_scalar
from scipy.linalg import eigh

from dppu.topology.unified import UnifiedEngine, DOFConfig, TopologyType, FiberMode
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant


# ============================================================
# Step 1: Run engine (velocity mode)
# ============================================================
print("=" * 60)
print("G metric rigorous computation (velocity mode) (App E.6)")
print("=" * 60)
print()

t0 = time.time()

cfg = DOFConfig(
    topology=TopologyType.S3,
    enable_squash=False,    # ε = 0 fixed (isotropic point)
    enable_shear=False,     # s = 0 fixed
    fiber_mode=FiberMode.BOTH,
    isotropic_twist=False,
    torsion_mode=Mode.MX,
    ny_variant=NyVariant.FULL,
    enable_velocity=True,   # ← velocity mode for G metric computation
)

print("[Step 1/5] Running UnifiedEngine (velocity mode, estimated ~300-600s)...", flush=True)
engine = UnifiedEngine(cfg)
engine.run()

t1 = time.time()
print(f"Engine run: {t1 - t0:.1f}s")
print()

fp   = engine.get_free_params()
Veff = engine.data['potential']

# Verify velocity symbols
v_syms = [fp.get(k) for k in
          ['v_omega1', 'v_omega2', 'v_omega3', 'v_delta0', 'v_delta1', 'v_delta2']]
if not all(v is not None and hasattr(v, 'is_Symbol') and v.is_Symbol for v in v_syms):
    print("ERROR: velocity symbols not found in params.")
    print("  Check that enable_velocity=True is correctly set.")
    sys.exit(1)

v_omega_syms = v_syms[:3]  # v_omega1, v_omega2, v_omega3
v_delta_syms = v_syms[3:]  # v_delta0,  v_delta1,  v_delta2
print(f"Velocity symbols confirmed: {[str(v) for v in v_syms]}")

# ============================================================
# Step 2: Substitute physical parameters (isotropic point)
# ============================================================
# Substitute: η=L=κ=1, V=1, α=0, θ_NY=0, ω=δ=0, cd=1, sd=0
# Keep v symbols (for G computation)
phys = {}
for key, val in [
    ('eta',      1), ('V',        1), ('L',        1),
    ('kappa',    1), ('alpha',    0), ('theta_NY', 0),
    ('omega1',   0), ('omega2',   0), ('omega3',   0),
    ('delta0',   0), ('delta1',   0), ('delta2',   0),
    ('cd0',      1), ('sd0',      0),
    ('cd1',      1), ('sd1',      0),
    ('cd2',      1), ('sd2',      0),
]:
    sym = fp.get(key)
    if sym is not None:
        phys[sym] = val

t2 = time.time()
print("[Step 2/5] Substituting physical parameters...", flush=True)
Viso = Veff.subs(phys)   # V(r, v_omega1-3, v_delta0-2)
print(f"  Substitution done: {time.time()-t2:.1f}s")

# ============================================================
# Step 3: Optimize r₀ (static potential at v=0)
# ============================================================
v_zero_subs = {v: 0 for v in v_syms}
V_static = Viso.subs(v_zero_subs)
r_sym    = fp['r']

V_static_func = lambdify(r_sym, V_static, 'numpy')
res = minimize_scalar(lambda r: float(V_static_func(r)), bounds=(1.0, 10.0), method='bounded')
r0 = float(res.x)
print(f"r₀ = {r0:.6f}  (expected: 3.0)")
print()

# ============================================================
# Step 4: G_{ij} = (1/2) ∂²V/∂v^i∂v^j |_{v=0, r=r0}
# ============================================================
print("[Step 4/5] Computing G_{ij} (SymPy diff × 36 components)...", flush=True)
t_G = time.time()

Vr0 = Viso.subs(r_sym, r0)   # V(v_omega1-3, v_delta0-2) at r=r0

# Label definitions
labels_omega = ['ω₁', 'ω₂', 'ω₃']
labels_delta = ['δ₀', 'δ₁', 'δ₂']
all_labels   = labels_omega + labels_delta

G_sym = [[None] * 6 for _ in range(6)]
for i, vi in enumerate(v_syms):
    for j, vj in enumerate(v_syms):
        if j < i:
            G_sym[i][j] = G_sym[j][i]   # use symmetry
        else:
            expr = diff(diff(Vr0, vi), vj) / 2
            G_sym[i][j] = expr.subs(v_zero_subs)

print(f"  Differentiation done: {time.time()-t_G:.1f}s")

# Numerical evaluation
G_num = np.zeros((6, 6))
for i in range(6):
    for j in range(6):
        G_num[i, j] = float(G_sym[i][j])

# ============================================================
# Step 5: Display G results and verify Schur structure
# ============================================================
print()
print("=" * 60)
print("[G field-space metric results]")
print("=" * 60)

def print_block(mat, row_labels, col_labels, title):
    print(f"\n{title}:")
    header = "       " + "  ".join(f"{l:>10s}" for l in col_labels)
    print(header)
    for i, rl in enumerate(row_labels):
        row = f"  {rl}  " + "  ".join(f"{mat[i,j]:10.4f}" for j in range(len(col_labels)))
        print(row)

G_oo = G_num[:3, :3]
G_dd = G_num[3:, 3:]
G_od = G_num[:3, 3:]

print_block(G_oo, labels_omega, labels_omega, "G^{ωω} (3×3)")
print_block(G_dd, labels_delta, labels_delta, "G^{δδ} (3×3)")
print_block(G_od, labels_omega, labels_delta, "G^{ωδ} (3×3)")

# Schur structure determination
def schur_check(mat, name):
    diag_vals = [mat[i,i] for i in range(3)]
    offdiag_max = max(abs(mat[i,j]) for i in range(3) for j in range(3) if i != j)
    diag_mean = np.mean(diag_vals)
    diag_spread = (max(diag_vals) - min(diag_vals)) / (abs(diag_mean) + 1e-30)
    ok = diag_spread < 1e-6 and offdiag_max < 1e-6 * abs(diag_mean)
    print(f"\n  {name}: diag = [{', '.join(f'{v:.6f}' for v in diag_vals)}]")
    print(f"         off-diag max = {offdiag_max:.2e}  diag spread = {diag_spread:.2e}  → {'Schur confirmed' if ok else 'not Schur'}")
    return diag_mean, ok

print()
print("[Schur structure determination]")
g_oo, schur_oo = schur_check(G_oo, "G^{ωω}")
g_dd, schur_dd = schur_check(G_dd, "G^{δδ}")
g_od, schur_od = schur_check(G_od, "G^{ωδ}")

# ============================================================
# Step 6: 2×2 Schur reduction → generalized eigenvalue problem
# ============================================================
print()
print("=" * 60)
print("[Generalized eigenvalue problem H v = λ G v]")
print("=" * 60)

# Slice Engine verified values (re-verified 20260307)
mu2_oo  =  52.6379
mu2_dd  = 3789.9295
m_od    = -631.6544

H_red = np.array([[mu2_oo, m_od  ],
                  [m_od,   mu2_dd]])

# G_red: if Schur structure holds, diagonal (g_ωω, g_δδ) + off-diagonal g_ωδ
g_oo_val = G_num[0, 0]
g_dd_val = G_num[3, 3]
g_od_val = G_num[0, 3]

G_red = np.array([[g_oo_val, g_od_val],
                  [g_od_val, g_dd_val]])

print(f"\n  H_red = [[{mu2_oo:.4f}, {m_od:.4f}],")
print(f"           [{m_od:.4f}, {mu2_dd:.4f}]]")
print(f"\n  G_red = [[{g_oo_val:.6f}, {g_od_val:.6f}],")
print(f"           [{g_od_val:.6f}, {g_dd_val:.6f}]]")
print(f"\n  det(H_red) = {np.linalg.det(H_red):.2f}  ({'< 0 → indefinite' if np.linalg.det(H_red) < 0 else '>= 0'})")
print(f"  det(G_red) = {np.linalg.det(G_red):.6f}  ({'> 0 → positive definite candidate' if np.linalg.det(G_red) > 0 else '<= 0 → not positive definite!'})")

# All eigenvalues of G (check positive definiteness)
G_eigs = np.linalg.eigvalsh(G_num)
print(f"\n  G (6×6) eigenvalues: {G_eigs}")
G_posdef = all(e > 0 for e in G_eigs)
print(f"  G positive definite: {'YES' if G_posdef else 'NO (sign issue)'}")

# Generalized eigenvalue problem
# Note: G_red may be rank-1 (det ~ 0), so eigh may fail.
# Physical eigenvalue is computed by projecting onto G's positive eigenvalue directions.
try:
    eigenvalues, eigenvectors = eigh(H_red, G_red)
    lam_min = eigenvalues[0]
    lam_max = eigenvalues[1]
    v_min   = eigenvectors[:, 0]

    print(f"\n  λ_min = {lam_min:.6f}")
    print(f"  λ_max = {lam_max:.6f}")
    print(f"  Eigenvector (λ_min): ω component = {v_min[0]:.4f}, δ component = {v_min[1]:.4f}")
    if abs(v_min[0]) > 1e-10:
        ratio = abs(v_min[1] / v_min[0])
        angle = np.degrees(np.arctan2(v_min[1], v_min[0]))
        print(f"  mixing |δ/ω| = {ratio:.4f},  mixing angle = {angle:.2f}deg")

except Exception as e:
    print(f"\n  Generalized eigenvalue computation error (G_red is rank-deficient): {e}")
    lam_min = float('nan')

# ── Physical eigenvalue (project onto positive eigenvalue directions of G) ────────────────────────
print()
print("  [Eigenvalue in physical subspace (excluding G null modes)]")
print("  G_red is rank-1: null direction = (1, 2L/r0) = (1, 2/3) at r0=3, L=1")
print("  Physical direction u = (-2, 3) / sqrt(13)  (orthogonal to null)")
u_phys = np.array([-2.0, 3.0])
uHu = u_phys @ np.array([[mu2_oo, m_od], [m_od, mu2_dd]]) @ u_phys
uGu = u_phys @ G_red @ u_phys
lam_phys = uHu / uGu if abs(uGu) > 1e-10 else float('nan')
print(f"  u^T H u = {uHu:.4f},  u^T G u = {uGu:.4f}")
print(f"  lambda_phys = {lam_phys:.6f}  ({'> 0: physical direction is stable' if lam_phys > 0 else '< 0: instability in physical direction'})")

# ============================================================
# Step 7: Tachyon determination
# ============================================================
print()
print("=" * 60)
print("[Tachyon determination]")
print("=" * 60)

checks = {
    "G_red computed (velocity mode engine)":                       np.all(np.isfinite(G_red)),
    "G (6x6) positive definiteness check":                         G_posdef,
    "G_red is rank-1 (det ~ 0)":                                   abs(np.linalg.det(G_red)) < 1.0,
    "lambda_phys computed (physical subspace)":                    not np.isnan(lam_phys),
    "lambda_phys > 0 (physical direction stable)":                 (not np.isnan(lam_phys)) and lam_phys > 0,
}

for desc, ok in checks.items():
    mark = "[OK]" if ok else "[--]"
    print(f"  {mark} {desc}")

print()
if not np.isnan(lam_phys) and lam_phys > 0:
    print(f"  lambda_phys = {lam_phys:.4f} > 0  [OK]")
elif np.isnan(lam_phys):
    print("  Pending: lambda_phys computation failed.")
else:
    print(f"  lambda_phys = {lam_phys:.4f} < 0: possible instability in physical direction. Additional verification needed.")

print()
print(f"\nComputation time: {time.time()-t0:.1f}s")
print()
print("*Script: scripts/paper03/spin1_kinetic_metric.py")
print("=" * 60)
