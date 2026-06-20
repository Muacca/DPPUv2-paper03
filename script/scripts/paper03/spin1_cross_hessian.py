"""
Cross-Hessian between physical spin-1 mode X_i and other sectors (App F.7).

Targets:
  1. spin-1 × spin-2  ∂²V/∂q_A ∂X_i  (A=1..5, i=0..2)  → 15 components
  2. spin-1 × radial  ∂²V/∂r  ∂X_i   (i=0..2)           →  3 components
  3. spin-1 × torsion ∂²V/∂η  ∂X_i   (i=0..2)           →  3 components
                                                  total 21 components

Physical spin-1 mode (null mode analysis result (App E.6)):
    X_k = (-2 omega_k + 3 delta_k) / sqrt(13),  k=0,1,2

Settings: r₀=3, L=κ=η=1, α=θ_NY=0  (known values (App E.6), Mode.AX)

Method: pure numpy V_eff function + central finite difference (h=1e-3)

Verification: before computing cross-Hessian, reproduce known values (App E.6)
    μ²_ω ≈ 52.64,  μ²_δ ≈ 3789.93,  m_ωδ ≈ -631.65

Author: Muacca
Date: 2026-03-15
"""

from __future__ import annotations
import io
import sys

import numpy as np
from itertools import permutations


def _configure_utf8_stdio() -> None:
    """Avoid cp932 console failures on Windows when printing Unicode labels."""
    for name in ("stdout", "stderr"):
        stream = getattr(sys, name)
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
            continue
        except Exception:
            pass

        buffer = getattr(stream, "buffer", None)
        if buffer is not None:
            setattr(sys, name, io.TextIOWrapper(buffer, encoding="utf-8", errors="replace"))


_configure_utf8_stdio()

# =============================================================================
# Physical parameters  (known values (App E.6))
# =============================================================================
R0       = 3.0
L_VAL    = 1.0
KAPPA    = 1.0
ETA      = 1.0
ALPHA    = 0.0
THETA_NY = 0.0

# Physical spin-1 mode direction in (omega, sd=sin(delta)) space
# X_k = (DX_OMEGA * omega_k + DX_DELTA * sd_k) / (perturbation amplitude)
DX_OMEGA = -2.0 / np.sqrt(13.0)   # d(omega_k)/d(X_k)
DX_DELTA = +3.0 / np.sqrt(13.0)   # d(sd_k)/d(X_k)  [sd ≈ delta at iso point]

# Finite difference step
H = 1e-3

# =============================================================================
# Levi-Civita symbols  (computed once at module load)
# =============================================================================

def _make_eps3() -> np.ndarray:
    eps3 = np.zeros((3, 3, 3))
    eps3[0, 1, 2] = eps3[1, 2, 0] = eps3[2, 0, 1] = +1.0
    eps3[0, 2, 1] = eps3[2, 1, 0] = eps3[1, 0, 2] = -1.0
    return eps3

def _make_eps4() -> np.ndarray:
    """4D Levi-Civita ε_{abcd} (all lower, Euclidean convention)."""
    eps4 = np.zeros((4, 4, 4, 4))
    for p in permutations(range(4)):
        inv = sum(1 for i in range(4) for j in range(i + 1, 4) if p[i] > p[j])
        eps4[p] = +1.0 if inv % 2 == 0 else -1.0
    return eps4

EPS3 = _make_eps3()
EPS4 = _make_eps4()

# =============================================================================
# Mixing rotation matrix  M = R03 @ R13 @ R23
# (see base_topology.py _apply_s3_mixing_rotation_to_C)
# =============================================================================

def _make_mixing_matrix(sd0: float, sd1: float, sd2: float) -> np.ndarray:
    """
    Build 4×4 mixing rotation M = R₀₃ @ R₁₃ @ R₂₃.

    sd_i = sin(θ_i) = δ_i / sqrt(1+δ_i²)
    cd_i = cos(θ_i) = 1   / sqrt(1+δ_i²) = sqrt(1-sd_i²)

    R_{i3}: rotation in (ê^i, ê^3) plane:
        row i: [...,  cd_i, ..., sd_i, ...]
        row 3: [..., -sd_i, ..., cd_i, ...]
    """
    cd0 = np.sqrt(max(0.0, 1.0 - sd0**2))
    cd1 = np.sqrt(max(0.0, 1.0 - sd1**2))
    cd2 = np.sqrt(max(0.0, 1.0 - sd2**2))

    R03 = np.eye(4)
    R03[0, 0] =  cd0;  R03[0, 3] =  sd0
    R03[3, 0] = -sd0;  R03[3, 3] =  cd0

    R13 = np.eye(4)
    R13[1, 1] =  cd1;  R13[1, 3] =  sd1
    R13[3, 1] = -sd1;  R13[3, 3] =  cd1

    R23 = np.eye(4)
    R23[2, 2] =  cd2;  R23[2, 3] =  sd2
    R23[3, 2] = -sd2;  R23[3, 3] =  cd2

    return R03 @ R13 @ R23

# =============================================================================
# Core V_eff computation (numpy)
# =============================================================================

def compute_Veff(
    r:       float,
    eps:     float,   # squash  (q₁ ~ T₁ mode)
    s:       float,   # shear   (q₂ ~ T₂ mode)
    q3:      float,   # off-diagonal shear T₃  (01-plane)
    q4:      float,   # off-diagonal shear T₄  (02-plane)
    q5:      float,   # off-diagonal shear T₅  (12-plane)
    omega1:  float,   # S¹ twist along σ⁰
    omega2:  float,   # S¹ twist along σ¹
    omega3:  float,   # S¹ twist along σ²
    sd0:     float,   # sin(δ₀): mixing in (ê⁰,ê³) plane
    sd1:     float,   # sin(δ₁): mixing in (ê¹,ê³) plane
    sd2:     float,   # sin(δ₂): mixing in (ê²,ê³) plane
    L:       float,
    kappa:   float,
    eta:     float,   # axial torsion parameter (Mode.AX)
    theta_NY: float = 0.0,
) -> float:
    """
    V_eff for S³×S¹ with full DOF set.

    Pipeline:
      1. Spatial frame factors f₀,f₁,f₂ from (ε,s)
      2. Off-diagonal rotation G = G₃ @ G₄ @ G₅ from (q₃,q₄,q₅)
      3. Full frame F = diag(f) @ G → C^i_{jk} = (4/r) Σ F[il] F⁻¹[jm] F⁻¹[kn] ε_{lmn}
      4. Twist terms: C^3_{jk} from (ω₁,ω₂,ω₃)
      5. Mixing rotation M(δ) applied to all C (tensor transformation)
      6. Torsion T^i_{jk} = (2η/r) ε_{ijk}  (Mode.AX)
      7. LC connection Γ_LC = (1/2)(C + C^T210 − C^T102)
      8. Contortion K = (1/2)(T + T^T201 − T^T120)
      9. EC connection Γ_EC = Γ_LC + K
     10. Riemann R^a_{bcd} = Γ^a_{ec}Γ^e_{bd} − Γ^a_{ed}Γ^e_{bc} + Γ^a_{be}C^e_{cd}
     11. Ricci scalar R = Σ_{a,b} R^a_{bab}
     12. (Optional) Nieh-Yan N  [skipped if theta_NY=0]
     13. V_eff = −(R/2κ² + θ_NY·N) × Vol,  Vol = 2π²Lr³
    """
    # ── Step 1: diagonal frame factors ───────────────────────────────────────
    f0 = (1.0 + eps)**(2.0/3.0) * (1.0 + s)**2
    f1 = (1.0 + eps)**(2.0/3.0) / (1.0 + s)**2
    f2 = (1.0 + eps)**(-4.0/3.0)

    # ── Step 2: off-diagonal rotation G = G₃ @ G₄ @ G₅ ─────────────────────
    sq2 = np.sqrt(2.0)
    c3, sh3 = np.cosh(q3 / sq2), np.sinh(q3 / sq2)
    c4, sh4 = np.cosh(q4 / sq2), np.sinh(q4 / sq2)
    c5, sh5 = np.cosh(q5 / sq2), np.sinh(q5 / sq2)

    G3 = np.array([[c3, sh3, 0.0],
                   [sh3, c3, 0.0],
                   [0.0, 0.0, 1.0]])
    G4 = np.array([[c4,  0.0, sh4],
                   [0.0, 1.0, 0.0],
                   [sh4, 0.0, c4 ]])
    G5 = np.array([[1.0, 0.0, 0.0 ],
                   [0.0, c5,  sh5 ],
                   [0.0, sh5, c5  ]])
    G = G3 @ G4 @ G5

    # ── Step 3: full frame F and spatial structure constants ─────────────────
    F     = np.diag([f0, f1, f2]) @ G
    F_inv = np.linalg.inv(F)

    C4 = np.zeros((4, 4, 4))
    # C^i_{jk} = (4/r) Σ_{l,m,n} F[i,l] F⁻¹[j,m] F⁻¹[k,n] ε_{lmn}
    C4[:3, :3, :3] = (4.0 / r) * np.einsum('il,jm,kn,lmn->ijk',
                                             F, F_inv, F_inv, EPS3)

    # ── Step 4: twist terms  C^3_{jk} from (ω₁,ω₂,ω₃) ─────────────────────
    # From base_topology.py _add_s3_twist_C:
    #   ê³ = L(dτ + ω₁σ⁰ + ω₂σ¹ + ω₃σ²)
    #   C³_{12} = +4L ω₁ (1+ε)^{1/3}(1+s)  / r²
    #   C³_{20} = +4L ω₂ (1+ε)^{1/3}/(1+s) / r²
    #   C³_{01} = +4L ω₃ (1+ε)^{−2/3}       / r²
    t12 = (1.0 + eps)**(1.0/3.0) * (1.0 + s)
    t20 = (1.0 + eps)**(1.0/3.0) / (1.0 + s)
    t01 = (1.0 + eps)**(-2.0/3.0)

    C4[3, 1, 2] = +4.0 * L * omega1 * t12 / r**2
    C4[3, 2, 1] = -4.0 * L * omega1 * t12 / r**2
    C4[3, 2, 0] = +4.0 * L * omega2 * t20 / r**2
    C4[3, 0, 2] = -4.0 * L * omega2 * t20 / r**2
    C4[3, 0, 1] = +4.0 * L * omega3 * t01 / r**2
    C4[3, 1, 0] = -4.0 * L * omega3 * t01 / r**2

    # ── Step 5: mixing rotation  M(δ₀,δ₁,δ₂) applied to C ──────────────────
    # C'[a,b,c] = Σ M[a,d] M[e,b] M[f,c] C_old[d,e,f]
    # (M is orthogonal: M⁻¹ = Mᵀ, so (M⁻¹)[b,e] = M[e,b])
    if not (sd0 == 0.0 and sd1 == 0.0 and sd2 == 0.0):
        M  = _make_mixing_matrix(sd0, sd1, sd2)
        C4 = np.einsum('ad,eb,fc,def->abc', M, M, M, C4)

    # ── Step 6: torsion T^a_{bc}  (Mode.AX: T^i_{jk} = 2η/r ε_{ijk}) ───────
    T4 = np.zeros((4, 4, 4))
    T4[:3, :3, :3] = (2.0 * eta / r) * EPS3

    # ── Step 7: LC connection Γ^a_{bc} = (1/2)(C^a_{bc} + C^c_{ba} − C^b_{ac}) ─
    Gamma_LC = 0.5 * (C4
                      + np.transpose(C4, (2, 1, 0))   # C^c_{ba}
                      - np.transpose(C4, (1, 0, 2)))   # C^b_{ac}

    # ── Step 8: contortion K^a_{bc} = (1/2)(T^a_{bc} + T^b_{ca} − T^c_{ab}) ─
    K = 0.5 * (T4
               + np.transpose(T4, (2, 0, 1))   # T^b_{ca}
               - np.transpose(T4, (1, 2, 0)))   # T^c_{ab}

    # ── Step 9: EC connection ─────────────────────────────────────────────────
    Gamma = Gamma_LC + K

    # ── Step 10: Riemann tensor R^a_{bcd} ─────────────────────────────────────
    # R^a_{bcd} = Γ^a_{ec}Γ^e_{bd} − Γ^a_{ed}Γ^e_{bc} + Γ^a_{be}C^e_{cd}
    Riem = (np.einsum('aec,ebd->abcd', Gamma, Gamma)
            - np.einsum('aed,ebc->abcd', Gamma, Gamma)
            + np.einsum('abe,ecd->abcd', Gamma, C4))

    # ── Step 11: Ricci scalar R = Σ_{a,b} R^a_{bab} ─────────────────────────
    R_scalar = np.einsum('abab->', Riem)

    # ── Step 12: Nieh-Yan (skipped when θ_NY=0) ──────────────────────────────
    N = 0.0
    if theta_NY != 0.0:
        N_TT  = 0.25 * np.einsum('abcd,eab,ecd->', EPS4, T4, T4)
        N_Ree = 0.25 * np.einsum('abcd,abcd->', EPS4, Riem)
        N = N_TT - N_Ree

    # ── Step 13: V_eff = −(R/2κ² + θ_NY·N) × Vol ────────────────────────────
    Vol  = 2.0 * np.pi**2 * L * r**3
    return -(R_scalar / (2.0 * kappa**2) + theta_NY * N) * Vol


# =============================================================================
# Helper: evaluate V_eff at isotropic point with override kwargs
# =============================================================================
_ISO = dict(
    r=R0, eps=0.0, s=0.0, q3=0.0, q4=0.0, q5=0.0,
    omega1=0.0, omega2=0.0, omega3=0.0,
    sd0=0.0, sd1=0.0, sd2=0.0,
    L=L_VAL, kappa=KAPPA, eta=ETA, theta_NY=THETA_NY,
)

def veff(**overrides) -> float:
    return compute_Veff(**{**_ISO, **overrides})

# =============================================================================
# X_k perturbation: when X_k amplitude changes by h,
#   omega_{k+1} changes by DX_OMEGA * h  (1-indexed: omega1,omega2,omega3)
#   sd_{k}      changes by DX_DELTA * h  (0-indexed: sd0, sd1, sd2)
# =============================================================================
_OMEGA_KEYS = ['omega1', 'omega2', 'omega3']
_SD_KEYS    = ['sd0',    'sd1',    'sd2'   ]

def _perturb_X(k: int, h: float) -> dict:
    """Return kwargs for X_k perturbation of amplitude h."""
    return {_OMEGA_KEYS[k]: DX_OMEGA * h,
            _SD_KEYS[k]:    DX_DELTA * h}

# =============================================================================
# Cross-Hessian by central finite difference
# ∂²V / ∂p ∂X_k  at isotropic point
#
#  p is perturbed around its isotropic-point value p₀:
#    p → p₀ ± h_p
#  X_k is perturbed by ± h_X (from 0):
#    omega_{k+1} → ± DX_OMEGA * h_X
#    sd_{k}      → ± DX_DELTA * h_X
#
# Formula:
#   [V(p₀+h,+hX) − V(p₀+h,−hX) − V(p₀−h,+hX) + V(p₀−h,−hX)] / (4h²)
# =============================================================================

def cross_hessian(p_key: str, p0: float, k: int, h: float = H) -> float:
    """Compute ∂²V/∂p ∂X_k via central FD around isotropic point."""
    xp = _perturb_X(k, +h)
    xm = _perturb_X(k, -h)
    vpp = veff(**{p_key: p0 + h}, **xp)
    vpm = veff(**{p_key: p0 + h}, **xm)
    vmp = veff(**{p_key: p0 - h}, **xp)
    vmm = veff(**{p_key: p0 - h}, **xm)
    return (vpp - vpm - vmp + vmm) / (4.0 * h * h)

# =============================================================================
# Main
# =============================================================================

def main():
    SEP = "=" * 72
    sep = "-" * 72

    print(SEP)
    print("Quadratic mixing of physical spin-1 with other sectors (App F.7)")
    print("X_k = (-2ω_k + 3δ_k)/√13  (null mode analysis result (App E.6))")
    print(f"Settings: r₀={R0}, L={L_VAL}, κ={KAPPA}, η={ETA}, α={ALPHA}, θ_NY={THETA_NY}")
    print(f"Finite-difference step h = {H}")
    print(SEP)

    # ── Sanity checks: reproduce known values (App E.6) ──────────────────────
    print("\n── Sanity check: reproducing known values (App E.6) ──")
    V0 = veff()
    print(f"V_eff (isotropic):        {V0:.6f}")

    vwp = veff(omega1=+H);  vwm = veff(omega1=-H)
    mu2_w = (vwp - 2*V0 + vwm) / H**2
    print(f"μ²_ω  (∂²V/∂ω₁²):        {mu2_w:.4f}   [expect ≈  52.64]")

    vdp = veff(sd0=+H);     vdm = veff(sd0=-H)
    mu2_d = (vdp - 2*V0 + vdm) / H**2
    print(f"μ²_δ  (∂²V/∂sd₀²):       {mu2_d:.4f}   [expect ≈ 3789.93]")

    vpp_ = veff(omega1=+H, sd0=+H)
    vpm_ = veff(omega1=+H, sd0=-H)
    vmp_ = veff(omega1=-H, sd0=+H)
    vmm_ = veff(omega1=-H, sd0=-H)
    m_wd = (vpp_ - vpm_ - vmp_ + vmm_) / (4*H**2)
    print(f"m_ωδ  (∂²V/∂ω₁∂sd₀):    {m_wd:.4f}   [expect ≈ -631.65]")

    ok = (abs(mu2_w -   52.64) < 0.1 and
          abs(mu2_d - 3789.93) < 1.0 and
          abs(m_wd  - (-631.65)) < 1.0)
    print(f"Sanity check: {'OK ✓' if ok else 'NG ✗ — numerical mismatch, check results'}")

    # ── §3: spin-1 × spin-2 ──────────────────────────────────────────────────
    print()
    print(sep)
    print("§3: spin-1 × spin-2  M_{Ai}^{(qX)} = ∂²V/∂q_A ∂X_i")
    print(sep)
    print(f"{'A':>2}  {'param':>15}  {'X₀ (k=0)':>12}  {'X₁ (k=1)':>12}  {'X₂ (k=2)':>12}")
    print("-" * 60)

    spin2_specs = [
        ('eps', 0.0, 'ε (T₁, squash)'),
        ('s',   0.0, 's (T₂, shear)'),
        ('q3',  0.0, 'q₃ (T₃, 01-shear)'),
        ('q4',  0.0, 'q₄ (T₄, 02-shear)'),
        ('q5',  0.0, 'q₅ (T₅, 12-shear)'),
    ]

    all_spin2: list[float] = []
    for A, (p_key, p0, label) in enumerate(spin2_specs, start=1):
        row = [cross_hessian(p_key, p0, k) for k in range(3)]
        all_spin2.extend(row)
        print(f"{A:>2}  {label:>15}  {row[0]:>12.4f}  {row[1]:>12.4f}  {row[2]:>12.4f}")

    max_s2 = max(abs(v) for v in all_spin2)
    print(f"\n  Max |M_Ai^(qX)| = {max_s2:.4e}")
    print(f"  → {'DECOUPLE confirmed' if max_s2 < 1.0 else 'NON-ZERO detected'}")

    # ── §4: spin-1 × radial ──────────────────────────────────────────────────
    print()
    print(sep)
    print("§4: spin-1 × radial  M_i^{(rX)} = ∂²V/∂r ∂X_i")
    print(sep)
    print(f"{'i':>3}  {'M_i^(rX)':>16}")
    print("-" * 22)

    all_radial: list[float] = []
    for k in range(3):
        val = cross_hessian('r', R0, k)
        all_radial.append(val)
        print(f"{k+1:>3}  {val:>16.4f}")

    max_r = max(abs(v) for v in all_radial)
    print(f"\n  Max |M_i^(rX)| = {max_r:.4e}")
    print(f"  → {'DECOUPLE confirmed' if max_r < 1.0 else 'NON-ZERO detected'}")

    # ── §5: spin-1 × torsion (η) ──────────────────────────────────────────────
    print()
    print(sep)
    print("§5: spin-1 × torsion (η)  M_αi^{(TX)} = ∂²V/∂η ∂X_i")
    print("    (Mode.AX: η only, since V=0)")
    print(sep)
    print(f"{'torsion':>10}  {'X₀ (k=0)':>12}  {'X₁ (k=1)':>12}  {'X₂ (k=2)':>12}")
    print("-" * 50)

    all_torsion: list[float] = []
    row_eta = [cross_hessian('eta', ETA, k) for k in range(3)]
    all_torsion.extend(row_eta)
    print(f"{'η':>10}  {row_eta[0]:>12.4f}  {row_eta[1]:>12.4f}  {row_eta[2]:>12.4f}")

    max_t = max(abs(v) for v in all_torsion)
    print(f"\n  Max |M_αi^(TX)| = {max_t:.4e}")
    print(f"  → {'DECOUPLE confirmed' if max_t < 1.0 else 'NON-ZERO detected'}")

    # ── Final judgment ─────────────────────────────────────────────────────────
    print()
    print(SEP)
    print("Final judgment")
    print(SEP)
    all_vals = all_spin2 + all_radial + all_torsion
    max_all  = max(abs(v) for v in all_vals)
    print(f"Max |cross-Hessian| over all {len(all_vals)} components = {max_all:.4e}")
    print()
    if max_all < 1.0:
        print("Case A: all mixed Hessians at zero level")
        print("  physical spin-1 mode X_k decouples from other sectors at quadratic order")
        print("  spin-1 stability is fixed by the separate Schur-reduced X-Y block")
    else:
        print("Case B: non-zero components found")
        print("  stability cannot be determined from spin-1 sector alone")
        print("  full coupled Hessian required")

    print()
    print("── Numerical details ──")
    labels_s2 = [f"M({spec[2][:2].strip(',')},X{k})" for spec in spin2_specs for k in range(3)]
    labels_r  = [f"M(r,X{k})"   for k in range(3)]
    labels_t  = [f"M(η,X{k})"   for k in range(3)]
    labels_all = labels_s2 + labels_r + labels_t
    for lbl, val in sorted(zip(labels_all, all_vals), key=lambda x: abs(x[1]), reverse=True)[:5]:
        print(f"  {lbl:>15}: {val:>12.4e}")


if __name__ == '__main__':
    import os as _os, sys as _sys
    _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..'))
    from dppu.utils.tee_logger import setup_log, teardown_log
    setup_log(__file__, log_dir='output')
    main()
    teardown_log()
