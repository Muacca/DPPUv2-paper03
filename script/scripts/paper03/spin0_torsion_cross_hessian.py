"""
Cross-Hessian in the spin-0 / torsion sector (Appendix H extension).

Targets  (paper03 mode dictionary, homogeneous minisuperspace, MX background):

  Priority A: spin-2 × torsion
      ∂²V/∂q_A ∂η   (A=1..5, 5 components)
      ∂²V/∂q_A ∂V   (A=1..5, 5 components)
                     ----------
                     10 total

  Priority B: radial × torsion
      ∂²V/∂r ∂η                →  1 component
      ∂²V/∂r ∂V                →  1 component
                     ----------
                      2 total

  Priority C: torsion block H^(η,V)
      ∂²V/∂η²,  ∂²V/∂η∂V,  ∂²V/∂V²   →  3 independent components

Background: r₀=3, L=κ=η=1, V=0 (MX mode, V=0 at isotropic equal-squash point)
Settings:   α=θ_NY=0  (pure EC, no Nieh-Yan correction)

Note on torsion parameterisation (MX mode):
  Axial component:        T^i_{jk} = (2η/r) ε_{ijk}          for i,j,k ∈ {0,1,2}
  Vector-trace component: T^a_{3a} = +V/3,  T^a_{a3} = −V/3  for a=0,1,2
    derived from  T_{abc} = (1/3)(δ_{ac} V_b − δ_{ab} V_c),  V_μ = V·δ^3_μ

Method: pure numpy V_eff + central 4-point finite difference (h=1e-3)

Sanity check: reproduce known spin-1 sector values (App E.6)
    μ²_ω ≈ 52.64,  μ²_δ ≈ 3789.93,  m_ωδ ≈ −631.65

Author: Muacca
Date: 2026-03-18
"""

from __future__ import annotations

import io
import sys
import time as _time

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
# Physical parameters  (isotropic equal-squash background, App E.6)
# =============================================================================

R0        = 3.0
L_VAL     = 1.0
KAPPA     = 1.0
ETA_BG    = 1.0    # axial torsion background
V_BG      = 0.0    # vector-trace torsion background (V=0 at iso, MX mode)
ALPHA     = 0.0
THETA_NY  = 0.0

# Physical spin-1 mode direction (for sanity check, App E.6)
DX_OMEGA = -2.0 / np.sqrt(13.0)
DX_DELTA = +3.0 / np.sqrt(13.0)

# Finite-difference step
H = 1e-3

# =============================================================================
# Levi-Civita symbols
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
# =============================================================================

def _make_mixing_matrix(sd0: float, sd1: float, sd2: float) -> np.ndarray:
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
# Core V_eff computation (numpy, Mode.MX: axial η + vector-trace V)
# =============================================================================

def compute_Veff(
    r:          float,
    eps:        float,    # squash         (q₁ ~ T₁ mode)
    s:          float,    # shear          (q₂ ~ T₂ mode)
    q3:         float,    # off-diag shear T₃  (01-plane)
    q4:         float,    # off-diag shear T₄  (02-plane)
    q5:         float,    # off-diag shear T₅  (12-plane)
    omega1:     float,    # S¹ twist along σ⁰
    omega2:     float,    # S¹ twist along σ¹
    omega3:     float,    # S¹ twist along σ²
    sd0:        float,    # sin(δ₀)
    sd1:        float,    # sin(δ₁)
    sd2:        float,    # sin(δ₂)
    L:          float,
    kappa:      float,
    eta:        float,          # axial torsion η    (Mode.AX / Mode.MX)
    V_torsion:  float = 0.0,   # vector-trace V     (Mode.VT / Mode.MX)
    theta_NY:   float = 0.0,
) -> float:
    """
    V_eff for S³×S¹, Mode.MX (axial η + vector-trace V torsion).

    Pipeline identical to spin1_cross_hessian.compute_Veff, extended at Step 6:
      Axial:         T^i_{jk} = (2η/r) ε_{ijk}          i,j,k ∈ {0,1,2}
      Vector-trace:  T^a_{3a} = +V/3,  T^a_{a3} = −V/3  a=0,1,2
        (identity-metric frame; derived from T_{abc}=(1/3)(δ_{ac}V_b−δ_{ab}V_c))
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

    G3 = np.array([[c3, sh3, 0.0], [sh3, c3, 0.0], [0.0, 0.0, 1.0]])
    G4 = np.array([[c4, 0.0, sh4], [0.0, 1.0, 0.0], [sh4, 0.0, c4]])
    G5 = np.array([[1.0, 0.0, 0.0], [0.0, c5, sh5], [0.0, sh5, c5]])
    G = G3 @ G4 @ G5

    # ── Step 3: full frame F and spatial structure constants ─────────────────
    F     = np.diag([f0, f1, f2]) @ G
    F_inv = np.linalg.inv(F)

    C4 = np.zeros((4, 4, 4))
    C4[:3, :3, :3] = (4.0 / r) * np.einsum('il,jm,kn,lmn->ijk',
                                             F, F_inv, F_inv, EPS3)

    # ── Step 4: twist terms  C^3_{jk} from (ω₁,ω₂,ω₃) ─────────────────────
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
    if not (sd0 == 0.0 and sd1 == 0.0 and sd2 == 0.0):
        M  = _make_mixing_matrix(sd0, sd1, sd2)
        C4 = np.einsum('ad,eb,fc,def->abc', M, M, M, C4)

    # ── Step 6: torsion T^a_{bc}  (Mode.MX) ─────────────────────────────────
    T4 = np.zeros((4, 4, 4))

    # Axial component: T^i_{jk} = (2η/r) ε_{ijk}
    T4[:3, :3, :3] = (2.0 * eta / r) * EPS3

    # Vector-trace component: T^a_{3a} = +V/3,  T^a_{a3} = −V/3  (a=0,1,2)
    # Derived from T_{abc} = (1/3)(δ_{ac} V_b − δ_{ab} V_c), V_μ = V·δ^3_μ
    if V_torsion != 0.0:
        v3 = V_torsion / 3.0
        for a in range(3):
            T4[a, 3, a] += v3
            T4[a, a, 3] -= v3

    # ── Step 7: LC connection Γ^a_{bc} = (1/2)(C^a_{bc} + C^c_{ba} − C^b_{ac}) ─
    Gamma_LC = 0.5 * (C4
                      + np.transpose(C4, (2, 1, 0))
                      - np.transpose(C4, (1, 0, 2)))

    # ── Step 8: contortion K^a_{bc} = (1/2)(T^a_{bc} + T^b_{ca} − T^c_{ab}) ─
    K = 0.5 * (T4
               + np.transpose(T4, (2, 0, 1))
               - np.transpose(T4, (1, 2, 0)))

    # ── Step 9: EC connection ─────────────────────────────────────────────────
    Gamma = Gamma_LC + K

    # ── Step 10: Riemann tensor R^a_{bcd} ─────────────────────────────────────
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
    Vol = 2.0 * np.pi**2 * L * r**3
    return -(R_scalar / (2.0 * kappa**2) + theta_NY * N) * Vol


# =============================================================================
# Helper: evaluate V_eff at background point with override kwargs
# =============================================================================

_BG = dict(
    r=R0, eps=0.0, s=0.0, q3=0.0, q4=0.0, q5=0.0,
    omega1=0.0, omega2=0.0, omega3=0.0,
    sd0=0.0, sd1=0.0, sd2=0.0,
    L=L_VAL, kappa=KAPPA,
    eta=ETA_BG, V_torsion=V_BG, theta_NY=THETA_NY,
)


def veff(**overrides) -> float:
    return compute_Veff(**{**_BG, **overrides})


# =============================================================================
# Finite-difference helpers
# =============================================================================

def hess_diag(p_key: str, p0: float, h: float = H) -> float:
    """∂²V/∂p² at background via 3-point central FD."""
    return (veff(**{p_key: p0 + h})
            - 2.0 * veff(**{p_key: p0})
            + veff(**{p_key: p0 - h})) / h**2


def hess_cross(p1_key: str, p1_0: float,
               p2_key: str, p2_0: float, h: float = H) -> float:
    """∂²V/∂p₁∂p₂ at background via 4-point central FD."""
    vpp = veff(**{p1_key: p1_0 + h, p2_key: p2_0 + h})
    vpm = veff(**{p1_key: p1_0 + h, p2_key: p2_0 - h})
    vmp = veff(**{p1_key: p1_0 - h, p2_key: p2_0 + h})
    vmm = veff(**{p1_key: p1_0 - h, p2_key: p2_0 - h})
    return (vpp - vpm - vmp + vmm) / (4.0 * h**2)


# =============================================================================
# Main
# =============================================================================

def main():
    t_start = _time.time()
    SEP = "=" * 72
    sep = "-" * 72

    print(SEP)
    print("Spin-0 / torsion cross-Hessian (Appendix H extension)")
    print("Priority A: spin-2 x torsion (10),  B: radial x torsion (2),")
    print("            C: torsion block H^(eta,V) (3)")
    print(f"Background: r0={R0}, L={L_VAL}, kappa={KAPPA}, eta={ETA_BG}, V={V_BG}")
    print(f"            alpha={ALPHA}, theta_NY={THETA_NY},  h={H}")
    print("Mode: MX (axial eta + vector-trace V torsion)")
    print(SEP)

    # ── Sanity check: reproduce known spin-1 values (App E.6) ─────────────────
    #    Confirms that the extended compute_Veff (V_torsion=0 at background)
    #    is identical to spin1_cross_hessian.compute_Veff.
    print("\n── Sanity check: reproducing known spin-1 values (App E.6) ──")
    V0 = veff()
    print(f"V_eff (background):       {V0:.6f}")

    vwp = veff(omega1=+H);  vwm = veff(omega1=-H)
    mu2_w = (vwp - 2 * V0 + vwm) / H**2
    print(f"mu2_omega (d2V/domega1^2):  {mu2_w:.4f}   [expect ~  52.64]")

    vdp = veff(sd0=+H);     vdm = veff(sd0=-H)
    mu2_d = (vdp - 2 * V0 + vdm) / H**2
    print(f"mu2_delta (d2V/dsd0^2):     {mu2_d:.4f}   [expect ~ 3789.93]")

    vpp_ = veff(omega1=+H, sd0=+H)
    vpm_ = veff(omega1=+H, sd0=-H)
    vmp_ = veff(omega1=-H, sd0=+H)
    vmm_ = veff(omega1=-H, sd0=-H)
    m_wd = (vpp_ - vpm_ - vmp_ + vmm_) / (4 * H**2)
    print(f"m_od   (d2V/domega1 dsd0):  {m_wd:.4f}   [expect ~ -631.65]")

    ok = (abs(mu2_w -    52.64) < 0.1 and
          abs(mu2_d - 3789.93) < 1.0 and
          abs(m_wd  - (-631.65)) < 1.0)
    print(f"Sanity check: {'OK ✓' if ok else 'NG ✗ — numerical mismatch, check results'}")
    if not ok:
        print("  WARNING: pipeline may not match App E.6 reference — proceed with caution.")

    # ── Spin-2 sector specs (q_A coordinates) ─────────────────────────────────
    spin2_specs = [
        ('eps', 0.0, 'eps  (T1, squash)'),
        ('s',   0.0, 's    (T2, shear)'),
        ('q3',  0.0, 'q3   (T3, 01-off-diag)'),
        ('q4',  0.0, 'q4   (T4, 02-off-diag)'),
        ('q5',  0.0, 'q5   (T5, 12-off-diag)'),
    ]

    # ── Priority A-1: spin-2 × η  (5 components) ──────────────────────────────
    print()
    print(sep)
    print("Priority A-1:  spin-2 x eta   d2V/dq_A deta   (A=1..5, 5 components)")
    print(sep)
    print(f"  {'A':>2}  {'parameter':>22}  {'d2V/dq_A deta':>16}  verdict")
    print("  " + "-" * 62)

    results_s2_eta: list[float] = []
    for A, (pk, p0, label) in enumerate(spin2_specs, start=1):
        val = hess_cross(pk, p0, 'eta', ETA_BG)
        results_s2_eta.append(val)
        verdict = "zero" if abs(val) < 1.0 else f"NON-ZERO"
        print(f"  {A:>2}  {label:>22}  {val:>16.6f}  {verdict}")

    max_s2_eta = max(abs(v) for v in results_s2_eta)
    print(f"\n  Max |d2V/dq_A deta| = {max_s2_eta:.4e}"
          f"  -> {'all zero ✓' if max_s2_eta < 1.0 else 'NON-ZERO ✗'}")

    # ── Priority A-2: spin-2 × V  (5 components) ──────────────────────────────
    print()
    print(sep)
    print("Priority A-2:  spin-2 x V   d2V/dq_A dV   (A=1..5, 5 components)")
    print(sep)
    print(f"  {'A':>2}  {'parameter':>22}  {'d2V/dq_A dV':>16}  verdict")
    print("  " + "-" * 62)

    results_s2_V: list[float] = []
    for A, (pk, p0, label) in enumerate(spin2_specs, start=1):
        val = hess_cross(pk, p0, 'V_torsion', V_BG)
        results_s2_V.append(val)
        verdict = "zero" if abs(val) < 1.0 else f"NON-ZERO"
        print(f"  {A:>2}  {label:>22}  {val:>16.6f}  {verdict}")

    max_s2_V = max(abs(v) for v in results_s2_V)
    print(f"\n  Max |d2V/dq_A dV| = {max_s2_V:.4e}"
          f"  -> {'all zero ✓' if max_s2_V < 1.0 else 'NON-ZERO ✗'}")

    # ── Priority B: radial × torsion  (2 components) ──────────────────────────
    print()
    print(sep)
    print("Priority B:  radial x torsion   d2V/dr deta,  d2V/dr dV   (2 components)")
    print(sep)
    print(f"  {'component':>16}  {'value':>18}  verdict")
    print("  " + "-" * 48)

    r_eta = hess_cross('r', R0, 'eta', ETA_BG)
    r_V   = hess_cross('r', R0, 'V_torsion', V_BG)

    for label, val in [('d2V/dr deta', r_eta), ('d2V/dr dV', r_V)]:
        verdict = "zero" if abs(val) < 1.0 else "NON-ZERO"
        print(f"  {label:>16}  {val:>18.6f}  {verdict}")

    max_radial = max(abs(r_eta), abs(r_V))
    print(f"\n  Max |d2V/dr d{{eta,V}}| = {max_radial:.4e}"
          f"  -> {'all zero ✓' if max_radial < 1.0 else 'mixing found ✗'}")

    # ── Priority C: torsion block H^(η,V)  (3 components) ─────────────────────
    print()
    print(sep)
    print("Priority C:  torsion block H^(eta,V)   (3 independent components)")
    print(sep)
    print(f"  {'component':>16}  {'value':>18}")
    print("  " + "-" * 38)

    H_etaeta = hess_diag('eta', ETA_BG)
    H_etaV   = hess_cross('eta', ETA_BG, 'V_torsion', V_BG)
    H_VV     = hess_diag('V_torsion', V_BG)

    for label, val in [('d2V/deta^2',   H_etaeta),
                       ('d2V/deta dV',  H_etaV),
                       ('d2V/dV^2',     H_VV)]:
        print(f"  {label:>16}  {val:>18.6f}")

    det_H = H_etaeta * H_VV - H_etaV**2
    tr_H  = H_etaeta + H_VV
    disc  = max(0.0, (H_etaeta - H_VV)**2 + 4.0 * H_etaV**2)
    eig1  = 0.5 * (tr_H + np.sqrt(disc))
    eig2  = 0.5 * (tr_H - np.sqrt(disc))
    print(f"\n  det H^(eta,V) = {det_H:.6e}")
    print(f"  eigenvalues   = {eig1:.4f},  {eig2:.4f}")
    if eig1 > 0 and eig2 > 0:
        stab_str = "stable (both eigenvalues > 0)"
    elif eig1 > 0 and eig2 < 0:
        stab_str = "SADDLE (mixed signs)"
    else:
        stab_str = "UNSTABLE (both < 0)"
    print(f"  torsion block: {stab_str}")

    # ── Grand summary ──────────────────────────────────────────────────────────
    all_cross = results_s2_eta + results_s2_V + [r_eta, r_V]
    max_cross = max(abs(v) for v in all_cross)

    print()
    print(SEP)
    print("Grand summary")
    print(SEP)
    print(f"  spin-2 x eta  (5 components)  max |val| = {max_s2_eta:.4e}"
          f"  -> {'zero ✓' if max_s2_eta < 1.0 else 'NON-ZERO ✗'}")
    print(f"  spin-2 x V    (5 components)  max |val| = {max_s2_V:.4e}"
          f"  -> {'zero ✓' if max_s2_V < 1.0 else 'NON-ZERO ✗'}")
    print(f"  radial  x eta (1 component)   val = {r_eta:.6f}"
          f"  -> {'zero ✓' if abs(r_eta) < 1.0 else 'NON-ZERO ✗'}")
    print(f"  radial  x V   (1 component)   val = {r_V:.6f}"
          f"  -> {'zero ✓' if abs(r_V) < 1.0 else 'NON-ZERO ✗'}")
    print(f"  torsion block H^(eta,V)  det = {det_H:.4e}"
          f"  eigenvalues = ({eig1:.2f}, {eig2:.2f})  ({stab_str})")
    print(f"\n  Total cross-Hessian components checked: {len(all_cross)}")
    print(f"  Max |cross-Hessian| = {max_cross:.4e}")
    print()
    if max_cross < 1.0:
        print("  -> Conclusion: spin-2 <-> torsion and radial <-> torsion")
        print("     decouple at quadratic order.")
        print("     Block-diagonal EFT structure confirmed in spin-0/torsion sector.")
    else:
        nz_s2eta = sum(1 for v in results_s2_eta if abs(v) >= 1.0)
        nz_s2V   = sum(1 for v in results_s2_V   if abs(v) >= 1.0)
        nz_r     = sum(1 for v in [r_eta, r_V]   if abs(v) >= 1.0)
        print("  -> Conclusion: mixing found in the following blocks:")
        if nz_s2eta > 0: print(f"       spin-2 x eta: {nz_s2eta}/5 nonzero")
        if nz_s2V   > 0: print(f"       spin-2 x V:   {nz_s2V}/5 nonzero")
        if nz_r     > 0: print(f"       radial  x {{eta,V}}: {nz_r}/2 nonzero")
        print("     Full coupled Hessian required for stability analysis.")

    # ── Appendix H draft ──────────────────────────────────────────────────────
    print()
    print(SEP)
    print("Appendix H draft  (cross-Hessian summary table extension)")
    print(SEP)

    def _res(vals, tol: float = 1.0) -> str:
        nz = sum(1 for v in vals if abs(v) >= tol)
        return "all zero" if nz == 0 else f"{nz} nonzero"

    rows = [
        ("spin-1 x spin-2", 15, "all zero",            "(App F.7, established)"),
        ("spin-1 x radial",  3, "all zero",             "(App F.7, established)"),
        ("spin-1 x torsion", 3, "all zero",             "(App F.7, established)"),
        ("spin-2 x eta",     5, _res(results_s2_eta),   "Priority A-1 (this script)"),
        ("spin-2 x V",       5, _res(results_s2_V),     "Priority A-2 (this script)"),
        ("radial  x eta",    1, _res([r_eta]),           "Priority B   (this script)"),
        ("radial  x V",      1, _res([r_V]),             "Priority B   (this script)"),
    ]

    print(f"\n  {'block':<22}  {'# components':>14}  {'result':<12}  note")
    print("  " + "-" * 72)
    for block, n, res, note in rows:
        print(f"  {block:<22}  {n:>14}  {res:<12}  {note}")

    print()
    print(f"  Torsion block H^(eta,V) eigenvalues: {eig1:.2f}, {eig2:.2f}"
          f"  ({stab_str})")

    print()
    elapsed = _time.time() - t_start
    print(f"Computation time: {elapsed:.1f}s")
    print()
    print("*Script: scripts/paper03/spin0_torsion_cross_hessian.py")


if __name__ == '__main__':
    import os as _os, sys as _sys
    _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..'))
    from dppu.utils.tee_logger import setup_log, teardown_log
    setup_log(__file__, log_dir='output')
    main()
    teardown_log()
