"""
Eigenvalue analysis of the coupled (r, eta) 2x2 Hessian block.

The spin-0/torsion cross-Hessian check (spin0_torsion_cross_hessian.py) showed:
  - d2V/dr deta  != 0  (~118.4)
  - d2V/dr dV    =  0
  - H^(eta,V) is diagonal and degenerate (eigenvalues ~355.31)

This script constructs the 2x2 block

    H^(r,eta) = [[H_rr,   H_reta],
                 [H_etar, H_etaeta]]

at the isotropic MX background and evaluates its eigenvalues to determine
whether the scalar/torsion subsector is locally stable when coupling is included.

Background  (same as spin0_torsion_cross_hessian.py):
  r0=3, L=kappa=eta=1, V=0, alpha=theta_NY=0, mode=MX

Bonus: AX / VT / MX mode comparison to check whether the r-eta mixing
       is specific to the axial torsion component.

Method: pure numpy V_eff + central finite difference (h=1e-3)

Author: Muacca
Date: 2026-03-19
"""

from __future__ import annotations

import io
import sys
import time as _time

import numpy as np
from itertools import permutations


def _configure_utf8_stdio() -> None:
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
# Physical parameters
# =============================================================================

R0       = 3.0
L_VAL    = 1.0
KAPPA    = 1.0
ETA_BG   = 1.0
V_BG     = 0.0
THETA_NY = 0.0

H = 1e-3   # finite-difference step

# =============================================================================
# Levi-Civita symbols
# =============================================================================

def _make_eps3() -> np.ndarray:
    eps3 = np.zeros((3, 3, 3))
    eps3[0, 1, 2] = eps3[1, 2, 0] = eps3[2, 0, 1] = +1.0
    eps3[0, 2, 1] = eps3[2, 1, 0] = eps3[1, 0, 2] = -1.0
    return eps3


def _make_eps4() -> np.ndarray:
    eps4 = np.zeros((4, 4, 4, 4))
    for p in permutations(range(4)):
        inv = sum(1 for i in range(4) for j in range(i + 1, 4) if p[i] > p[j])
        eps4[p] = +1.0 if inv % 2 == 0 else -1.0
    return eps4


EPS3 = _make_eps3()
EPS4 = _make_eps4()

# =============================================================================
# Mixing rotation matrix
# =============================================================================

def _make_mixing_matrix(sd0: float, sd1: float, sd2: float) -> np.ndarray:
    cd0 = np.sqrt(max(0.0, 1.0 - sd0**2))
    cd1 = np.sqrt(max(0.0, 1.0 - sd1**2))
    cd2 = np.sqrt(max(0.0, 1.0 - sd2**2))
    R03 = np.eye(4); R03[0,0]=cd0; R03[0,3]=sd0; R03[3,0]=-sd0; R03[3,3]=cd0
    R13 = np.eye(4); R13[1,1]=cd1; R13[1,3]=sd1; R13[3,1]=-sd1; R13[3,3]=cd1
    R23 = np.eye(4); R23[2,2]=cd2; R23[2,3]=sd2; R23[3,2]=-sd2; R23[3,3]=cd2
    return R03 @ R13 @ R23

# =============================================================================
# Core V_eff computation (Mode.MX: axial eta + vector-trace V)
# =============================================================================

def compute_Veff(
    r:          float,
    eps:        float = 0.0,
    s:          float = 0.0,
    q3:         float = 0.0,
    q4:         float = 0.0,
    q5:         float = 0.0,
    omega1:     float = 0.0,
    omega2:     float = 0.0,
    omega3:     float = 0.0,
    sd0:        float = 0.0,
    sd1:        float = 0.0,
    sd2:        float = 0.0,
    L:          float = L_VAL,
    kappa:      float = KAPPA,
    eta:        float = ETA_BG,
    V_torsion:  float = V_BG,
    theta_NY:   float = THETA_NY,
) -> float:
    f0 = (1.0 + eps)**(2.0/3.0) * (1.0 + s)**2
    f1 = (1.0 + eps)**(2.0/3.0) / (1.0 + s)**2
    f2 = (1.0 + eps)**(-4.0/3.0)

    sq2 = np.sqrt(2.0)
    c3, sh3 = np.cosh(q3/sq2), np.sinh(q3/sq2)
    c4, sh4 = np.cosh(q4/sq2), np.sinh(q4/sq2)
    c5, sh5 = np.cosh(q5/sq2), np.sinh(q5/sq2)
    G3 = np.array([[c3,sh3,0.],[sh3,c3,0.],[0.,0.,1.]])
    G4 = np.array([[c4,0.,sh4],[0.,1.,0.],[sh4,0.,c4]])
    G5 = np.array([[1.,0.,0.],[0.,c5,sh5],[0.,sh5,c5]])
    G  = G3 @ G4 @ G5

    F     = np.diag([f0, f1, f2]) @ G
    F_inv = np.linalg.inv(F)

    C4 = np.zeros((4, 4, 4))
    C4[:3,:3,:3] = (4.0/r) * np.einsum('il,jm,kn,lmn->ijk', F, F_inv, F_inv, EPS3)

    t12 = (1.0+eps)**(1.0/3.0)*(1.0+s)
    t20 = (1.0+eps)**(1.0/3.0)/(1.0+s)
    t01 = (1.0+eps)**(-2.0/3.0)
    C4[3,1,2]=+4.*L*omega1*t12/r**2; C4[3,2,1]=-4.*L*omega1*t12/r**2
    C4[3,2,0]=+4.*L*omega2*t20/r**2; C4[3,0,2]=-4.*L*omega2*t20/r**2
    C4[3,0,1]=+4.*L*omega3*t01/r**2; C4[3,1,0]=-4.*L*omega3*t01/r**2

    if not (sd0==0. and sd1==0. and sd2==0.):
        M  = _make_mixing_matrix(sd0, sd1, sd2)
        C4 = np.einsum('ad,eb,fc,def->abc', M, M, M, C4)

    T4 = np.zeros((4, 4, 4))
    T4[:3,:3,:3] = (2.0*eta/r) * EPS3
    if V_torsion != 0.0:
        v3 = V_torsion / 3.0
        for a in range(3):
            T4[a,3,a] += v3; T4[a,a,3] -= v3

    Gamma_LC = 0.5*(C4
                    + np.transpose(C4,(2,1,0))
                    - np.transpose(C4,(1,0,2)))
    K = 0.5*(T4
             + np.transpose(T4,(2,0,1))
             - np.transpose(T4,(1,2,0)))
    Gamma = Gamma_LC + K

    Riem = (np.einsum('aec,ebd->abcd', Gamma, Gamma)
            - np.einsum('aed,ebc->abcd', Gamma, Gamma)
            + np.einsum('abe,ecd->abcd', Gamma, C4))

    R_scalar = np.einsum('abab->', Riem)

    N = 0.0
    if theta_NY != 0.0:
        N_TT  = 0.25*np.einsum('abcd,eab,ecd->', EPS4, T4, T4)
        N_Ree = 0.25*np.einsum('abcd,abcd->', EPS4, Riem)
        N = N_TT - N_Ree

    Vol = 2.0*np.pi**2*L*r**3
    return -(R_scalar/(2.0*kappa**2) + theta_NY*N)*Vol

# =============================================================================
# FD helpers (background fixed per call)
# =============================================================================

def hess_diag_at(p_key: str, p0: float, bg: dict, h: float = H) -> float:
    def v(**kw): return compute_Veff(**{**bg, **kw})
    return (v(**{p_key: p0+h}) - 2.*v(**{p_key: p0}) + v(**{p_key: p0-h})) / h**2


def hess_cross_at(p1: str, p1_0: float, p2: str, p2_0: float,
                  bg: dict, h: float = H) -> float:
    def v(**kw): return compute_Veff(**{**bg, **kw})
    vpp = v(**{p1: p1_0+h, p2: p2_0+h})
    vpm = v(**{p1: p1_0+h, p2: p2_0-h})
    vmp = v(**{p1: p1_0-h, p2: p2_0+h})
    vmm = v(**{p1: p1_0-h, p2: p2_0-h})
    return (vpp - vpm - vmp + vmm) / (4.*h**2)


def analyze_block(label: str, H_rr: float, H_reta: float, H_etaeta: float) -> dict:
    """Compute det, trace, eigenvalues, eigenvectors, stability verdict."""
    mat  = np.array([[H_rr, H_reta], [H_reta, H_etaeta]])
    eigvals, eigvecs = np.linalg.eigh(mat)
    det  = H_rr * H_etaeta - H_reta**2
    tr   = H_rr + H_etaeta
    if eigvals[0] > 0:
        verdict = "positive definite (stable)"
    elif eigvals[0] >= -1e-3 * abs(eigvals[1]):
        verdict = "positive semidefinite (near-marginal)"
    elif eigvals[1] < 0:
        verdict = "negative definite (unstable)"
    else:
        verdict = "indefinite (saddle)"
    return dict(label=label, H_rr=H_rr, H_reta=H_reta, H_etaeta=H_etaeta,
                det=det, tr=tr, eig1=eigvals[0], eig2=eigvals[1],
                vec1=eigvecs[:,0], vec2=eigvecs[:,1], verdict=verdict)


# =============================================================================
# Main
# =============================================================================

def main():
    t_start = _time.time()
    SEP = "=" * 72
    sep = "-" * 72

    print(SEP)
    print("Coupled (r, eta) 2x2 Hessian block — eigenvalue analysis")
    print("Instruction: 10_logs/20260318-05_r-eta_kotyu-chi_hyoka.md")
    print(f"Background: r0={R0}, L={L_VAL}, kappa={KAPPA}, eta={ETA_BG}, V={V_BG}")
    print(f"            theta_NY={THETA_NY},  h={H}")
    print(SEP)

    # ── Sanity check ───────────────────────────────────────────────────────────
    print("\n── Sanity check (App E.6 / cross-Hessian script reference values) ──")
    bg_mx = dict(r=R0, eps=0., s=0., q3=0., q4=0., q5=0.,
                 omega1=0., omega2=0., omega3=0.,
                 sd0=0., sd1=0., sd2=0.,
                 L=L_VAL, kappa=KAPPA,
                 eta=ETA_BG, V_torsion=V_BG, theta_NY=THETA_NY)

    eta_ref  = hess_diag_at('eta', ETA_BG, bg_mx)
    reta_ref = hess_cross_at('r', R0, 'eta', ETA_BG, bg_mx)
    print(f"  H_etaeta (reference) = {eta_ref:.4f}   [expect ~355.31]")
    print(f"  H_reta   (reference) = {reta_ref:.4f}   [expect ~118.44]")
    ok = abs(eta_ref - 355.31) < 1.0 and abs(reta_ref - 118.44) < 1.0
    print(f"  Sanity check: {'OK ✓' if ok else 'NG ✗ — check pipeline'}")

    # =========================================================================
    # Mode comparison: AX / VT / MX
    # =========================================================================
    print()
    print(SEP)
    print("Mode comparison: AX / VT / MX")
    print("(r-eta mixing arises from axial torsion; radial-V coupling expected zero)")
    print(SEP)

    mode_configs = [
        # label,  eta,    V_torsion
        ("AX (eta=1, V=0)",  ETA_BG, 0.0),
        ("VT (eta=0, V=1)",  0.0,    1.0),
        ("MX (eta=1, V=1)",  ETA_BG, 1.0),
    ]

    results = []
    for mode_label, eta_val, V_val in mode_configs:
        bg = dict(r=R0, eps=0., s=0., q3=0., q4=0., q5=0.,
                  omega1=0., omega2=0., omega3=0.,
                  sd0=0., sd1=0., sd2=0.,
                  L=L_VAL, kappa=KAPPA,
                  eta=eta_val, V_torsion=V_val, theta_NY=THETA_NY)

        # For VT mode, eta=0 but we still differentiate w.r.t. eta around 0
        eta_diff_bg = eta_val  # background value for differentiation

        H_rr      = hess_diag_at('r',   R0,           bg)
        H_reta    = hess_cross_at('r', R0, 'eta', eta_diff_bg, bg)
        H_etaeta  = hess_diag_at('eta', eta_diff_bg,  bg)

        res = analyze_block(mode_label, H_rr, H_reta, H_etaeta)
        results.append(res)

    # ── Primary result: MX (main case) ────────────────────────────────────────
    main_res = results[0]   # AX corresponds to paper03 main case (eta=1, V=0)

    print()
    print(sep)
    print("Primary result: AX mode  (eta=1, V=0)  — paper03 main case")
    print(sep)
    print(f"  H_rr      = {main_res['H_rr']:.6f}")
    print(f"  H_reta    = {main_res['H_reta']:.6f}")
    print(f"  H_etaeta  = {main_res['H_etaeta']:.6f}")
    print(f"  det H     = {main_res['det']:.6e}")
    print(f"  tr  H     = {main_res['tr']:.6f}")
    print(f"  eigenvalue 1 (lambda_1) = {main_res['eig1']:.6f}")
    print(f"  eigenvalue 2 (lambda_2) = {main_res['eig2']:.6f}")
    print(f"  eigenvector 1 = [{main_res['vec1'][0]:.4f},  {main_res['vec1'][1]:.4f}]  (r, eta)")
    print(f"  eigenvector 2 = [{main_res['vec2'][0]:.4f},  {main_res['vec2'][1]:.4f}]  (r, eta)")
    print(f"  verdict: {main_res['verdict']}")

    # ── Numerical table ────────────────────────────────────────────────────────
    print()
    print(sep)
    print("H^(r,eta) numerical table  [AX mode, primary result]")
    print(sep)
    rows = [
        ("H_rr",      main_res['H_rr']),
        ("H_reta",    main_res['H_reta']),
        ("H_etaeta",  main_res['H_etaeta']),
        ("det H",     main_res['det']),
        ("tr  H",     main_res['tr']),
        ("lambda_1",  main_res['eig1']),
        ("lambda_2",  main_res['eig2']),
    ]
    print(f"  {'component':>12}  {'value':>18}")
    print("  " + "-" * 34)
    for label, val in rows:
        print(f"  {label:>12}  {val:>18.6f}")

    # ── Mode comparison table ──────────────────────────────────────────────────
    print()
    print(sep)
    print("Mode comparison table")
    print(sep)
    print(f"  {'mode':<22}  {'H_rr':>12}  {'H_reta':>10}  {'H_etaeta':>10}"
          f"  {'lambda_1':>10}  {'lambda_2':>10}  verdict")
    print("  " + "-" * 90)
    for res in results:
        print(f"  {res['label']:<22}  {res['H_rr']:>12.4f}  {res['H_reta']:>10.4f}"
              f"  {res['H_etaeta']:>10.4f}  {res['eig1']:>10.4f}  {res['eig2']:>10.4f}"
              f"  {res['verdict']}")

    # ── r-eta mixing origin summary ───────────────────────────────────────────
    ax_reta  = results[0]['H_reta']
    vt_reta  = results[1]['H_reta']
    mx_reta  = results[2]['H_reta']
    print()
    print(f"  r-eta cross-Hessian by mode:")
    print(f"    AX (eta=1, V=0): H_reta = {ax_reta:.4f}")
    print(f"    VT (eta=0, V=1): H_reta = {vt_reta:.4f}")
    print(f"    MX (eta=1, V=1): H_reta = {mx_reta:.4f}")
    if abs(ax_reta) > 1.0 and abs(vt_reta) < 1.0:
        print("  -> r-eta mixing is specific to the axial torsion component (AX/MX).")
    elif abs(ax_reta) < 1.0 and abs(vt_reta) < 1.0:
        print("  -> r-eta mixing vanishes in both AX and VT; likely decoupled in this background.")
    else:
        print("  -> r-eta mixing present in both AX and VT modes.")

    # ── Grand verdict ─────────────────────────────────────────────────────────
    print()
    print(SEP)
    print("Grand verdict")
    print(SEP)
    print(f"  Coupled (r, eta) 2x2 block [AX mode, paper03 main case]:")
    print(f"    lambda_1 = {main_res['eig1']:.4f}")
    print(f"    lambda_2 = {main_res['eig2']:.4f}")
    print(f"    {main_res['verdict']}")
    print()
    if "positive definite" in main_res['verdict']:
        print("  -> Scalar/torsion subsector is locally stable even with r-eta mixing.")
        print("     Paper03 §6.5 / Appendix H can state:")
        print("     'The coupled (r,eta) block has positive-definite Hessian")
        print(f"      (eigenvalues approximately {main_res['eig1']:.1f} and {main_res['eig2']:.1f}),")
        print("      confirming local stability of the scalar/torsion subsector.'")
    elif "indefinite" in main_res['verdict'] or "semidefinite" in main_res['verdict']:
        print("  -> Near-marginal or saddle direction found.")
        print("     Investigate eigenvector direction and physical interpretation.")
    else:
        print("  -> Instability detected. Re-check background and parameter choice.")

    # ── Appendix H draft ──────────────────────────────────────────────────────
    print()
    print(SEP)
    print("Appendix H draft  (H^(r,eta) block supplement)")
    print(SEP)
    print()
    print("  At the isotropic MX point, the coupled (r, eta) 2x2 Hessian block is")
    print("  (numerically):")
    print()
    print(f"    H^(r,eta) = [[{main_res['H_rr']:.2f},  {main_res['H_reta']:.2f}],")
    print(f"                  [{main_res['H_reta']:.2f},  {main_res['H_etaeta']:.2f}]]")
    print()
    print(f"  eigenvalues: lambda_1 = {main_res['eig1']:.2f},  lambda_2 = {main_res['eig2']:.2f}")
    print(f"  det = {main_res['det']:.4e},  tr = {main_res['tr']:.4f}")
    print(f"  -> {main_res['verdict']}")

    elapsed = _time.time() - t_start
    print()
    print(f"Computation time: {elapsed:.1f}s")
    print()
    print("*Script: scripts/paper03/spin0_reta_block_eigenvalues.py")


if __name__ == '__main__':
    main()
