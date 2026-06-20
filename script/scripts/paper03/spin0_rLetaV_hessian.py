"""
4x4 Hessian of V_eff over (r, L, eta, V) — extended scalar sector analysis.

Instruction: 10_logs/20260319-01_spin-0 拡張セクターヘッシアンの確認評価.md

Purpose:
  Previous checks showed the (r, eta) 2x2 block is saddle-type (H_rr ~ 0 in
  AX background).  This script temporarily releases the S^1 scale L and
  evaluates the full 4x4 Hessian H(r, L, eta, V) to determine:
    (1) Is L a spectator, or does it mix significantly with r/eta/V?
    (2) How do the known (r,eta,V) cross terms change?
    (3) Does including L alter the local stability verdict per mode (AX/VT/MX)?

Background: r0=3, L0=1, kappa=1, alpha=theta_NY=0 (pure EC gravity)
Modes evaluated: AX (eta=1, V=0), VT (eta=0, V=1), MX (eta=1, V=1)
Method: central FD with multiple step sizes h=1e-2, 1e-3, 5e-4, 1e-4

Author: Muacca
Date: 2026-03-19
"""

from __future__ import annotations

import io
import sys
import time as _time
from itertools import permutations

import numpy as np

# ── UTF-8 stdout ──────────────────────────────────────────────────────────────
def _configure_utf8_stdio() -> None:
    for name in ("stdout", "stderr"):
        stream = getattr(sys, name)
        try:
            stream.reconfigure(encoding="utf-8", errors="replace"); continue
        except Exception:
            pass
        buf = getattr(stream, "buffer", None)
        if buf is not None:
            setattr(sys, name, io.TextIOWrapper(buf, encoding="utf-8", errors="replace"))

_configure_utf8_stdio()

# =============================================================================
# Physical parameters (background)
# =============================================================================
R0        = 3.0
L_VAL     = 1.0
KAPPA     = 1.0
ETA_BG    = 1.0
V_BG      = 0.0
THETA_NY  = 0.0

H_STEPS   = [1e-2, 1e-3, 5e-4, 1e-4]   # step sizes to compare
H_PRIMARY = 1e-3                         # primary step for tables

# =============================================================================
# Levi-Civita symbols
# =============================================================================
def _make_eps3() -> np.ndarray:
    e = np.zeros((3, 3, 3))
    e[0,1,2]=e[1,2,0]=e[2,0,1]=+1.; e[0,2,1]=e[2,1,0]=e[1,0,2]=-1.
    return e

def _make_eps4() -> np.ndarray:
    e = np.zeros((4, 4, 4, 4))
    for p in permutations(range(4)):
        inv = sum(1 for i in range(4) for j in range(i+1,4) if p[i]>p[j])
        e[p] = +1. if inv%2==0 else -1.
    return e

EPS3 = _make_eps3()
EPS4 = _make_eps4()

# =============================================================================
# V_eff  (identical to spin0_torsion_cross_hessian.py)
# =============================================================================
def compute_Veff(
    r: float, eps: float = 0., s: float = 0.,
    q3: float = 0., q4: float = 0., q5: float = 0.,
    omega1: float = 0., omega2: float = 0., omega3: float = 0.,
    sd0: float = 0., sd1: float = 0., sd2: float = 0.,
    L: float = L_VAL, kappa: float = KAPPA,
    eta: float = ETA_BG, V_torsion: float = V_BG,
    theta_NY: float = THETA_NY,
) -> float:
    f0 = (1.+eps)**(2./3.)*(1.+s)**2
    f1 = (1.+eps)**(2./3.)/(1.+s)**2
    f2 = (1.+eps)**(-4./3.)
    sq2 = np.sqrt(2.)
    c3,sh3=np.cosh(q3/sq2),np.sinh(q3/sq2)
    c4,sh4=np.cosh(q4/sq2),np.sinh(q4/sq2)
    c5,sh5=np.cosh(q5/sq2),np.sinh(q5/sq2)
    G3=np.array([[c3,sh3,0.],[sh3,c3,0.],[0.,0.,1.]])
    G4=np.array([[c4,0.,sh4],[0.,1.,0.],[sh4,0.,c4]])
    G5=np.array([[1.,0.,0.],[0.,c5,sh5],[0.,sh5,c5]])
    G = G3@G4@G5
    F     = np.diag([f0,f1,f2])@G
    F_inv = np.linalg.inv(F)

    C4 = np.zeros((4,4,4))
    C4[:3,:3,:3] = (4./r)*np.einsum('il,jm,kn,lmn->ijk',F,F_inv,F_inv,EPS3)
    t12=(1.+eps)**(1./3.)*(1.+s)
    t20=(1.+eps)**(1./3.)/(1.+s)
    t01=(1.+eps)**(-2./3.)
    C4[3,1,2]=+4.*L*omega1*t12/r**2; C4[3,2,1]=-4.*L*omega1*t12/r**2
    C4[3,2,0]=+4.*L*omega2*t20/r**2; C4[3,0,2]=-4.*L*omega2*t20/r**2
    C4[3,0,1]=+4.*L*omega3*t01/r**2; C4[3,1,0]=-4.*L*omega3*t01/r**2

    T4 = np.zeros((4,4,4))
    T4[:3,:3,:3] = (2.*eta/r)*EPS3
    if V_torsion != 0.:
        v3 = V_torsion/3.
        for a in range(3): T4[a,3,a]+=v3; T4[a,a,3]-=v3

    Gamma_LC = 0.5*(C4 + np.transpose(C4,(2,1,0)) - np.transpose(C4,(1,0,2)))
    K        = 0.5*(T4 + np.transpose(T4,(2,0,1)) - np.transpose(T4,(1,2,0)))
    Gamma    = Gamma_LC + K

    Riem = (np.einsum('aec,ebd->abcd',Gamma,Gamma)
            -np.einsum('aed,ebc->abcd',Gamma,Gamma)
            +np.einsum('abe,ecd->abcd',Gamma,C4))
    R_scalar = np.einsum('abab->',Riem)

    N = 0.
    if theta_NY != 0.:
        N_TT  = 0.25*np.einsum('abcd,eab,ecd->',EPS4,T4,T4)
        N_Ree = 0.25*np.einsum('abcd,abcd->',EPS4,Riem)
        N = N_TT - N_Ree

    Vol = 2.*np.pi**2*L*r**3
    return -(R_scalar/(2.*kappa**2) + theta_NY*N)*Vol

# =============================================================================
# Finite-difference helpers (background dict)
# =============================================================================
VARS = ["r", "L", "eta", "V_torsion"]    # canonical order for 4x4 matrix
VAR_LABELS = ["r", "L", "eta", "V"]      # display labels

def _veff(bg: dict, **overrides) -> float:
    return compute_Veff(**{**bg, **overrides})

def hdiag(key: str, p0: float, bg: dict, h: float) -> float:
    return (_veff(bg, **{key: p0+h}) - 2.*_veff(bg, **{key: p0}) + _veff(bg, **{key: p0-h})) / h**2

def hcross(k1: str, p1: float, k2: str, p2: float, bg: dict, h: float) -> float:
    vpp = _veff(bg, **{k1: p1+h, k2: p2+h})
    vpm = _veff(bg, **{k1: p1+h, k2: p2-h})
    vmp = _veff(bg, **{k1: p1-h, k2: p2+h})
    vmm = _veff(bg, **{k1: p1-h, k2: p2-h})
    return (vpp - vpm - vmp + vmm) / (4.*h**2)

def build_hessian_4x4(bg: dict, bg_vals: dict, h: float) -> np.ndarray:
    """Build the 4x4 Hessian H[i,j] for vars=(r, L, eta, V_torsion)."""
    H = np.zeros((4, 4))
    for i, (ki, pi) in enumerate(zip(VARS, [bg_vals[v] for v in VARS])):
        for j, (kj, pj) in enumerate(zip(VARS, [bg_vals[v] for v in VARS])):
            if i == j:
                H[i, i] = hdiag(ki, pi, bg, h)
            elif j > i:
                H[i, j] = hcross(ki, pi, kj, pj, bg, h)
            else:
                H[i, j] = H[j, i]
    return H

# =============================================================================
# Helpers for verdict / display
# =============================================================================
ZERO_THRESH = 1.0   # absolute threshold: smaller → zero-compatible candidate

def _zero_verdict(val: float, vals_across_h: list[float]) -> str:
    """Classify a Hessian entry as zero-compatible / nonzero / inconclusive."""
    maxabs = max(abs(v) for v in vals_across_h)
    minabs = min(abs(v) for v in vals_across_h)
    stable_sign = all(v > 0 for v in vals_across_h) or all(v < 0 for v in vals_across_h)
    if maxabs < ZERO_THRESH:
        if maxabs / minabs < 5 if minabs > 1e-10 else True:
            return "zero-compatible"
        return "zero-compatible (noisy)"
    if stable_sign and (maxabs / max(minabs, 1e-12)) < 3.:
        return "clearly nonzero"
    if maxabs >= ZERO_THRESH and not stable_sign:
        return "inconclusive"
    return "clearly nonzero"

def _stability(eigvals: np.ndarray) -> str:
    if eigvals[0] > 1e-3:
        return "positive definite (stable)"
    if eigvals[0] > -1e-3 * max(abs(eigvals[-1]), 1.):
        return "positive semidefinite (near-marginal)"
    if eigvals[-1] < 0:
        return "negative definite (unstable)"
    return "indefinite / saddle"

def _print_matrix(M: np.ndarray, labels: list[str], indent: str = "  ") -> None:
    n = len(labels)
    header = indent + f"{'':>8}" + "".join(f"  {lbl:>10}" for lbl in labels)
    print(header)
    print(indent + "-" * (10 + 13*n))
    for i, row_lbl in enumerate(labels):
        row = indent + f"{row_lbl:>8}"
        for j in range(n):
            row += f"  {M[i,j]:>10.4f}"
        print(row)

def _print_eigvals(eigvals: np.ndarray, eigvecs: np.ndarray,
                   labels: list[str], indent: str = "  ") -> None:
    for k, (lam, vec) in enumerate(zip(eigvals, eigvecs.T)):
        dominant = labels[int(np.argmax(np.abs(vec)))]
        comp = ", ".join(f"{labels[i]}:{vec[i]:+.3f}" for i in range(len(labels)))
        print(f"{indent}  lambda_{k+1} = {lam:>10.4f}   dominant: {dominant}   [{comp}]")

# =============================================================================
# Main
# =============================================================================
SEP = "=" * 76
sep = "-" * 76

def main() -> None:
    t0 = _time.time()

    print(SEP)
    print("4x4 Hessian H(r, L, eta, V) — extended scalar sector")
    print("Instruction: 10_logs/20260319-01_spin-0 拡張セクターヘッシアンの確認評価.md")
    print(f"Background: r0={R0}, L0={L_VAL}, kappa={KAPPA},")
    print(f"            theta_NY={THETA_NY},  step sizes: {H_STEPS}")
    print(SEP)

    # ── Mode definitions ────────────────────────────────────────────────────
    mode_configs = [
        ("AX  (eta=1, V=0)", dict(r=R0, L=L_VAL, kappa=KAPPA,
                                   eta=1., V_torsion=0., theta_NY=THETA_NY,
                                   eps=0., s=0., q3=0., q4=0., q5=0.,
                                   omega1=0., omega2=0., omega3=0.,
                                   sd0=0., sd1=0., sd2=0.)),
        ("VT  (eta=0, V=1)", dict(r=R0, L=L_VAL, kappa=KAPPA,
                                   eta=0., V_torsion=1., theta_NY=THETA_NY,
                                   eps=0., s=0., q3=0., q4=0., q5=0.,
                                   omega1=0., omega2=0., omega3=0.,
                                   sd0=0., sd1=0., sd2=0.)),
        ("MX  (eta=1, V=1)", dict(r=R0, L=L_VAL, kappa=KAPPA,
                                   eta=1., V_torsion=1., theta_NY=THETA_NY,
                                   eps=0., s=0., q3=0., q4=0., q5=0.,
                                   omega1=0., omega2=0., omega3=0.,
                                   sd0=0., sd1=0., sd2=0.)),
    ]

    # ── Sanity check (AX with H_PRIMARY) ────────────────────────────────────
    print()
    print(sep)
    print("Sanity check — AX mode, h = H_PRIMARY = 1e-3")
    print(sep)
    ax_bg  = mode_configs[0][1]
    ax_bg_vals = {"r": R0, "L": L_VAL, "eta": 1., "V_torsion": 0.}
    H_rr_san   = hdiag("r",   R0,   ax_bg, H_PRIMARY)
    H_ee_san   = hdiag("eta", 1.,   ax_bg, H_PRIMARY)
    H_reta_san = hcross("r", R0, "eta", 1., ax_bg, H_PRIMARY)
    H_VV_san   = hdiag("V_torsion", 0., ax_bg, H_PRIMARY)
    H_rV_san   = hcross("r", R0, "V_torsion", 0., ax_bg, H_PRIMARY)
    H_eV_san   = hcross("eta", 1., "V_torsion", 0., ax_bg, H_PRIMARY)
    print(f"  H_rr       = {H_rr_san:>12.6f}   [expect ~0]")
    print(f"  H_etaeta   = {H_ee_san:>12.6f}   [expect ~355.31]")
    print(f"  H_reta     = {H_reta_san:>12.6f}   [expect ~118.44]")
    print(f"  H_VV       = {H_VV_san:>12.6f}   [expect ~355.31]")
    print(f"  H_rV       = {H_rV_san:>12.6f}   [expect ~0]")
    print(f"  H_etaV     = {H_eV_san:>12.6f}   [expect ~0]")
    ok = (abs(H_ee_san - 355.31) < 1. and abs(H_reta_san - 118.44) < 1.
          and abs(H_rr_san) < 1. and abs(H_rV_san) < 1. and abs(H_eV_san) < 1.)
    print(f"  Sanity: {'OK ✓' if ok else 'NG ✗'}")

    # ── Step-size stability for L-related cross terms ────────────────────────
    print()
    print(sep)
    print("Step-size stability: L-related cross terms (AX mode)")
    print(sep)
    cross_pairs = [
        ("r",   R0,  "L",         L_VAL,  "H_rL"),
        ("L",   L_VAL, "eta",     1.,     "H_Leta"),
        ("L",   L_VAL, "V_torsion", 0.,   "H_LV"),
        ("L",   L_VAL, "L",       L_VAL,  "H_LL"),   # diagonal
    ]
    step_table: dict[str, list[float]] = {pair[-1]: [] for pair in cross_pairs}
    print(f"  {'term':<10}" + "".join(f"  {h:>12.0e}" for h in H_STEPS))
    print("  " + "-" * (12 + 14*len(H_STEPS)))
    for k1, p1, k2, p2, name in cross_pairs:
        vals = []
        for h in H_STEPS:
            if k1 == k2:
                v = hdiag(k1, p1, ax_bg, h)
            else:
                v = hcross(k1, p1, k2, p2, ax_bg, h)
            vals.append(v)
        step_table[name] = vals
        row = f"  {name:<10}" + "".join(f"  {v:>12.4f}" for v in vals)
        print(row)

    # verdict based on multi-h stability
    print()
    print("  Zero verdicts (threshold |val| < 1.0):")
    for name, vals in step_table.items():
        verdict = _zero_verdict(vals[-1], vals)   # primary = last h value
        print(f"    {name:<10}  {verdict}")

    # ── Full 4x4 for each mode ───────────────────────────────────────────────
    mode_results: list[tuple[str, np.ndarray, np.ndarray, np.ndarray]] = []

    for mode_label, bg in mode_configs:
        bg_vals = {"r": R0, "L": L_VAL,
                   "eta": bg["eta"], "V_torsion": bg["V_torsion"]}

        print()
        print(SEP)
        print(f"Mode: {mode_label}")
        print(SEP)

        # --- build primary and symmetrized matrix ---------------------------
        H4 = build_hessian_4x4(bg, bg_vals, H_PRIMARY)
        H4_sym = 0.5 * (H4 + H4.T)

        # symmetry check
        asym_max = np.max(np.abs(H4 - H4.T))
        print(f"\n  Asymmetry max |H_ab - H_ba| = {asym_max:.2e}")

        # --- print full matrix -----------------------------------------------
        print(f"\n  4x4 Hessian H(r, L, eta, V)  [h={H_PRIMARY:.0e}, symmetrized]")
        _print_matrix(H4_sym, VAR_LABELS)

        # --- cross-term verdict ----------------------------------------------
        print()
        print("  Cross-term verdicts  (threshold |val| < 1.0):")
        cross_entries = [
            ("H_rL",  "r", R0, "L", L_VAL),
            ("H_rV",  "r", R0, "V_torsion", bg["V_torsion"]),
            ("H_Leta","L", L_VAL, "eta", bg["eta"]),
            ("H_LV",  "L", L_VAL, "V_torsion", bg["V_torsion"]),
            ("H_reta","r", R0, "eta", bg["eta"]),
            ("H_etaV","eta", bg["eta"], "V_torsion", bg["V_torsion"]),
        ]
        for name, k1, p1, k2, p2 in cross_entries:
            vals = [hcross(k1,p1,k2,p2,bg,h) for h in H_STEPS]
            verdict = _zero_verdict(vals[1], vals)
            print(f"    {name:<10}  h-values: " +
                  ", ".join(f"{v:.4f}" for v in vals) +
                  f"   -> {verdict}")

        # --- eigenvalues ----------------------------------------------------
        eigvals, eigvecs = np.linalg.eigh(H4_sym)
        stability = _stability(eigvals)
        print()
        print(f"  Eigenvalues and eigenvectors:")
        _print_eigvals(eigvals, eigvecs, VAR_LABELS)
        print(f"  Stability verdict: {stability}")

        mode_results.append((mode_label, H4_sym, eigvals, eigvecs))

    # ── Partial block analysis (AX mode) ────────────────────────────────────
    ax_label, ax_H, ax_eig, ax_vec = mode_results[0]
    # index map: r=0, L=1, eta=2, V=3
    print()
    print(SEP)
    print("Partial block analysis — AX mode")
    print(SEP)

    block_defs = [
        ("(r, eta)",     [0, 2]),
        ("(r, L)",       [0, 1]),
        ("(L, V)",       [1, 3]),
        ("(r, eta, V)",  [0, 2, 3]),
        ("(r, L, V)",    [0, 1, 3]),
    ]
    for blk_name, idx in block_defs:
        blk = ax_H[np.ix_(idx, idx)]
        blk_labels = [VAR_LABELS[i] for i in idx]
        eig_blk = np.linalg.eigvalsh(blk)
        print(f"\n  Block {blk_name}:")
        _print_matrix(blk, blk_labels, indent="    ")
        print(f"    eigenvalues: " + ", ".join(f"{e:.4f}" for e in eig_blk))
        print(f"    stability:   {_stability(eig_blk)}")

    # ── Mode comparison summary ─────────────────────────────────────────────
    print()
    print(SEP)
    print("Mode comparison summary")
    print(SEP)
    print(f"  {'mode':<22}  {'lambda_min':>12}  {'lambda_max':>12}  {'H_rL':>10}  {'H_LL':>10}  verdict")
    print("  " + "-" * 80)
    for mode_label, H4s, eigvals, _ in mode_results:
        bg = [b for lbl, b in mode_configs if lbl == mode_label][0]
        H_LL  = hdiag("L", L_VAL, bg, H_PRIMARY)
        H_rL  = hcross("r", R0, "L", L_VAL, bg, H_PRIMARY)
        print(f"  {mode_label:<22}  {eigvals[0]:>12.4f}  {eigvals[-1]:>12.4f}"
              f"  {H_rL:>10.4f}  {H_LL:>10.4f}  {_stability(eigvals)}")

    # ── Physical summary ────────────────────────────────────────────────────
    print()
    print(SEP)
    print("Physical summary")
    print(SEP)
    _, ax_H4, ax_eig4, ax_vec4 = mode_results[0]
    H_LL   = ax_H4[1, 1]
    H_rL   = ax_H4[0, 1]
    H_Leta = ax_H4[1, 2]
    H_LV   = ax_H4[1, 3]

    def _spectator_verdict(H_rL, H_Leta, H_LV):
        cross_max = max(abs(H_rL), abs(H_Leta), abs(H_LV))
        H_diag_mag = 355.  # ~ H_etaeta for reference
        ratio = cross_max / H_diag_mag
        if ratio < 0.01:
            return "L is spectator (couplings < 1% of diagonal scale)"
        if ratio < 0.1:
            return "L mixes weakly but measurably"
        return "L mixes significantly"

    print(f"  H_LL   = {H_LL:.4f}")
    print(f"  H_rL   = {H_rL:.4f}")
    print(f"  H_Leta = {H_Leta:.4f}")
    print(f"  H_LV   = {H_LV:.4f}")
    print(f"  -> {_spectator_verdict(H_rL, H_Leta, H_LV)}")
    print()
    print(f"  AX mode 4x4 eigenvalues: " +
          ", ".join(f"{e:.4f}" for e in ax_eig4))
    print(f"  AX mode stability: {_stability(ax_eig4)}")

    print()
    print(f"Computation time: {_time.time()-t0:.1f}s")
    print()
    print("*Script: scripts/paper03/spin0_rLetaV_hessian.py")


if __name__ == "__main__":
    import os as _os, sys as _sys
    _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..'))
    from dppu.utils.tee_logger import setup_log, teardown_log
    setup_log(__file__, log_dir='output')
    main()
    teardown_log()
