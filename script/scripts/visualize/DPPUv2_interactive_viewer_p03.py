# =============================================================================
# DPPUv2 Paper 03 Interactive Phase Diagram, Potential & Spin-2 Spectrum Viewer
# =============================================================================
# Phase diagram in (epsilon, eta) space
#   — extends paper02's (V, eta) diagram to the squash dimension.
#
# Features:
#   - Phase diagram: epsilon (squash) x-axis, eta (axial torsion) y-axis
#   - Fixed sliders: V (vector torsion), theta_NY, alpha (Weyl)
#   - V_eff(r) plot at up to 3 click-selected points
#   - Spin-2 H_eps(epsilon) spectrum at selected eta values
#   - Same color scheme as paper02 (Type I / II / III / G)
#   - Engine: S3xS1 with squash epsilon, Mode.MX (eta and V both active)
#   - shear (s) and twist (omega) excluded — see notebook Markdown cell
#
# Requirements: %matplotlib widget, ipywidgets, scipy, sympy, dppu package
#
# Usage (from notebook):
#   vfunc, info = DPPUv2Paper03Viewer.setup_engine(engine_root)
#   viewer = DPPUv2Paper03Viewer(vfunc)
#   viewer.display()
#
# =============================================================================

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
from scipy.optimize import minimize_scalar
from itertools import permutations


# =============================================================================
# Pure-numpy V_eff for spin-2 branch computation
# (identical physics to the lambdified vfunc; allows s/q3/q4/q5 perturbation)
# =============================================================================

def _make_eps3_():
    e = np.zeros((3, 3, 3))
    e[0,1,2]=e[1,2,0]=e[2,0,1]=+1.; e[0,2,1]=e[2,1,0]=e[1,0,2]=-1.
    return e

def _make_eps4_():
    e = np.zeros((4, 4, 4, 4))
    for p in permutations(range(4)):
        inv = sum(1 for i in range(4) for j in range(i+1,4) if p[i]>p[j])
        e[p] = +1. if inv%2==0 else -1.
    return e

_EPS3_NP = _make_eps3_()
_EPS4_NP = _make_eps4_()


def _compute_Veff_np(r, eps=0., s=0., q3=0., q4=0., q5=0.,
                     eta=1., V_torsion=0., theta_NY=0., alpha=0., L=1., kappa=1.):
    """Pure-numpy V_eff for S3xS1 EC+torsion (Mode.MX, squash+shear+mixing).

    Lagrangian: L = R/(2κ²) + θ_NY·N + α·C²
    V_eff = -L × Vol

    Frame vierbeins (not metric components):
      f0 = (1+eps)^(1/3) * (1+s),   f1 = (1+eps)^(1/3) / (1+s),   f2 = (1+eps)^(-2/3)
    These match the engine's spatial_factors:
      C^i_{jk} = (4/r) * f_i/(f_j*f_k) * eps_{ijk}
              => factor_0 = (1+eps)^(2/3)*(1+s)^2   [= f0/(f1*f2)]  ✓

    Frame metric is flat (η_{ab} = δ_{ab}), so R_{abcd} = Riem[a,b,c,d].
    Weyl tensor: C_{abcd} = R_{abcd}
                  - (1/2)(δ_{ac}R_{bd} - δ_{ad}R_{bc} - δ_{bc}R_{ad} + δ_{bd}R_{ac})
                  + (R/6)(δ_{ac}δ_{bd} - δ_{ad}δ_{bc})
    C² = C_{abcd}C^{abcd} = Σ C_{abcd}² (flat metric)
    """
    f0 = (1.+eps)**(1./3.)*(1.+s)
    f1 = (1.+eps)**(1./3.)/(1.+s)
    f2 = (1.+eps)**(-2./3.)
    sq2 = np.sqrt(2.)
    c3,sh3 = np.cosh(q3/sq2), np.sinh(q3/sq2)
    c4,sh4 = np.cosh(q4/sq2), np.sinh(q4/sq2)
    c5,sh5 = np.cosh(q5/sq2), np.sinh(q5/sq2)
    G3=np.array([[c3,sh3,0.],[sh3,c3,0.],[0.,0.,1.]])
    G4=np.array([[c4,0.,sh4],[0.,1.,0.],[sh4,0.,c4]])
    G5=np.array([[1.,0.,0.],[0.,c5,sh5],[0.,sh5,c5]])
    G = G3@G4@G5
    F     = np.diag([f0,f1,f2])@G
    F_inv = np.linalg.inv(F)

    C4 = np.zeros((4,4,4))
    C4[:3,:3,:3] = (4./r)*np.einsum('il,jm,kn,lmn->ijk',F,F_inv,F_inv,_EPS3_NP)

    T4 = np.zeros((4,4,4))
    T4[:3,:3,:3] = (2.*eta/r)*_EPS3_NP
    if V_torsion != 0.:
        v3 = V_torsion/3.
        for a in range(3): T4[a,3,a]+=v3; T4[a,a,3]-=v3

    Gamma_LC = 0.5*(C4+np.transpose(C4,(2,1,0))-np.transpose(C4,(1,0,2)))
    # Contortion: K_{abc} = (1/2)(T_{abc} + T_{bca} - T_{cab})
    # np.transpose(T4,(1,2,0))[a,b,c] = T4[b,c,a]  = T_{bca}
    # np.transpose(T4,(2,0,1))[a,b,c] = T4[c,a,b]  = T_{cab}
    K        = 0.5*(T4+np.transpose(T4,(1,2,0))-np.transpose(T4,(2,0,1)))
    Gamma    = Gamma_LC + K

    Riem = (np.einsum('aec,ebd->abcd',Gamma,Gamma)
            -np.einsum('aed,ebc->abcd',Gamma,Gamma)
            +np.einsum('abe,ecd->abcd',Gamma,C4))
    R_scalar = np.einsum('abab->',Riem)

    N = 0.
    if theta_NY != 0.:
        N_TT  = 0.25*np.einsum('abcd,eab,ecd->',_EPS4_NP,T4,T4)
        N_Ree = 0.25*np.einsum('abcd,abcd->',_EPS4_NP,Riem)
        N = N_TT - N_Ree

    # Weyl scalar C² = C_{abcd}C^{abcd}
    # Frame metric is flat (δ_{ab}), so C^{abcd} = C_{abcd} and C² = Σ C_{abcd}².
    # Ricci tensor: R_{bd} = Σ_a R_{abad} = einsum('abac->bc', Riem)
    C_sq = 0.
    if alpha != 0.:
        I4     = np.eye(4)
        Ricci  = np.einsum('abac->bc', Riem)
        C_abcd = (Riem
                  - 0.5 * (  np.einsum('ac,bd->abcd', I4, Ricci)
                            - np.einsum('ad,bc->abcd', I4, Ricci)
                            - np.einsum('bc,ad->abcd', I4, Ricci)
                            + np.einsum('bd,ac->abcd', I4, Ricci))
                  + (R_scalar / 6.) * (  np.einsum('ac,bd->abcd', I4, I4)
                                       - np.einsum('ad,bc->abcd', I4, I4)))
        C_sq = np.einsum('abcd,abcd->', C_abcd, C_abcd)

    Vol = 2.*np.pi**2*L*r**3
    return -(R_scalar/(2.*kappa**2) + theta_NY*N + alpha*C_sq)*Vol


def _d2V_np_branch(key, val0, bg_np, h=5e-3):
    """d²V/d(key)² at val0 using pure-numpy Veff."""
    def v(x): return _compute_Veff_np(**{**bg_np, key: x})
    return (v(val0+h) - 2.*v(val0) + v(val0-h)) / h**2


# =============================================================================
# Stability helpers
# =============================================================================

_ALPHA_THRESHOLD       = 1e-9
_COUNT_MINIMA_R_LOW    = 0.02
_COUNT_MINIMA_R_HIGH   = 1e5
_COUNT_MINIMA_N        = 300


def _count_local_minima(v_func, r_low=_COUNT_MINIMA_R_LOW,
                        r_high=_COUNT_MINIMA_R_HIGH, n=_COUNT_MINIMA_N) -> int:
    """
    Count local minima of v_func(r) in [r_low, r_high] using log-spaced sampling.
    A point r_k is a local minimum when v(r_k) < v(r_{k-1}) AND v(r_k) < v(r_{k+1}).
    Used for the α ≠ 0 classification branches (Type G / I-W / II-W).
    """
    r_arr = np.logspace(np.log10(r_low), np.log10(r_high), n)
    vals  = np.empty(n)
    for k, r in enumerate(r_arr):
        try:
            vals[k] = float(v_func(r))
        except Exception:
            vals[k] = np.nan

    count = 0
    for k in range(1, n - 1):
        if (np.isfinite(vals[k - 1]) and np.isfinite(vals[k])
                and np.isfinite(vals[k + 1])
                and vals[k] < vals[k - 1] and vals[k] < vals[k + 1]):
            count += 1
    return count


def _classify_veff(f, r_lo=0.05, r_hi=100.0, alpha=0.0, epsilon=0.0):
    """
    Classify V_eff(r) = f(r) into stability types (matches paper02 / v4 scheme).

    Returns
    -------
    code : int
        α ≈ 0 (standard):
          1 = Type I    (stable well — barrier near r=0)
          2 = Type II   (rolling — no barrier near r=0)
          3 = Type III  (no interior local minimum)
        α > 0 (ghost regime):
          4 = Type G    (has local min — metastable)
          7 = Type III-G (no local min — unbounded below at r→0)
        α < 0 (repulsive-wall regime):
          6 = Type I-W  (double well — two or more local minima)
          5 = Type II-W (single well — exactly one local minimum)
          3 = Type III  (no local minimum)
        0 = Error
    r0 : float  (location of (deepest) minimum, or nan)

    Note on effective α:
        S³×S¹ is conformally flat when ε=0 (round sphere) → C²=0 → α has no
        physical effect.  Only for ε≠0 (squashed) does the Weyl term contribute.
        This matches paper02 v4's get_weyl_c2_status() logic.
    """
    NAN = float('nan')

    # ── C² = 0 check: α is physically inactive when ε = 0 (conformally flat)
    effective_alpha = 0.0 if abs(epsilon) < 1e-10 else alpha

    # ── α ≈ 0: standard Type I / II / III logic ───────────────────────────
    if abs(effective_alpha) < _ALPHA_THRESHOLD:
        try:
            res = minimize_scalar(f, bounds=(r_lo, r_hi), method='bounded',
                                  options={'xatol': 1e-4})
            r0 = float(res.x)

            if r0 < r_lo * 1.02 or r0 > r_hi * 0.98:
                return 3, NAN

            V0 = float(f(r0))
            h  = r0 * 1e-3
            d2V = (f(r0 + h) - 2.0 * V0 + f(r0 - h)) / h ** 2
            if d2V <= 0:
                return 3, NAN

            # Slope near r=0 distinguishes I (barrier) vs II (rolling)
            r_test = r_lo * 2
            slope  = f(r_test + r_lo * 0.1) - f(r_test)
            return (1 if slope > 0 else 2), r0

        except Exception:
            return 0, NAN

    # ── α > 0: ghost regime — check for any local minimum ─────────────────
    elif effective_alpha > 0:
        try:
            n_min = _count_local_minima(f)
            if n_min >= 1:
                res = minimize_scalar(f, bounds=(r_lo, r_hi), method='bounded',
                                      options={'xatol': 1e-4})
                r0 = float(res.x) if res.success else NAN
                return 4, r0
            else:
                return 7, NAN
        except Exception:
            return 0, NAN

    # ── α < 0: repulsive-wall regime — count local minima ─────────────────
    else:
        try:
            n_min = _count_local_minima(f)
            if n_min >= 1:
                res = minimize_scalar(f, bounds=(r_lo, r_hi), method='bounded',
                                      options={'xatol': 1e-4})
                r0 = float(res.x) if res.success else NAN
                return (6 if n_min >= 2 else 5), r0
            else:
                return 3, NAN
        except Exception:
            return 0, NAN


def _find_r0(f, r_lo=0.05, r_hi=100.0):
    """Return (r0, V0) minimising f(r), or (nan, nan) on failure."""
    try:
        res = minimize_scalar(f, bounds=(r_lo, r_hi), method='bounded',
                              options={'xatol': 1e-4})
        return float(res.x), float(res.fun)
    except Exception:
        return float('nan'), float('nan')


def _d2V_deps2(vfunc, eps, eta, V_val, theta, alpha, r0, h=5e-4):
    """
    Finite-difference d^2 V_eff / d epsilon^2 at fixed r = r0.

    Proportional to the spin-2 T1 (m=0) Hessian:
        H11^q ≈ (3/2)(1+epsilon)^2 * H_eps_eps
    (see paper03 Sec. 4 for normalization details).
    """
    try:
        Vp = float(vfunc(r0, eps + h, eta, V_val, theta, alpha))
        Vm = float(vfunc(r0, eps - h, eta, V_val, theta, alpha))
        V0 = float(vfunc(r0, eps,     eta, V_val, theta, alpha))
        return (Vp - 2.0 * V0 + Vm) / h ** 2
    except Exception:
        return float('nan')


# =============================================================================
# Main Viewer Class
# =============================================================================

class DPPUv2Paper03Viewer:
    """
    Interactive viewer for DPPUv2 Paper 03.

    Phase diagram in (epsilon, eta) space.
    V_eff(r) and spin-2 H_eps(epsilon) at click-selected points.
    """

    # ── Slider defaults (modelled on paper02 v4) ──────────────────────────────
    SLIDER_THETA_NY_MIN     = 0.0
    SLIDER_THETA_NY_MAX     = 10.0
    SLIDER_THETA_NY_DEFAULT = 1.0

    SLIDER_V_MIN     = 0.0
    SLIDER_V_MAX     = 15.0
    SLIDER_V_DEFAULT = 4.0

    SLIDER_ALPHA_MIN     = -5.0
    SLIDER_ALPHA_MAX     =  3.0
    SLIDER_ALPHA_DEFAULT =  0.0

    # epsilon range sliders
    SLIDER_EPS_MIN_MIN     = -0.5
    SLIDER_EPS_MIN_MAX     =  0.0
    SLIDER_EPS_MIN_DEFAULT = -0.35

    SLIDER_EPS_MAX_MIN     =  0.0
    SLIDER_EPS_MAX_MAX     =  2.0
    SLIDER_EPS_MAX_DEFAULT =  0.60

    # eta range sliders (paper02 defaults)
    SLIDER_ETA_MIN_MIN     = -50.0
    SLIDER_ETA_MIN_MAX     =   0.0
    SLIDER_ETA_MIN_DEFAULT = -15.0

    SLIDER_ETA_MAX_MIN     =  0.0
    SLIDER_ETA_MAX_MAX     = 50.0
    SLIDER_ETA_MAX_DEFAULT =  5.0

    # Potential plot range
    SLIDER_R_MAX_MIN     = 0.01    # 0.01
    SLIDER_R_MAX_MAX     = 100     # 100
    SLIDER_R_MAX_DEFAULT = 10      # 10

    # ── Phase colors (matches paper02 v4) ─────────────────────────────────────
    COLOR_I    = '#90EE90'   # light green  — Type I    (stable well)
    COLOR_II   = '#FFFF99'   # light yellow — Type II   (rolling)
    COLOR_III  = '#E8E8E8'   # light gray   — Type III  (no well)
    COLOR_G    = '#FFCCCC'   # light red    — Type G    (metastable, α>0)
    COLOR_IIIG = '#E8E8E8'   # light gray   — Type III-G (no well, α>0; hatched)
    COLOR_IIW  = '#FFE0B2'   # light orange — Type II-W (single well, α<0)
    COLOR_IW   = '#BBDEFB'   # light blue   — Type I-W  (double well, α<0)
    COLOR_ERR  = '#FFFFFF'   # white        — error / unclassified
    # Integer → colour mapping:
    #  0: err  1: I  2: II  3: III  4: G  5: II-W  6: I-W  7: III-G

    # ── Branch styles ─────────────────────────────────────────────────────────
    BRANCH_LINESTYLES = ['-', '--', ':']   # m=0, m=±2, m=±1

    # ── Grid resolution ───────────────────────────────────────────────────────
    GRID_N_EPS = 40
    GRID_N_ETA = 40

    # =========================================================================
    # Engine setup (class method)
    # =========================================================================

    @classmethod
    def setup_engine(cls, engine_root: str):
        """
        Run S3xS1 squash engine (Mode.MX) and return lambdified V_eff.

        Parameters
        ----------
        engine_root : str
            Path containing the `dppu` package (two levels up from this file).

        Returns
        -------
        vfunc : callable
            f(r, eps, eta, V_param, theta_NY, alpha) -> float
        info : dict
            {'engine': engine, 'V_eff_sym': sympy_expr}
        """
        if engine_root not in sys.path:
            sys.path.insert(0, engine_root)

        from dppu.topology.unified import UnifiedEngine, DOFConfig, TopologyType, FiberMode
        from dppu.torsion.mode import Mode
        from dppu.torsion.nieh_yan import NyVariant
        from dppu.engine.logger import NullLogger
        from sympy import lambdify, S

        print("Initialising S3xS1 squash engine (Mode.MX, NyVariant.FULL) ...")
        cfg = DOFConfig(
            topology=TopologyType.S3,
            enable_squash=True,
            enable_shear=False,
            fiber_mode=FiberMode.NONE,
            torsion_mode=Mode.MX,
            ny_variant=NyVariant.FULL,
        )
        engine = UnifiedEngine(cfg, NullLogger())
        engine.run()
        print("Engine complete. Lambdifying V_eff ...")

        params    = engine.data['params']
        V_eff_sym = engine.data['potential']

        # Fix L=1, kappa=1  (user can adjust via wrapper defaults later)
        V_sub = V_eff_sym.subs({params['L']: S.One, params['kappa']: S.One})

        r_s      = params['r']
        eps_s    = params['epsilon']
        eta_s    = params['eta']
        V_s      = params['V']
        theta_s  = params['theta_NY']
        alpha_s  = params['alpha']

        _raw = lambdify(
            [r_s, eps_s, eta_s, V_s, theta_s, alpha_s],
            V_sub,
            modules='numpy',
            cse=True,
        )

        def vfunc(r, eps, eta, V_val, theta, alpha):
            try:
                v = float(_raw(r, eps, eta, V_val, theta, alpha))
                return v if np.isfinite(v) else 1e20
            except Exception:
                return 1e20

        print("Lambdify done. Viewer is ready.\n")
        return vfunc, {'engine': engine, 'V_eff_sym': V_eff_sym}

    # =========================================================================
    # Constructor
    # =========================================================================

    def __init__(self, vfunc, info=None):
        self.vfunc = vfunc
        self.info  = info or {}

        # Internal state
        self.phase_data = None   # 2D int array of stability codes
        self.Heps_data  = None   # 2D float array of d^2V/deps^2  (m=0)
        self.Hs_data    = None   # 2D float array of d^2V/ds^2 ×0.5 (m=±2)
        self.eps_grid   = None
        self.eta_grid   = None

        self.selected_point = {'eps': 0.0, 'eta': -4.0}

        # Current parameter values (kept in sync with sliders)
        self.theta_NY = self.SLIDER_THETA_NY_DEFAULT
        self.V_val    = self.SLIDER_V_DEFAULT
        self.alpha    = self.SLIDER_ALPHA_DEFAULT
        self.eps_min  = self.SLIDER_EPS_MIN_DEFAULT
        self.eps_max  = self.SLIDER_EPS_MAX_DEFAULT
        self.eta_min  = self.SLIDER_ETA_MIN_DEFAULT
        self.eta_max  = self.SLIDER_ETA_MAX_DEFAULT
        self.r_max    = self.SLIDER_R_MAX_DEFAULT

        # Suppress cascading widget callbacks during programmatic updates
        self._suppress = False

        # Matplotlib handles
        self.fig            = None
        self.ax_phase       = None
        self.ax_pot         = None
        self.ax_spec        = None
        self.click_cid      = None
        self._spec_eps_line = None   # red ε marker on spectrum plot

        self._build_widgets()
        self._build_figure()

    # =========================================================================
    # Widget construction
    # =========================================================================

    def _build_widgets(self):
        import ipywidgets as W
        from ipywidgets import Layout

        sl  = {'description_width': '110px'}
        sl2 = {'description_width': '90px'}
        lay_s  = Layout(width='280px')
        lay_s2 = Layout(width='240px')
        lay_drop = Layout(width='200px')
        lay_btn  = Layout(width='165px')
        lay_full = Layout(width='100%')
        lay_txt  = Layout(width='160px')

        # ── Fixed parameter sliders ───────────────────────────────────────────
        self.w_theta_NY = W.FloatSlider(
            value=self.SLIDER_THETA_NY_DEFAULT,
            min=self.SLIDER_THETA_NY_MIN, max=self.SLIDER_THETA_NY_MAX, step=0.1,
            description='theta_NY:', style=sl, layout=lay_s, readout_format='.2f')

        self.w_V_val = W.FloatSlider(
            value=self.SLIDER_V_DEFAULT,
            min=self.SLIDER_V_MIN, max=self.SLIDER_V_MAX, step=0.1,
            description='V (torsion):', style=sl, layout=lay_s, readout_format='.2f')

        self.w_alpha = W.FloatSlider(
            value=self.SLIDER_ALPHA_DEFAULT,
            min=self.SLIDER_ALPHA_MIN, max=self.SLIDER_ALPHA_MAX, step=0.05,
            description='alpha (Weyl):', style=sl, layout=lay_s, readout_format='.3f')

        # ── Phase diagram range ───────────────────────────────────────────────
        self.w_eps_min = W.FloatSlider(
            value=self.SLIDER_EPS_MIN_DEFAULT,
            min=self.SLIDER_EPS_MIN_MIN, max=self.SLIDER_EPS_MIN_MAX, step=0.05,
            description='eps min:', style=sl2, layout=lay_s2, readout_format='.2f')

        self.w_eps_max = W.FloatSlider(
            value=self.SLIDER_EPS_MAX_DEFAULT,
            min=self.SLIDER_EPS_MAX_MIN, max=self.SLIDER_EPS_MAX_MAX, step=0.05,
            description='eps max:', style=sl2, layout=lay_s2, readout_format='.2f')

        self.w_eta_min = W.FloatSlider(
            value=self.SLIDER_ETA_MIN_DEFAULT,
            min=self.SLIDER_ETA_MIN_MIN, max=self.SLIDER_ETA_MIN_MAX, step=1.0,
            description='eta min:', style=sl2, layout=lay_s2, readout_format='.1f')

        self.w_eta_max = W.FloatSlider(
            value=self.SLIDER_ETA_MAX_DEFAULT,
            min=self.SLIDER_ETA_MAX_MIN, max=self.SLIDER_ETA_MAX_MAX, step=1.0,
            description='eta max:', style=sl2, layout=lay_s2, readout_format='.1f')

        # ── Potential plot range ──────────────────────────────────────────────
        self.w_r_max = W.FloatSlider(
            value=self.SLIDER_R_MAX_DEFAULT,
            min=self.SLIDER_R_MAX_MIN, max=self.SLIDER_R_MAX_MAX, step=0.01,
            description='r max:', style=sl2, layout=lay_s2, readout_format='.2f')

        self.w_veff_ymax = W.FloatLogSlider(
            value=10000., base=10, min=0., max=8., step=0.1,
            description='Veff y-max:', style=sl2, layout=lay_s2, readout_format='.1e')
        self.w_veff_ymax_auto = W.Checkbox(
            value=True, description='auto',
            style={'description_width': '10px'},
            layout=Layout(width='80px'))

        self.w_r_scale = W.ToggleButtons(
            options=['Linear', 'Log'], value='Linear',
            description='r axis:', style=sl2)

        self.w_Veff_scale = W.ToggleButtons(
            options=['Linear', 'Symlog'], value='Linear',
            description='V_eff axis:', style=sl2)

        # ── Spectrum range ────────────────────────────────────────────────────
        self.w_spec_ymax = W.FloatLogSlider(
            value=1000., base=10, min=0., max=8., step=0.1,
            description='+y-max:', style=sl2, layout=lay_s2, readout_format='.1e')
        self.w_spec_ymin = W.FloatLogSlider(
            value=1000., base=10, min=0., max=8., step=0.1,
            description='-y-min:', style=sl2, layout=lay_s2, readout_format='.1e')
        self.w_spec_ymax_auto = W.Checkbox(
            value=True, description='auto',
            style={'description_width': '10px'},
            layout=Layout(width='80px'))

        self.w_spec_show_branches = W.Checkbox(
            value=True, description='Show 3 branches (m=0,±1,±2)',
            style={'description_width': '10px'},
            layout=Layout(width='260px'))

        # ── Point selection ───────────────────────────────────────────────────
        self.w_point_eps = W.FloatSlider(
            value=self.selected_point['eps'],
            min=self.SLIDER_EPS_MIN_DEFAULT,
            max=self.SLIDER_EPS_MAX_DEFAULT,
            step=0.01, description='epsilon:',
            style=sl2, layout=lay_s2, readout_format='.3f')

        self.w_point_eta = W.FloatSlider(
            value=self.selected_point['eta'],
            min=self.SLIDER_ETA_MIN_DEFAULT,
            max=self.SLIDER_ETA_MAX_DEFAULT,
            step=0.1, description='eta:',
            style=sl2, layout=lay_s2, readout_format='.2f')

        # ── Buttons ───────────────────────────────────────────────────────────
        self.w_draw_btn  = W.Button(
            description='Draw Phase Diagram',
            button_style='primary', layout=lay_btn)
        self.w_pot_btn   = W.Button(
            description='Update Potential',
            button_style='info', layout=lay_btn)
        self.w_spec_btn  = W.Button(
            description='Update Spectrum',
            button_style='success', layout=lay_btn)

        # ── Status ────────────────────────────────────────────────────────────
        self.w_status = W.HTML(
            value='<i style="color:gray;">Ready<br/>- click "Draw Phase Diagram".</i>',
            layout=lay_full)

        # ── Wire observers ────────────────────────────────────────────────────
        self.w_draw_btn .on_click(self._on_draw_clicked)
        self.w_pot_btn  .on_click(self._on_potential_clicked)
        self.w_spec_btn .on_click(self._on_spectrum_clicked)

        self.w_point_eps.observe(self._on_point_coords_changed, names='value')
        self.w_point_eta.observe(self._on_point_coords_changed, names='value')

        # ── Assemble control panel (top HBox layout) ──────────────────────────
        sep = W.HTML('<hr style="margin:4px 0; border-color:#ccc;">')

        hdrstyle = 'font-size:12px; font-weight:bold; margin:2px 0 4px 0;'

        grp_fixed = W.VBox([
            W.HTML(f'<div style="{hdrstyle}">Fixed Parameters</div>'),
            self.w_theta_NY, self.w_V_val, self.w_alpha,
        ], layout=Layout(padding='0 12px 0 0'))

        grp_range = W.VBox([
            W.HTML(f'<div style="{hdrstyle}">Phase Diagram Range</div>'),
            self.w_eps_min, self.w_eps_max,
            self.w_eta_min, self.w_eta_max,
            self.w_draw_btn,
        ], layout=Layout(padding='0 12px 0 0'))

        grp_point = W.VBox([
            W.HTML(f'<div style="{hdrstyle}">Point Selection</div>'),
            self.w_point_eps, self.w_point_eta,
            self.w_status,
        ], layout=Layout(padding='0 12px 0 0'))

        grp_pot = W.VBox([
            W.HTML(f'<div style="{hdrstyle}">Potential Range</div>'),
            self.w_r_max,
            W.HBox([self.w_veff_ymax, self.w_veff_ymax_auto]),
            self.w_r_scale, self.w_Veff_scale,
            self.w_pot_btn,
        ], layout=Layout(padding='0 12px 0 0'))

        grp_spec = W.VBox([
            W.HTML(f'<div style="{hdrstyle}">Spectrum Range</div>'),
            self.w_spec_show_branches,
            W.HBox([self.w_spec_ymax, self.w_spec_ymax_auto]),
            self.w_spec_ymin,
            self.w_spec_btn,
        ], layout=Layout(padding='0 12px 0 0'))

        self.w_controls = W.HBox([
            grp_fixed, grp_range, grp_point, grp_pot, grp_spec,
        ], layout=W.Layout(
            flex_flow='row wrap',
            padding='6px 8px',
            border='1px solid #ddd',
            background_color='#fafafa',
        ))

    # =========================================================================
    # Figure construction
    # =========================================================================

    def _build_figure(self):
        # 3 plots side by side: phase (wider) | potential | spectrum
        # Use ioff() to prevent automatic display in %matplotlib widget mode
        # (figure is rendered via viewer.display() instead)
        with plt.ioff():
            self.fig = plt.figure(figsize=(17, 5))
        gs = GridSpec(1, 3, figure=self.fig,
                      width_ratios=[2.2, 1.5, 1.5],
                      wspace=0.32,
                      left=0.06, right=0.97, top=0.90, bottom=0.13)

        self.ax_phase = self.fig.add_subplot(gs[0, 0])
        self.ax_pot   = self.fig.add_subplot(gs[0, 1])
        self.ax_spec  = self.fig.add_subplot(gs[0, 2])

        self.ax_phase.set_title(
            'Phase Diagram  (ε vs η)\n'
            'Click "Draw Phase Diagram" to compute', fontsize=10)
        self.ax_pot  .set_title('V_eff(r)', fontsize=10)
        self.ax_spec .set_title('Spin-2 spectrum  H_εε(ε)', fontsize=10)

    # =========================================================================
    # Public entry point
    # =========================================================================

    def display(self):
        """Render the interactive UI in a Jupyter notebook."""
        import ipywidgets as W
        from IPython.display import display as _display

        # Controls on top, figure canvas below
        ui = W.VBox([self.w_controls, self.fig.canvas])
        _display(ui)

        # (Re-)connect click handler
        if self.click_cid is not None:
            self.fig.canvas.mpl_disconnect(self.click_cid)
        self.click_cid = self.fig.canvas.mpl_connect(
            'button_press_event', self._on_click)

    # =========================================================================
    # Internal helpers
    # =========================================================================

    def _read_sliders(self):
        self.theta_NY = self.w_theta_NY.value
        self.V_val    = self.w_V_val.value
        self.alpha    = self.w_alpha.value
        self.eps_min  = self.w_eps_min.value
        self.eps_max  = self.w_eps_max.value
        self.eta_min  = self.w_eta_min.value
        self.eta_max  = self.w_eta_max.value
        self.r_max    = self.w_r_max.value

        # Keep point sliders bounded within diagram range
        for w, lo, hi in [
            (self.w_point_eps, self.eps_min, self.eps_max),
            (self.w_point_eta, self.eta_min, self.eta_max),
        ]:
            if w.min != lo or w.max != hi:
                self._suppress = True
                w.min, w.max = lo, hi
                w.value = float(np.clip(w.value, lo, hi))
                self._suppress = False

    def _vfunc_at(self, r, eps, eta):
        return self.vfunc(r, eps, eta, self.V_val, self.theta_NY, self.alpha)

    def _set_status(self, msg, color='black'):
        self.w_status.value = (
            f'<div style="color:{color}; font-size:11px; '
            f'white-space:normal; overflow-wrap:break-word; word-break:break-word;">'
            f'{msg}</div>'
        )

    # =========================================================================
    # Event handlers
    # =========================================================================

    def _on_draw_clicked(self, _b):
        self._read_sliders()
        self._set_status(
            'Computing phase diagram ...<br/>  '
            f'- theta_NY={self.theta_NY:.2f}<br/>- V={self.V_val:.2f}<br/>- alpha={self.alpha:.3f})',
            color='orange')
        self.fig.canvas.draw_idle()

        self._run_phase_scan()
        self._draw_phase_diagram()
        self._draw_potential_plot()
        self._set_status(
            'Computing spin-2 spectrum ...<br/>'
            f'- theta_NY={self.theta_NY:.2f}<br/>- V={self.V_val:.2f}<br/>- alpha={self.alpha:.3f}',
            color='orange')
        self.fig.canvas.draw_idle()
        self._draw_spectrum_plot()
        self._set_status(
            f'Done<br/>- theta_NY={self.theta_NY:.2f}<br/>- V={self.V_val:.2f}<br/>- alpha={self.alpha:.3f}',
            color='green')
        self.fig.canvas.draw_idle()

    def _on_potential_clicked(self, _b):
        self._read_sliders()
        self._draw_potential_plot()
        self.fig.canvas.draw_idle()

    def _on_spectrum_clicked(self, _b):
        self._read_sliders()
        self._set_status(
            f'Computing spin-2 spectrum ...<br/>- theta_NY={self.theta_NY:.2f}<br/>- V={self.V_val:.2f}<br/>- alpha={self.alpha:.3f}',
            color='orange')
        self.fig.canvas.draw_idle()
        self._draw_spectrum_plot()
        self._set_status(
            f'Spectrum updated<br/>- theta_NY={self.theta_NY:.2f}<br/>- V={self.V_val:.2f}<br/>- alpha={self.alpha:.3f}',
            color='green')
        self.fig.canvas.draw_idle()

    def _on_point_coords_changed(self, _change):
        if self._suppress:
            return
        eta_changed = (self.w_point_eta.value != self.selected_point['eta'])
        self.selected_point['eps'] = self.w_point_eps.value
        self.selected_point['eta'] = self.w_point_eta.value
        if self.phase_data is not None:
            self._draw_phase_diagram()
        self._draw_potential_plot()
        if eta_changed:
            self._draw_spectrum_plot()
        else:
            self._update_spec_eps_line()
        self.fig.canvas.draw_idle()

    def _on_click(self, event):
        if event.button != 1:
            return

        if event.inaxes == self.ax_phase:
            eps_c = event.xdata
            eta_c = event.ydata
            if eps_c is None or eta_c is None:
                return
            eps_c = float(np.clip(eps_c, self.eps_min, self.eps_max))
            eta_c = float(np.clip(eta_c, self.eta_min, self.eta_max))

            self.selected_point.update(eps=eps_c, eta=eta_c)
            self._suppress = True
            self.w_point_eps.value = eps_c
            self.w_point_eta.value = eta_c
            self._suppress = False

            if self.phase_data is not None:
                self._draw_phase_diagram()
            self._draw_potential_plot()
            self._draw_spectrum_plot()
            self.fig.canvas.draw_idle()

        elif event.inaxes == self.ax_spec:
            eps_c = event.xdata
            if eps_c is None:
                return
            eps_c = float(np.clip(eps_c, self.eps_min, self.eps_max))

            self.selected_point['eps'] = eps_c
            self._suppress = True
            self.w_point_eps.value = eps_c
            self._suppress = False

            if self.phase_data is not None:
                self._draw_phase_diagram()
            self._draw_potential_plot()
            self._update_spec_eps_line()
            self.fig.canvas.draw_idle()

    def _update_spec_eps_line(self):
        """Move the red ε marker on the spectrum plot without redrawing the full spectrum."""
        if self._spec_eps_line is not None:
            self._spec_eps_line.set_xdata([self.selected_point['eps'],
                                           self.selected_point['eps']])

    # =========================================================================
    # Phase map scan
    # =========================================================================

    def _run_phase_scan(self):
        n_e = self.GRID_N_EPS
        n_n = self.GRID_N_ETA

        eps_arr = np.linspace(self.eps_min, self.eps_max, n_e)
        eta_arr = np.linspace(self.eta_min, self.eta_max, n_n)
        self.eps_grid, self.eta_grid = np.meshgrid(eps_arr, eta_arr)

        phase = np.zeros((n_n, n_e), dtype=int)
        Heps  = np.full((n_n, n_e), np.nan)
        Hs    = np.full((n_n, n_e), np.nan)

        r_lo, r_hi = 0.05, 100.0
        for j, eta_v in enumerate(eta_arr):
            for i, eps_v in enumerate(eps_arr):
                f = lambda r, _e=eps_v, _n=eta_v: self._vfunc_at(r, _e, _n)
                code, r0 = _classify_veff(f, r_lo, r_hi, alpha=self.alpha, epsilon=eps_v)
                phase[j, i] = code
                # Compute spin-2 Hessians at any point that has a local minimum
                if code in (1, 2, 4, 5, 6) and np.isfinite(r0):
                    Heps[j, i] = _d2V_deps2(
                        self.vfunc, eps_v, eta_v,
                        self.V_val, self.theta_NY, self.alpha, r0)
                    bg_np = dict(r=r0, eps=eps_v, s=0., q3=0., q4=0., q5=0.,
                                 eta=eta_v, V_torsion=self.V_val,
                                 theta_NY=self.theta_NY, alpha=self.alpha, L=1., kappa=1.)
                    Hs[j, i] = _d2V_np_branch('s', 0., bg_np) * 0.5

        self.phase_data = phase
        self.Heps_data  = Heps
        self.Hs_data    = Hs

    # =========================================================================
    # Phase diagram drawing
    # =========================================================================

    def _draw_phase_diagram(self):
        ax = self.ax_phase
        ax.clear()

        if self.phase_data is None:
            ax.set_title('Phase Diagram  (click "Draw Phase Diagram")', fontsize=11)
            return

        # Discrete colour map (codes 0–7, matches paper02 v4)
        #  0:err  1:I  2:II  3:III  4:G  5:II-W  6:I-W  7:III-G
        cmap = mcolors.ListedColormap([
            self.COLOR_ERR, self.COLOR_I,   self.COLOR_II, self.COLOR_III,
            self.COLOR_G,   self.COLOR_IIW, self.COLOR_IW, self.COLOR_IIIG,
        ])
        bounds = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5]
        norm = mcolors.BoundaryNorm(bounds, cmap.N)

        eps_arr = self.eps_grid[0, :]
        eta_arr = self.eta_grid[:, 0]
        ax.pcolormesh(eps_arr, eta_arr, self.phase_data,
                      cmap=cmap, norm=norm, shading='auto')

        # ── Type III-G hatching overlay (α>0, no local minimum) ──────────────
        if np.any(self.phase_data == 7):
            _iiig_bin = (self.phase_data == 7).astype(float)
            _n_before = len(ax.collections)
            ax.contourf(eps_arr, eta_arr, _iiig_bin,
                        levels=[0.5, 1.5], colors=['none'], hatches=['///'])
            for _coll in ax.collections[_n_before:]:
                _coll.set_edgecolor('#808080')
                _coll.set_linewidth(0.3)

        # Background-phase boundary contours (same colours as paper02)
        try:
            ax.contour(eps_arr, eta_arr, self.phase_data,
                       levels=[1.5, 2.5],
                       colors=['darkgreen', 'darkorange'],
                       linewidths=1.5)
        except Exception:
            pass

        # Spin-2 instability boundary  H_eps = 0  (m=0, navy dashed)
        if self.Heps_data is not None:
            _H = np.where(np.isnan(self.Heps_data), np.nanmedian(self.Heps_data),
                          self.Heps_data)
            if not np.all(np.isnan(_H)):
                try:
                    cs = ax.contour(eps_arr, eta_arr, _H,
                                    levels=[0.0], colors=['navy'],
                                    linewidths=2.2, linestyles='--')
                    ax.clabel(cs, fmt=r'$H_{\varepsilon\varepsilon}=0$',
                              fontsize=8, inline=True)
                except Exception:
                    pass

        # Spin-2 instability boundary  H_s = 0  (m=±2, brown dashed)
        if self.Hs_data is not None:
            _Hs = np.where(np.isnan(self.Hs_data), np.nanmedian(self.Hs_data),
                           self.Hs_data)
            if not np.all(np.isnan(_Hs)):
                try:
                    cs2 = ax.contour(eps_arr, eta_arr, _Hs,
                                     levels=[0.0], colors=['saddlebrown'],
                                     linewidths=2.0, linestyles='--')
                    ax.clabel(cs2, fmt=r'$H_{ss}=0$',
                              fontsize=8, inline=True)
                except Exception:
                    pass

        self._draw_point_markers(ax)

        # Legend — adapt to current α value (matches paper02 v4)
        # Note: within the scan, ε=0 cells are always classified as I/II/III
        # (effective_alpha=0 there), but ε≠0 cells use the full α.
        alpha_now = self.alpha
        if abs(alpha_now) < _ALPHA_THRESHOLD:
            phase_patches = [
                mpatches.Patch(color=self.COLOR_I,   label='Type I  (stable well)'),
                mpatches.Patch(color=self.COLOR_II,  label='Type II  (rolling)'),
                mpatches.Patch(color=self.COLOR_III, label='Type III  (no well)'),
            ]
        elif alpha_now > 0:
            phase_patches = [
                mpatches.Patch(color=self.COLOR_G,   label='Type G  (metastable, α>0)'),
                mpatches.Patch(color=self.COLOR_IIIG, label='Type III-G  (no well, α>0)',
                               hatch='///', edgecolor='#808080'),
            ]
        else:
            phase_patches = [
                mpatches.Patch(color=self.COLOR_IW,  label='Type I-W  (double well, α<0)'),
                mpatches.Patch(color=self.COLOR_IIW, label='Type II-W  (single well, α<0)'),
            ]

        legend_handles = phase_patches + [
            Line2D([0], [0], color='navy', lw=2, ls='--',
                   label=r'$H_{\varepsilon\varepsilon}=0$ ($m{=}0$)'),
            Line2D([0], [0], color='saddlebrown', lw=2, ls='--',
                   label=r'$H_{ss}=0$ ($m{=}{\pm}2$)'),
            Line2D([0], [0], marker='o', color='w',
                   markerfacecolor='none', markeredgecolor='black',
                   markersize=9, markeredgewidth=2, label='Point'),
        ]

        ax.legend(handles=legend_handles, loc='upper right',
                  fontsize=8, framealpha=0.85)

        ax.set_xlabel(r'$\varepsilon$  (squash)', fontsize=11)
        ax.set_ylabel(r'$\eta$  (axial torsion)', fontsize=11)
        ax.set_title(
            r'Phase Diagram  ($\varepsilon$, $\eta$)'
            f'  |  theta_NY={self.theta_NY:.2f}'
            f'  V={self.V_val:.2f}'
            f'  alpha={self.alpha:.3f}',
            fontsize=11)
        ax.grid(True, alpha=0.18, linestyle=':')

    def _draw_point_markers(self, ax):
        pt = self.selected_point
        ax.scatter(pt['eps'], pt['eta'],
                   marker='o', s=120, facecolors='none',
                   edgecolors='black', linewidths=2.2, zorder=10)

    # =========================================================================
    # Potential plot
    # =========================================================================

    def _draw_potential_plot(self):
        ax = self.ax_pot
        ax.clear()

        log_r  = (self.w_r_scale.value == 'Log')
        symlog = (self.w_Veff_scale.value == 'Symlog')
        r_lo   = 0.01
        r_hi   = self.w_r_max.value

        r_arr = (np.logspace(np.log10(r_lo), np.log10(r_hi), 600)
                 if log_r else np.linspace(r_lo, r_hi, 600))

        has_plot = False
        pt = self.selected_point
        V_arr = np.array(
            [self._vfunc_at(r, pt['eps'], pt['eta']) for r in r_arr])
        ax.plot(r_arr, V_arr, color='black', linestyle='-', linewidth=2,
                label=(r'$\varepsilon$='
                       f"{pt['eps']:.3f}, "
                       r'$\eta$='
                       f"{pt['eta']:.2f}"))
        has_plot = True

        if not has_plot:
            ax.text(0.5, 0.5, 'No active points selected',
                    ha='center', va='center', transform=ax.transAxes,
                    color='gray', fontsize=10)
            ax.set_title('Effective Potential  V_eff(r)', fontsize=10)
            return

        ax.axhline(0, color='gray', linestyle='--', linewidth=0.8)
        if log_r:
            ax.set_xscale('log')
        if symlog:
            ax.set_yscale('symlog', linthresh=1.0)

        # Apply y-axis upper limit (only when auto is unchecked)
        if not self.w_veff_ymax_auto.value:
            ax.set_ylim(top=self.w_veff_ymax.value)

        ax.set_xlabel(r'$r$', fontsize=10)
        ax.set_ylabel(r'$V_{\mathrm{eff}}(r)$', fontsize=10)
        ax.set_title(
            f'V_eff(r)  |  theta_NY={self.theta_NY:.2f}'
            f'  V={self.V_val:.2f}  alpha={self.alpha:.3f}',
            fontsize=10)
        ax.legend(fontsize=8, loc='best')
        ax.grid(True, alpha=0.3, linestyle='--')

    # =========================================================================
    # Spin-2 spectrum plot
    # =========================================================================

    def _draw_spectrum_plot(self):
        """
        Plot spin-2 Hessian branches vs epsilon at each selected point's eta.

        3 branches (when 'Show 3 branches' is checked):
          m=0   (T1, squash):  d2V/deps2   — same as original
          m=±2  (T2, shear):   d2V/ds2     — pure-numpy computation
          m=±1  (T3, mixing):  d2V/dq3^2  — pure-numpy computation
        """
        ax = self.ax_spec
        ax.clear()

        n_eps        = 80
        eps_arr      = np.linspace(self.eps_min, self.eps_max, n_eps)
        r_lo, r_hi   = 0.05, 100.0
        show_branches = self.w_spec_show_branches.value

        # Branch style: linestyle distinguishes branch; color distinguishes point
        branch_ls    = ['-', '--', ':']
        branch_names = [r'$m{=}0$ (squash)', r'$m{=}{\pm}2$ (shear)',
                        r'$m{=}{\pm}1$ (mixing)']

        has_plot = False
        pt    = self.selected_point
        eta_v = pt['eta']
        lbl_prefix = r'$\eta$=' + f'{eta_v:.2f}'

        # ── Branch 0: m=0 (squash, d2V/deps2) — via lambdified vfunc ────────
        H0_arr = []
        r0_arr = []
        for eps_v in eps_arr:
            f = lambda r, _e=eps_v, _n=eta_v: self._vfunc_at(r, _e, _n)
            r0, _ = _find_r0(f, r_lo, r_hi)
            if np.isfinite(r0):
                h0 = _d2V_deps2(self.vfunc, eps_v, eta_v,
                                self.V_val, self.theta_NY, self.alpha, r0)
            else:
                h0 = float('nan')
            H0_arr.append(h0)
            r0_arr.append(r0)
        # Canonical normalization: q1 = sqrt(2/3)*log(1+eps)
        # => d2V/dq1^2 = (3/2)(1+eps)^2 * d2V/deps^2  (at dV/deps=0)
        H0_arr = np.array(H0_arr) * (1.5 * (1. + eps_arr)**2)

        ax.plot(eps_arr, H0_arr, color='black', linestyle=self.BRANCH_LINESTYLES[0],
                linewidth=2,
                label=(lbl_prefix + '  ' + branch_names[0]
                       if show_branches else lbl_prefix))
        has_plot = True

        if show_branches:
            # ── Branches 1 & 2: m=±2 and m=±1 — pure-numpy Veff ─────────────
            # _compute_Veff_np now includes the Weyl (α·C²) term, so all three
            # branches are consistent. r0_arr (engine equilibrium) is shared.
            for br_idx, br_key, br_val0 in [(1, 's', 0.), (2, 'q4', 0.)]:
                H_arr = []
                for i_e, eps_v in enumerate(eps_arr):
                    # Shared r0 with m=0 branch (engine equilibrium, includes α)
                    r0_np = r0_arr[i_e]
                    if not np.isfinite(r0_np):
                        H_arr.append(float('nan'))
                        continue
                    bg_np = dict(r=r0_np, eps=eps_v, s=0., q3=0., q4=0., q5=0.,
                                 eta=eta_v, V_torsion=self.V_val,
                                 theta_NY=self.theta_NY, alpha=self.alpha, L=1., kappa=1.)
                    hval = _d2V_np_branch(br_key, br_val0, bg_np)
                    H_arr.append(hval)
                # Canonical normalization:
                #   m=±2 (s): s = exp(q2/√2)−1 → ds/dq2=(1+s)/√2 → factor 1/2
                #   m=±1 (q4): G4 = exp(q4/√2 * M4), M4=[[0,0,1],[0,0,0],[1,0,0]] → factor = 1
                #   Note: q3 (G3, 0-1 mixing) is T3 which is m=±2 like s — using q4 (0-2 mixing)
                #         gives the correct m=±1 mode (T4) under residual SO(2) in the 0-1 plane.
                norm = 0.5 if br_key == 's' else 1.0
                H_arr = np.array(H_arr) * norm
                ax.plot(eps_arr, H_arr, color='black',
                        linestyle=self.BRANCH_LINESTYLES[br_idx],
                        linewidth=1.8,
                        label=lbl_prefix + '  ' + branch_names[br_idx])

        if not has_plot:
            ax.text(0.5, 0.5, 'No active points selected',
                    ha='center', va='center', transform=ax.transAxes,
                    color='gray', fontsize=10)
        else:
            ax.axhline(0, color='gray', linestyle='--', linewidth=1.0)
            ax.axvline(0, color='lightgray', linestyle=':', linewidth=0.8)
            self._spec_eps_line = ax.axvline(
                pt['eps'], color='red', linestyle='-', linewidth=0.8)
            ax.set_xlabel(r'$\varepsilon$  (squash)', fontsize=10)
            ax.set_ylabel(r'$\partial^2 V_{\rm eff}$', fontsize=9)

            # Apply y-axis range (only when auto is unchecked)
            if not self.w_spec_ymax_auto.value:
                ax.set_ylim(-self.w_spec_ymin.value, self.w_spec_ymax.value)

            ax.legend(fontsize=7, loc='upper left')
            ax.grid(True, alpha=0.3, linestyle='--')

        title_suffix = ' (3 branches)' if show_branches else ''
        ax.set_title(
            r'Spin-2 spectrum$^{(\varepsilon)}$' + title_suffix + '\n'
            r'solid=$m{=}0$, dashed=$m{=}{\pm}2$, dot=$m{=}{\pm}1$',
            fontsize=9)
