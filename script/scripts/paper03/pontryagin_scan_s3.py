"""
S³ isotropic twist — Pontryagin scan (App B.7)
=================================================

Numerical implementation of isotropic twist analysis (highest priority, analysis frontier).

Strategy: numerics first → symbolic verification
    1. Run TwistS3S1Engine symbolically (keeping ω as a SymPy symbol)
    2. Lambdify R_{abcd}(r, ω, η, ...) for numerical evaluation
    3. Scan P(ω) = <R, *R> over ω ∈ [-1, 1]
    4. Check ω → -ω symmetry (expected: odd function P(-ω) = -P(ω))
    5. Scan V_eff(r, ω) in r → determine sign of m²(r*) = ∂²V/∂ω²|_{ω=0}

Done conditions:
    □ ω → -ω symmetry test passed
    □ Polynomial order structure (odd/even) of P(ω) confirmed
    □ Stability of ω=0 in V_eff(r, ω) determined (sign of m²(r*))

Author: Muacca
Date: 2026-03-15
"""

import sys
import os
import io

# Windows cp932 workaround
if sys.stdout.encoding and sys.stdout.encoding.lower() not in ('utf-8', 'utf8'):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

import argparse
import time
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from dppu.topology.unified import UnifiedEngine, DOFConfig

from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant
from dppu.curvature.hodge import compute_hodge_dual
from dppu.engine.logger import ComputationLogger, NullLogger


# ===========================================================================
# Default settings
# ===========================================================================

DEFAULT_PARAMS = {
    'L': 1.0,
    'kappa': 1.0,
    'theta_NY': 1.0,
    'alpha': -1.0,     # α≤0 fixed
    'V': 0.0,          # AX mode: V=0
}

OMEGA_SCAN = np.linspace(-1.0, 1.0, 41)   # ω ∈ [-1, 1], 41 points
R_SCAN = np.array([0.5, 1.0, 2.0, 3.0, 5.0])    # representative r values
ETA_DEFAULT = -2.0     # η reference value


# ===========================================================================
# Utilities
# ===========================================================================

def compute_pontryagin(R_abcd: np.ndarray) -> float:
    """P = <R, *R> = Σ_{abcd} R_{abcd} (*R)^{abcd}"""
    R_dual = compute_hodge_dual(R_abcd)
    return float(np.dot(R_abcd.ravel(), R_dual.ravel()))


def check_odd_symmetry(omega_vals: np.ndarray, P_vals: np.ndarray,
                       tol: float = 1e-6) -> dict:
    """
    Numerically verify the odd function property P(-ω) = -P(ω).

    Returns
    -------
    dict: {'is_odd': bool, 'max_asymmetry': float, 'details': list}
    """
    # Find pairs of ω and -ω
    pairs = []
    for i, w in enumerate(omega_vals):
        if w > 0:
            j = np.argmin(np.abs(omega_vals + w))
            if np.abs(omega_vals[j] + w) < 1e-8:
                pairs.append((i, j, w, P_vals[i], P_vals[j]))

    if not pairs:
        return {'is_odd': None, 'max_asymmetry': None, 'details': []}

    asymmetries = []
    for i_pos, i_neg, w, P_pos, P_neg in pairs:
        # Odd function: P(+ω) + P(-ω) ≈ 0
        asym = abs(P_pos + P_neg)
        # Relative error (when P is non-zero)
        scale = max(abs(P_pos), abs(P_neg), 1e-30)
        rel_asym = asym / scale
        asymmetries.append((w, P_pos, P_neg, asym, rel_asym))

    max_rel_asym = max(a[4] for a in asymmetries)
    is_odd = max_rel_asym < tol

    return {
        'is_odd': is_odd,
        'max_asymmetry': max_rel_asym,
        'details': asymmetries[:5],  # first 5 pairs only
    }


def fit_polynomial_order(omega_vals: np.ndarray, P_vals: np.ndarray) -> dict:
    """
    Numerically estimate the polynomial order structure of P(ω).

    If odd function: P(ω) ≈ c₁ω + c₃ω³ + ...
    Sign of c₁ matters (positive may mean V_eff ω-dependence does not stabilize)

    Returns
    -------
    dict: {'c1': float, 'c3': float, 'r_squared': float}
    """
    # Use only the ω > 0 part (coefficient fit for odd function)
    mask = omega_vals > 0.01
    w = omega_vals[mask]
    P = P_vals[mask]

    if len(w) < 3:
        return {'c1': None, 'c3': None, 'r_squared': None}

    # Fit P = c1*ω + c3*ω³
    A = np.column_stack([w, w**3])
    try:
        coeffs, residuals, rank, sv = np.linalg.lstsq(A, P, rcond=None)
        c1, c3 = coeffs

        # Fit quality (R²)
        P_pred = A @ coeffs
        ss_res = np.sum((P - P_pred)**2)
        ss_tot = np.sum((P - np.mean(P))**2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0

        return {'c1': float(c1), 'c3': float(c3), 'r_squared': float(r2)}
    except Exception:
        return {'c1': None, 'c3': None, 'r_squared': None}


def compute_mass_squared(V_func, r_val: float, eta_val: float,
                         params: dict, delta_omega: float = 0.01) -> float:
    """
    Compute m²(r) = ∂²V_eff/∂ω²|_{ω=0} by numerical differentiation.

    m² > 0 → ω=0 stable (chiral equilibrium maintained)
    m² < 0 → ω=0 unstable (spontaneous twist generation)

    V_func: callable f(r, omega, eta, theta_NY, L, kappa, alpha) -> float
    """
    w = delta_omega
    V_plus = V_func(r_val, w, eta_val,
                    params['theta_NY'], params['L'],
                    params['kappa'], params['alpha'])
    V_zero = V_func(r_val, 0.0, eta_val,
                    params['theta_NY'], params['L'],
                    params['kappa'], params['alpha'])
    V_minus = V_func(r_val, -w, eta_val,
                     params['theta_NY'], params['L'],
                     params['kappa'], params['alpha'])

    m2 = (V_plus - 2 * V_zero + V_minus) / w**2
    return float(m2)


# ===========================================================================
# Main analysis
# ===========================================================================

def run_phase3_0a(r_val: float, eta_val: float, verbose: bool = True):
    """
    Pontryagin scan for S³ isotropic twist (App B.7).

    Parameters
    ----------
    r_val : float
        Fixed r value
    eta_val : float
        Torsion coefficient η
    """
    print(f"\n{'='*65}")
    print(f"S³ isotropic twist scan  (r={r_val}, η={eta_val})")
    print(f"{'='*65}")

    # ------------------------------------------------------------------
    # Step 1: Run engine (symbolic computation)
    # ------------------------------------------------------------------
    print("\n[Step 1] Running TwistS3S1Engine symbolically...")
    t0 = time.time()

    mode = Mode.AX  # AX: axial torsion (V=0, eta≠0)
    _cfg_engine = DOFConfig.from_engine('TwistS3S1Engine')
    _cfg_engine.torsion_mode = mode
    _cfg_engine.ny_variant = NyVariant.FULL
    engine = UnifiedEngine(_cfg_engine, NullLogger())
    engine.run()

    t1 = time.time()
    print(f"  Done ({t1-t0:.1f}s)")
    print(f"  V_eff = {engine.data['potential']}")

    # ------------------------------------------------------------------
    # Step 2: Lambdify R_{abcd} and V_eff
    # ------------------------------------------------------------------
    print("\n[Step 2] Lambdifying R_{abcd} and V_eff...")

    R_func, arg_names = engine.get_riemann_lambdified()
    print(f"  R argument order: {arg_names}")

    # Lambdify V_eff
    from sympy import lambdify, S as Szero
    params_sym = engine.data['params']
    V_sym = engine.data['potential']

    # AX mode: V=0, so arguments are (r, omega, eta, theta_NY, L, kappa, alpha)
    v_arg_syms = [
        params_sym['r'], params_sym['omega'], params_sym['eta'],
        params_sym['theta_NY'], params_sym['L'],
        params_sym['kappa'], params_sym['alpha']
    ]
    V_func = lambdify(v_arg_syms, V_sym, modules='numpy')

    t2 = time.time()
    print(f"  Done ({t2-t1:.1f}s)")

    # ------------------------------------------------------------------
    # Step 3: Numerical scan of P(ω)
    # ------------------------------------------------------------------
    print(f"\n[Step 3] P(ω) scan (r={r_val}, η={eta_val})")
    print(f"  ω ∈ [{OMEGA_SCAN[0]:.1f}, {OMEGA_SCAN[-1]:.1f}], {len(OMEGA_SCAN)} points")

    P_vals = []
    for omega in OMEGA_SCAN:
        # argument order follows arg_names
        arg_map = {
            'r': r_val, 'omega': omega, 'eta': eta_val,
            'L': DEFAULT_PARAMS['L'], 'kappa': DEFAULT_PARAMS['kappa'],
            'theta_NY': DEFAULT_PARAMS['theta_NY'],
            'alpha': DEFAULT_PARAMS['alpha']
        }
        args = [arg_map[name] for name in arg_names]
        R = R_func(*args)
        P = compute_pontryagin(R)
        P_vals.append(P)

    P_vals = np.array(P_vals)
    max_P_scan = np.max(np.abs(P_vals))

    # Display results
    print(f"\n  {'ω':>8}  {'P(ω)':>16}")
    print(f"  {'-'*8}  {'-'*16}")
    for i in range(0, len(OMEGA_SCAN), 5):
        print(f"  {OMEGA_SCAN[i]:>8.3f}  {P_vals[i]:>16.6e}")

    # ------------------------------------------------------------------
    # Step 4: ω → -ω symmetry check (skip if P≡0)
    # ------------------------------------------------------------------
    print(f"\n[Step 4] ω → -ω symmetry check")
    if max_P_scan < 1e-10:
        print(f"  max|P| = {max_P_scan:.2e} → P is numerically zero; odd/even test skipped")
        sym_result = {'is_odd': None, 'max_asymmetry': None, 'details': []}
    else:
        print("  Checking odd function property (P(-ω) = -P(ω) ?)...")
        sym_result = check_odd_symmetry(OMEGA_SCAN, P_vals)
        if sym_result['is_odd'] is not None:
            print(f"  P(-ω) = -P(ω)? {'YES (odd function)' if sym_result['is_odd'] else 'NO'}")
            print(f"  Max relative asymmetry: {sym_result['max_asymmetry']:.2e}")
            if sym_result['details']:
                print("  Sample pairs:")
                for w, Ppos, Pneg, asym, rel_asym in sym_result['details'][:3]:
                    print(f"    ω={w:.2f}: P(+ω)={Ppos:.4e}, P(-ω)={Pneg:.4e}, "
                          f"|P++P-|/scale={rel_asym:.2e}")
        else:
            print("  [WARNING] Insufficient ω pairs found")

    # ------------------------------------------------------------------
    # Step 5: Polynomial order fit for P(ω) (skip if P≡0)
    # ------------------------------------------------------------------
    print(f"\n[Step 5] Polynomial order fit of P(ω) (P ≈ c₁ω + c₃ω³)")
    if max_P_scan < 1e-10:
        print(f"  max|P| = {max_P_scan:.2e} → P is numerically zero; polynomial fit skipped")
        fit = {'c1': None, 'c3': None, 'r_squared': None}
    else:
        fit = fit_polynomial_order(OMEGA_SCAN, P_vals)
        if fit['c1'] is not None:
            print(f"  c₁ = {fit['c1']:.6e}  (1st order coefficient)")
            print(f"  c₃ = {fit['c3']:.6e}  (3rd order coefficient)")
            print(f"  R² = {fit['r_squared']:.6f}  (fit quality)")
            if fit['r_squared'] > 0.99:
                print("  → P(ω) is well described by an odd polynomial")
            if abs(fit['c1']) > 1e-8:
                print(f"  → 1st order coefficient c₁ ≠ 0: P ≠ 0 for ω ≠ 0 confirmed")
            else:
                print(f"  → 1st order coefficient c₁ ≈ 0: P is higher order or zero")

    # ------------------------------------------------------------------
    # Step 6: V_eff(r, ω) — compute m²(r)
    # ------------------------------------------------------------------
    print(f"\n[Step 6] Scan of m²(r) = ∂²V_eff/∂ω²|_{{ω=0}}")
    print(f"  {'r':>6}  {'m²(r)':>14}  {'stability':>12}")
    print(f"  {'-'*6}  {'-'*14}  {'-'*12}")

    m2_results = []
    for r in R_SCAN:
        try:
            m2 = compute_mass_squared(V_func, r, eta_val, DEFAULT_PARAMS)
            stability = "stable (ω=0)" if m2 > 0 else "unstable → twist generation"
            print(f"  {r:>6.2f}  {m2:>14.6e}  {stability:>12}")
            m2_results.append((r, m2))
        except Exception as e:
            print(f"  {r:>6.2f}  [ERROR: {e}]")

    # Summary of m² signs
    m2_vals = [m2 for _, m2 in m2_results]
    if m2_vals:
        has_negative = any(m2 < 0 for m2 in m2_vals)
        min_m2_r, min_m2 = min(m2_results, key=lambda x: x[1])
        print(f"\n  min m²(r) = {min_m2:.4e} at r={min_m2_r:.2f}")
        if has_negative:
            print("  → m²(r*) < 0 exists: spontaneous twist generation possible!")
        else:
            print("  → m²(r) > 0 for all r: ω=0 is stable (chiral equilibrium maintained)")

    # ------------------------------------------------------------------
    # Done condition summary
    # ------------------------------------------------------------------
    print(f"\n{'='*65}")
    print("Done condition check (plan v5.0 §4.5)")
    print(f"{'='*65}")

    p_nonzero = np.max(np.abs(P_vals)) > 1e-8
    p_odd = sym_result.get('is_odd', False) or False
    m2_computed = len(m2_results) > 0

    print(f"  [{'OK' if p_nonzero else 'N/A (P≡0)'}] P(ω) ≠ 0 for ω ≠ 0 numerically confirmed")
    print(f"  [{'OK' if p_odd else ('N/A (P≡0)' if not p_nonzero else '--')}] ω → -ω symmetry (odd function property) verified")
    print(f"  [{'OK' if m2_computed else '--'}] sign of m²(r*) computed")
    print()
    print(f"  max|P(ω)|: {np.max(np.abs(P_vals)):.4e}")
    if fit['c1'] is not None:
        print(f"  c₁ (1st order coefficient): {fit['c1']:.4e}")

    return {
        'r_val': r_val, 'eta_val': eta_val,
        'omega_vals': OMEGA_SCAN, 'P_vals': P_vals,
        'symmetry': sym_result, 'fit': fit,
        'm2_results': m2_results,
    }


# ===========================================================================
# Entry point
# ===========================================================================

def main():
    parser = argparse.ArgumentParser(description='S³ isotropic twist Pontryagin scan (App B.7)')
    parser.add_argument('--r-val', type=float, default=1.0,
                        help='Fixed r value (default: 1.0)')
    parser.add_argument('--eta', type=float, default=ETA_DEFAULT,
                        help=f'Torsion coefficient η (default: {ETA_DEFAULT})')
    parser.add_argument('--all-r', action='store_true',
                        help='Scan all r values in R_SCAN')
    args = parser.parse_args()

    print("=" * 65)
    print("S³ isotropic twist — Pontryagin & V_eff scan (App B.7)")
    print("=" * 65)

    if args.all_r:
        results = []
        for r in R_SCAN:
            res = run_phase3_0a(r_val=r, eta_val=args.eta)
            results.append(res)
    else:
        results = [run_phase3_0a(r_val=args.r_val, eta_val=args.eta)]

    print("\n" + "=" * 65)
    print("Overall summary")
    print("=" * 65)
    print()
    print(f"{'r':>6}  {'η':>6}  {'max|P|':>12}  {'c₁':>12}  {'odd?':>6}")
    print(f"  {'-'*6}  {'-'*6}  {'-'*12}  {'-'*12}  {'-'*6}")
    for res in results:
        max_P = np.max(np.abs(res['P_vals']))
        c1 = res['fit'].get('c1')
        c1_str = f"{c1:.3e}" if c1 is not None else "N/A"
        odd_str = "YES" if res['symmetry'].get('is_odd') else "NO"
        print(f"  {res['r_val']:>6.2f}  {res['eta_val']:>6.2f}  "
              f"{max_P:>12.4e}  {c1_str:>12}  {odd_str:>6}")


if __name__ == '__main__':
    _t_start = time.time()
    main()
    print(f"\nComputation time: {time.time() - _t_start:.1f}s")
    print()
    print("*Script: scripts/paper03/pontryagin_scan_s3.py")
