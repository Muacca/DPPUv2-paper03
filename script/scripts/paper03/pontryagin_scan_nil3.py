"""
Nil³ isotropic twist — Pontryagin scan (App B.7)
===================================================

In the S³ isotropic twist result, P ≡ 0 was established symbolically.
Does the same cancellation occur in Nil³?

Features of Nil³:
    - Non-bi-invariant (no SO(3) symmetry)
    - Structure constants: only C^2_{01} is non-zero (off-diagonal)
    - Twist: C^3_{01} = -Lω/R² (single pair only, negative sign)

Expectations:
    - Cancellation in S³ is due to SO(3) symmetry or algebraic structure
    - With lower symmetry in Nil³, P ≠ 0 may occur
    - However, if P = 0 this provides further evidence for the "protection theorem"


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
from dppu.engine.logger import NullLogger
from sympy import lambdify, simplify, S as Szero


# ===========================================================================
# Default settings
# ===========================================================================

DEFAULT_PARAMS = {
    'L': 1.0,
    'kappa': 1.0,
    'theta_NY': 1.0,
    'alpha': -1.0,     # α≤0 fixed
}

OMEGA_SCAN = np.linspace(-1.0, 1.0, 41)
R_SCAN = np.array([0.5, 1.0, 2.0, 3.0, 5.0])
ETA_DEFAULT = -2.0


# ===========================================================================
# Utilities
# ===========================================================================

def compute_pontryagin(R_abcd: np.ndarray) -> float:
    """P = <R, *R>"""
    R_dual = compute_hodge_dual(R_abcd)
    return float(np.dot(R_abcd.ravel(), R_dual.ravel()))


def check_p_symmetry(omega_vals: np.ndarray, P_vals: np.ndarray,
                     tol: float = 1e-6) -> dict:
    """
    Check ω → -ω symmetry of P(ω).
    If P is an even function: P(+ω) - P(-ω) ≈ 0
    If P is an odd function: P(+ω) + P(-ω) ≈ 0
    """
    pairs = []
    for i, w in enumerate(omega_vals):
        if w > 1e-8:
            j = np.argmin(np.abs(omega_vals + w))
            if np.abs(omega_vals[j] + w) < 1e-8:
                pairs.append((w, P_vals[i], P_vals[j]))

    if not pairs:
        return {'is_even': None, 'is_odd': None}

    even_asym = max(abs(Ppos - Pneg) / max(abs(Ppos), abs(Pneg), 1e-30)
                    for _, Ppos, Pneg in pairs)
    odd_asym = max(abs(Ppos + Pneg) / max(abs(Ppos), abs(Pneg), 1e-30)
                   for _, Ppos, Pneg in pairs)

    return {
        'is_even': even_asym < tol,
        'is_odd': odd_asym < tol,
        'even_asym': even_asym,
        'odd_asym': odd_asym,
        'sample_pairs': pairs[:3],
    }


def compute_mass_squared(V_func, R_val: float, eta_val: float,
                         params: dict, delta_omega: float = 0.01) -> float:
    """m²(R) = ∂²V_eff/∂ω²|_{ω=0}"""
    w = delta_omega
    V_plus = V_func(R_val, w, eta_val,
                    params['theta_NY'], params['L'],
                    params['kappa'], params['alpha'])
    V_zero = V_func(R_val, 0.0, eta_val,
                    params['theta_NY'], params['L'],
                    params['kappa'], params['alpha'])
    V_minus = V_func(R_val, -w, eta_val,
                     params['theta_NY'], params['L'],
                     params['kappa'], params['alpha'])
    return float((V_plus - 2 * V_zero + V_minus) / w**2)


# ===========================================================================
# Symbolic verification
# ===========================================================================

def symbolic_verify(engine) -> dict:
    """
    Symbolically compute P = ⟨R, ∗R⟩ and check whether it is zero.

    Returns
    -------
    dict: {'P_symbolic': expr, 'P_is_zero': bool}
    """
    from sympy import S
    R_sym = engine.data['riemann_abcd']
    dim = engine.data['dim']

    print("\n  [Symbolic computation] Computing P = ⟨R, ∗R⟩...")

    # Hodge dual (second-pair convention: ε on last two indices)
    # P = Σ_{abcd} R_{abcd} * (ε_{cdef} R^{ab}_{ef})
    # Simplified: P = Σ R_{abcd} (*R)^{abcd}
    # Use numpy-based symbolic sum
    from dppu.utils.levi_civita import levi_civita_4d

    P_sym = S.Zero
    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                for d in range(dim):
                    R_abcd = R_sym[a, b, c, d]
                    if R_abcd == S.Zero:
                        continue
                    # (*R)^{abcd} = (1/2) ε^{abef} R_{ef}^{cd}
                    # For Euclidean metric: indices raised/lowered trivially
                    # (*R)_{abcd} = (1/2) ε_{abef} R^{ef}_{cd}
                    for e in range(dim):
                        for f in range(dim):
                            eps_val = levi_civita_4d(a, b, e, f)
                            if eps_val == 0:
                                continue
                            R_efcd = R_sym[e, f, c, d]
                            if R_efcd == S.Zero:
                                continue
                            P_sym += R_abcd * eps_val * R_efcd

    P_sym = P_sym / 2  # factor from the Hodge dual definition

    print(f"  P (before simplify): {'zero' if P_sym == S.Zero else 'non-zero expression'}")

    P_simplified = simplify(P_sym)
    P_is_zero = (P_simplified == S.Zero)

    print(f"  P (after simplify):  {'0 [CONFIRMED]' if P_is_zero else str(P_simplified)[:120]}")

    return {'P_symbolic': P_simplified, 'P_is_zero': P_is_zero}


# ===========================================================================
# Main analysis
# ===========================================================================

def run_phase3_0b(R_val: float, eta_val: float, do_symbolic: bool = True):
    """
    Pontryagin scan for Nil³ isotropic twist (App B.7).
    """
    print(f"\n{'='*65}")
    print(f"Nil³ isotropic twist scan  (R={R_val}, η={eta_val})")
    print(f"{'='*65}")

    # ------------------------------------------------------------------
    # Step 1: Run engine (symbolic computation)
    # ------------------------------------------------------------------
    print("\n[Step 1] Running TwistNil3S1Engine symbolically...")
    t0 = time.time()

    _cfg_engine = DOFConfig.from_engine('TwistNil3S1Engine')
    _cfg_engine.torsion_mode = Mode.AX
    _cfg_engine.ny_variant = NyVariant.FULL
    engine = UnifiedEngine(_cfg_engine, NullLogger())
    engine.run()

    t1 = time.time()
    print(f"  Done ({t1-t0:.1f}s)")
    print(f"  V_eff = {engine.data['potential']}")

    # ------------------------------------------------------------------
    # Step 2: Lambdify
    # ------------------------------------------------------------------
    print("\n[Step 2] Lambdifying R_{abcd} and V_eff...")

    R_func, arg_names = engine.get_riemann_lambdified()
    print(f"  R argument order: {arg_names}")

    # V_eff lambdify (AX mode: V=0, args are R, omega, eta, theta_NY, L, kappa, alpha)
    params_sym = engine.data['params']
    V_sym = engine.data['potential']
    v_arg_syms = [
        params_sym['R'], params_sym['omega'], params_sym['eta'],
        params_sym['theta_NY'], params_sym['L'],
        params_sym['kappa'], params_sym['alpha']
    ]
    V_func = lambdify(v_arg_syms, V_sym, modules='numpy')

    t2 = time.time()
    print(f"  Done ({t2-t1:.1f}s)")

    # ------------------------------------------------------------------
    # Step 3: Numerical scan of P(ω)
    # ------------------------------------------------------------------
    print(f"\n[Step 3] Numerical scan of P(ω) (R={R_val}, η={eta_val})")

    P_vals = []
    for omega in OMEGA_SCAN:
        arg_map = {
            'R': R_val, 'omega': omega, 'eta': eta_val,
            'L': DEFAULT_PARAMS['L'], 'kappa': DEFAULT_PARAMS['kappa'],
            'theta_NY': DEFAULT_PARAMS['theta_NY'],
            'alpha': DEFAULT_PARAMS['alpha']
        }
        args = [arg_map[name] for name in arg_names]
        R_num = R_func(*args)
        P = compute_pontryagin(R_num)
        P_vals.append(P)

    P_vals = np.array(P_vals)
    max_P = np.max(np.abs(P_vals))

    print(f"\n  {'ω':>8}  {'P(ω)':>18}")
    print(f"  {'-'*8}  {'-'*18}")
    for i in range(0, len(OMEGA_SCAN), 5):
        print(f"  {OMEGA_SCAN[i]:>8.3f}  {P_vals[i]:>18.6e}")

    print(f"\n  max|P(ω)| = {max_P:.4e}")
    if max_P < 1e-10:
        print("  → P ≡ 0 (machine precision)")
    else:
        print(f"  → P ≠ 0 confirmed! (P≠0 at ω≠0)")

    # ------------------------------------------------------------------
    # Step 4: ω → -ω symmetry check (skip if P≡0)
    # ------------------------------------------------------------------
    print(f"\n[Step 4] Symmetry check of P(ω)")
    if max_P < 1e-10:
        print(f"  P is numerically zero; symmetry test skipped")
        sym_result = {'is_even': None, 'is_odd': None}
    else:
        sym_result = check_p_symmetry(OMEGA_SCAN, P_vals)
        print(f"  P(-ω) = P(ω)?  (even function): "
              f"{'YES' if sym_result['is_even'] else 'NO'}  "
              f"(max_asym={sym_result.get('even_asym', 'N/A'):.2e})")
        print(f"  P(-ω) = -P(ω)? (odd function): "
              f"{'YES' if sym_result['is_odd'] else 'NO'}  "
              f"(max_asym={sym_result.get('odd_asym', 'N/A'):.2e})")
        if sym_result['sample_pairs']:
            print("  Sample pairs:")
            for w, Ppos, Pneg in sym_result['sample_pairs'][:3]:
                print(f"    ω={w:.2f}: P(+ω)={Ppos:.4e}, P(-ω)={Pneg:.4e}")

    # ------------------------------------------------------------------
    # Step 5: V_eff → m²(R) scan
    # ------------------------------------------------------------------
    print(f"\n[Step 5] Scan of m²(R) = ∂²V_eff/∂ω²|_{{ω=0}}")
    print(f"  {'R':>6}  {'m²(R)':>14}  {'stability':>14}")
    print(f"  {'-'*6}  {'-'*14}  {'-'*14}")

    m2_results = []
    for Rv in R_SCAN:
        try:
            m2 = compute_mass_squared(V_func, Rv, eta_val, DEFAULT_PARAMS)
            stability = "stable (ω=0)" if m2 > 0 else "unstable → twist generation"
            print(f"  {Rv:>6.2f}  {m2:>14.6e}  {stability:>14}")
            m2_results.append((Rv, m2))
        except Exception as e:
            print(f"  {Rv:>6.2f}  [ERROR: {e}]")

    if m2_results:
        m2_vals = [m2 for _, m2 in m2_results]
        has_negative = any(m2 < 0 for m2 in m2_vals)
        min_m2_R, min_m2 = min(m2_results, key=lambda x: x[1])
        print(f"\n  min m²(R) = {min_m2:.4e} at R={min_m2_R:.2f}")
        if has_negative:
            print("  → m²(R*) < 0: spontaneous twist generation possible!")
        else:
            print("  → m²(R) > 0 for all R: ω=0 stable (chiral equilibrium maintained)")

    # ------------------------------------------------------------------
    # Step 6: Symbolic verification (exact proof or counterexample for P = 0)
    # ------------------------------------------------------------------
    sym_P_result = None
    if do_symbolic:
        print(f"\n[Step 6] Symbolic computation of P = ⟨R, ∗R⟩")
        t3 = time.time()
        sym_P_result = symbolic_verify(engine)
        t4 = time.time()
        print(f"  Symbolic computation done ({t4-t3:.1f}s)")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print(f"\n{'='*65}")
    print("Nil³ isotropic twist — summary")
    print(f"{'='*65}")
    print(f"  max|P(ω)| (numerical) = {max_P:.4e}  {'→ P ≡ 0' if max_P < 1e-10 else '→ P ≠ 0!'}")

    if sym_P_result is not None:
        p_zero = sym_P_result['P_is_zero']
        print(f"  P (symbolic)          = {'0 [exact]' if p_zero else 'non-zero!'}")

    if m2_results:
        all_positive = all(m2 > 0 for _, m2 in m2_results)
        print(f"  sign of m²(R)         = {'all positive (ω=0 stable)' if all_positive else 'negative exists (twist generation)'}")

    print()
    if max_P < 1e-10:
        print("  Conclusion: P ≡ 0 (same cancellation as S³ occurs in Nil³ as well)")
        print("  → Evidence for 'protection theorem' that P≠0 is impossible for C^3_{ij}-type twist is strengthened")
    else:
        print("  Conclusion: P ≠ 0 confirmed! (violation found in Nil³)")
        print("  → Cancellation in S³ may be a special phenomenon due to SO(3) symmetry")

    return {
        'R_val': R_val, 'eta_val': eta_val,
        'omega_vals': OMEGA_SCAN, 'P_vals': P_vals,
        'max_P': max_P,
        'sym_result': sym_result,
        'sym_P_result': sym_P_result,
        'm2_results': m2_results,
    }


# ===========================================================================
# Entry point
# ===========================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Nil³ isotropic twist Pontryagin scan (App B.7)')
    parser.add_argument('--R-val', type=float, default=1.0,
                        help='Fixed R value (default: 1.0)')
    parser.add_argument('--eta', type=float, default=ETA_DEFAULT,
                        help=f'Torsion coefficient η (default: {ETA_DEFAULT})')
    parser.add_argument('--no-symbolic', action='store_true',
                        help='Skip symbolic computation (faster)')
    parser.add_argument('--all-R', action='store_true',
                        help='Scan all R values in R_SCAN')
    args = parser.parse_args()

    print("=" * 65)
    print("Nil³ isotropic twist — Pontryagin & V_eff scan (App B.7)")
    print("DPPUv2 plan v5.0 Option B")
    print("=" * 65)
    print()
    print("Features of Nil³ twist:")
    print("  S³ twist: C^3_{01}=C^3_{12}=C^3_{20}=+4Lω/r² (3 pairs, positive)")
    print("  Nil³ twist: C^3_{01}=-Lω/R² only (1 pair, negative)")
    print("  → Lower symmetry; cancellation may break")
    print()

    do_symbolic = not args.no_symbolic

    if args.all_R:
        results = []
        for Rv in R_SCAN:
            res = run_phase3_0b(R_val=Rv, eta_val=args.eta,
                                do_symbolic=(Rv == R_SCAN[0]))  # symbolic only for first point
            results.append(res)
    else:
        results = [run_phase3_0b(R_val=args.R_val, eta_val=args.eta,
                                 do_symbolic=do_symbolic)]

    print("\n" + "=" * 65)
    print("Overall summary")
    print("=" * 65)
    print()
    print(f"  {'R':>6}  {'η':>6}  {'max|P|':>14}  {'P=0?':>8}  {'m²>0?':>8}")
    print(f"  {'-'*6}  {'-'*6}  {'-'*14}  {'-'*8}  {'-'*8}")
    for res in results:
        p_zero = res['max_P'] < 1e-10
        m2_pos = all(m2 > 0 for _, m2 in res['m2_results']) if res['m2_results'] else None
        m2_str = 'YES' if m2_pos else ('NO' if m2_pos is False else 'N/A')
        print(f"  {res['R_val']:>6.2f}  {res['eta_val']:>6.2f}  "
              f"{res['max_P']:>14.4e}  {'YES' if p_zero else 'NO':>8}  {m2_str:>8}")


if __name__ == '__main__':
    _t_start = time.time()
    main()
    print(f"\nComputation time: {time.time() - _t_start:.1f}s")
    print()
    print("*Script: scripts/paper03/pontryagin_scan_nil3.py")
