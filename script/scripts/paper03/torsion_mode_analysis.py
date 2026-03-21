"""
MX/VT/AX mode NyVariant analysis — Search for vacuum and spontaneous twist generation (App B.7)
========================================================================

Task 1: Verify initial topological charge at ω=0
Task 2: Symbolic expansion of MX mode V_eff and check cross terms
Task 3: Search for m² < 0 (spontaneous twist generation)

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

import time

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from dppu.topology.unified import UnifiedEngine, DOFConfig

from sympy import (S, simplify, symbols, Matrix, factor, expand,
                   collect, sqrt, Rational, pprint, oo)
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant
from dppu.engine.logger import NullLogger


def compute_pfaffian_sum(R):
    """Compute P ∝ Σ_{c<d} Pf(ω^{(cd)})."""
    P = S.Zero
    for c in range(4):
        for d in range(4):
            if c < d:
                pf = (R[0, 1, c, d] * R[2, 3, c, d]
                      - R[0, 2, c, d] * R[1, 3, c, d]
                      + R[0, 3, c, d] * R[1, 2, c, d])
                P += simplify(pf)
    return simplify(P)


# ===========================================================================
# Task 1: Initial topological charge at ω=0
# ===========================================================================

def task1_p_at_omega_zero():
    """
    Verify P at ω₁=ω₂=ω₃=0 (no twist) in MX mode.
    Check whether "the vacuum is chirally asymmetric from the start".
    """
    print("=" * 68)
    print("  Task 1: Verify initial topological charge at ω=0")
    print("=" * 68)

    # --- S³×S¹ ---
    print("\n  [S³×S¹ MX mode] Running engine...")
    t0 = time.time()
    _cfg_eng_s3 = DOFConfig.from_engine('VectorTwistS3S1Engine')
    _cfg_eng_s3.torsion_mode = Mode.MX
    _cfg_eng_s3.ny_variant = NyVariant.FULL
    eng_s3 = UnifiedEngine(_cfg_eng_s3, NullLogger())
    eng_s3.run()
    print(f"  Done ({time.time()-t0:.1f}s)")

    R_s3 = eng_s3.data['riemann_abcd']
    params_s3 = eng_s3.data['params']
    omega1 = params_s3['omega1']
    omega2 = params_s3['omega2']
    omega3 = params_s3['omega3']
    zero_subs = {omega1: S.Zero, omega2: S.Zero, omega3: S.Zero}

    # General P expression
    P_s3 = compute_pfaffian_sum(R_s3)
    print(f"\n  P(ω₁,ω₂,ω₃) [S³] =")
    print(f"    {P_s3}")

    # Substitute ω=0
    P_s3_zero = simplify(P_s3.subs(zero_subs))
    print(f"\n  P(ω=0) [S³] =")
    print(f"    {P_s3_zero}")

    if P_s3_zero == S.Zero:
        print("  → P(ω=0) = 0 (vacuum without twist is chirally symmetric)")
        s3_vacuum_chiral = False
    else:
        print("  → P(ω=0) ≠ 0 (vacuum is chirally asymmetric from the start)")
        s3_vacuum_chiral = True

    # --- Nil³×S¹ ---
    print("\n  [Nil³×S¹ MX mode] Running engine...")
    t0 = time.time()
    _cfg_eng_nil3 = DOFConfig.from_engine('TwistNil3S1Engine')
    _cfg_eng_nil3.torsion_mode = Mode.MX
    _cfg_eng_nil3.ny_variant = NyVariant.FULL
    eng_nil3 = UnifiedEngine(_cfg_eng_nil3, NullLogger())
    eng_nil3.run()
    print(f"  Done ({time.time()-t0:.1f}s)")

    R_nil3 = eng_nil3.data['riemann_abcd']
    params_nil3 = eng_nil3.data['params']
    omega_nil3 = params_nil3.get('omega', None)

    P_nil3 = compute_pfaffian_sum(R_nil3)
    print(f"\n  P(ω) [Nil³] =")
    print(f"    {P_nil3}")

    if omega_nil3 is not None:
        P_nil3_zero = simplify(P_nil3.subs(omega_nil3, S.Zero))
    else:
        P_nil3_zero = P_nil3
    print(f"\n  P(ω=0) [Nil³] =")
    print(f"    {P_nil3_zero}")

    if P_nil3_zero == S.Zero:
        print("  → P(ω=0) = 0 (Nil³: vacuum without twist is chirally symmetric)")
        nil3_vacuum_chiral = False
    else:
        print("  → P(ω=0) ≠ 0 (Nil³: vacuum is chirally asymmetric from the start)")
        nil3_vacuum_chiral = True

    return eng_s3, P_s3, P_s3_zero, eng_nil3, P_nil3, P_nil3_zero


# ===========================================================================
# Task 2: Check cross terms in V_eff
# ===========================================================================

def task2_veff_cross_terms(eng_s3):
    """
    Symbolically obtain V_eff(r, η, V, ω₁, ω₂, ω₃) in MX mode and analyze its structure.
    1. Can it be fully expressed via I₂ = ω₁²+ω₂²+ω₃²?
    2. Are there cross terms between ω sector (I₂) and torsion sector (V, η)?
    """
    print("\n" + "=" * 68)
    print("  Task 2: Symbolic expansion of V_eff and check cross terms")
    print("=" * 68)

    V_eff = eng_s3.data.get('potential')
    if V_eff is None:
        print("  Cannot obtain V_eff. data keys:", list(eng_s3.data.keys()))
        return None

    params = eng_s3.data['params']
    omega1 = params['omega1']
    omega2 = params['omega2']
    omega3 = params['omega3']
    eta = params['eta']
    V = params['V']
    r = params['r']
    L = params['L']
    kappa = params['kappa']

    # Expand and display V_eff
    V_s = simplify(V_eff)
    print(f"\n  V_eff (before expansion) =")
    print(f"    {V_s}")

    # Replace I₂ = ω₁²+ω₂²+ω₃² and check
    I2 = symbols('I2', real=True, nonnegative=True)
    w = sqrt(I2 / 3)
    V_iso_check = simplify(V_eff.subs({omega1: w, omega2: w, omega3: w}))
    print(f"\n  V(√(I₂/3), ...) =")
    print(f"    {V_iso_check}")

    # Test whether full replacement ω₁²+ω₂²+ω₃² → I₂ is possible
    # Test with ω₂=ω₃=0, ω₁²=I₂
    V_i2_test = simplify(V_eff.subs({omega2: S.Zero, omega3: S.Zero,
                                      omega1: sqrt(I2)}))
    print(f"\n  V(√I₂, 0, 0) =")
    print(f"    {V_i2_test}")

    # If both are equal, V = f(I₂)
    diff_check = simplify(V_iso_check - V_i2_test)
    print(f"\n  V(√(I₂/3),...) - V(√I₂,0,0) = {diff_check}")
    if diff_check == S.Zero:
        print("  → V_eff = f(I₂) holds (fully factorable by I₂)")
        i2_factored = True
    else:
        print("  → V_eff = f(I₂) does not hold (not factorable by I₂)")
        i2_factored = False

    # Analysis of cross terms: expand V_eff as polynomial in ω₁², ω₂², ω₃²
    print(f"\n  --- Analysis of cross terms ---")

    # Value at ω=0 (zeroth order term)
    zero_subs = {omega1: S.Zero, omega2: S.Zero, omega3: S.Zero}
    V0 = simplify(V_s.subs(zero_subs))
    print(f"\n  V_eff(ω=0) =")
    print(f"    {V0}")

    # Coefficient of I₂ (linear term) = m²/2
    # ∂V/∂(I₂) |_{I₂=0} via ω₁ (chain rule: ∂/∂I₂ = (1/2ω₁) ∂/∂ω₁)
    # More directly: 1st order coefficient of I₂ in V_eff(I₂) - V_eff(0)
    dV_dI2_at0 = simplify(V_i2_test.diff(I2).subs(I2, S.Zero))
    print(f"\n  ∂V_eff/∂I₂ |_{{I₂=0}} =")
    print(f"    {dV_dI2_at0}")

    # Check if the coefficient of I₂ contains V, η
    print(f"\n  Symbols in ∂V/∂I₂(ω=0): {dV_dI2_at0.free_symbols}")

    # Second order coefficient of I₂
    d2V_dI2sq_at0 = simplify(V_i2_test.diff(I2, 2).subs(I2, S.Zero))
    print(f"\n  ∂²V_eff/∂I₂² |_{{I₂=0}} =")
    print(f"    {d2V_dI2sq_at0}")

    # Determine presence of cross terms
    has_V_cross = V in dV_dI2_at0.free_symbols
    has_eta_cross = eta in dV_dI2_at0.free_symbols
    print(f"\n  Cross term check (in coefficient of ∂V/∂I₂):")
    print(f"    Contribution of V (vector trace): {'present' if has_V_cross else 'absent'}")
    print(f"    Contribution of η (axial):        {'present' if has_eta_cross else 'absent'}")

    return V_s, V0, dV_dI2_at0, i2_factored


# ===========================================================================
# Task 3: Hessian and search for sign of m²
# ===========================================================================

def task3_hessian_and_stability(eng_s3):
    """
    Compute the Hessian (mass term m²) at ω=0 and search for conditions where m² < 0.
    """
    print("\n" + "=" * 68)
    print("  Task 3: Compute m²(r, V, η) and search for m² < 0")
    print("=" * 68)

    V_eff = eng_s3.data.get('potential')
    if V_eff is None:
        print("  Cannot obtain V_eff.")
        return None

    params = eng_s3.data['params']
    omega1 = params['omega1']
    omega2 = params['omega2']
    omega3 = params['omega3']
    eta = params['eta']
    V = params['V']
    r = params['r']
    L = params['L']
    kappa = params['kappa']

    zero_subs = {omega1: S.Zero, omega2: S.Zero, omega3: S.Zero}

    # Hessian H_{ij} = ∂²V_eff/∂ωᵢ∂ωⱼ |_{ω=0}
    print(f"\n  Computing Hessian H_{{ij}} = ∂²V/∂ωᵢ∂ωⱼ |_{{ω=0}}...")
    omegas = [omega1, omega2, omega3]
    t0 = time.time()
    H = Matrix.zeros(3, 3)
    for i, oi in enumerate(omegas):
        for j, oj in enumerate(omegas):
            if j < i:
                H[i, j] = H[j, i]
            else:
                h_ij = simplify(V_eff.diff(oi).diff(oj).subs(zero_subs))
                H[i, j] = h_ij
    print(f"  Computation time: {time.time()-t0:.1f}s")

    print(f"\n  H =")
    for i in range(3):
        row = [str(H[i, j]) for j in range(3)]
        print(f"    [{', '.join(row)}]")

    # Check diagonality and isotropy
    off_diag_zero = all(
        H[i, j] == S.Zero
        for i in range(3) for j in range(3) if i != j
    )
    print(f"\n  H is diagonal: {'yes' if off_diag_zero else 'no'}")

    m2 = simplify(H[0, 0])
    print(f"\n  m²(r, V, η) = H[0,0] =")
    print(f"    {m2}")
    print(f"\n  Symbols in m²: {m2.free_symbols}")

    # Analyze V, η dependence of m²
    m2_at_V0 = simplify(m2.subs(V, S.Zero))
    m2_at_eta0 = simplify(m2.subs(eta, S.Zero))
    print(f"\n  m²(V=0, η, r) [AX-like] = {m2_at_V0}")
    print(f"  m²(V, η=0, r) [VT-like] = {m2_at_eta0}")

    # Condition for m² = 0 (boundary of sign change)
    print(f"\n  --- Search for m² < 0 ---")

    # Expand m²
    m2_expanded = expand(m2)
    print(f"\n  m² (expanded) = {m2_expanded}")

    # Try to collect m² in terms of V²
    # m² = A(r,L,κ) + B(r,L,κ) * V² + C(r,L,κ) * η² + ...
    # Check if positive at V=0, becomes negative as V increases

    try:
        m2_collected_V = collect(m2_expanded, V)
        print(f"\n  m² (collected in V) = {m2_collected_V}")
    except Exception:
        m2_collected_V = m2_expanded

    try:
        m2_collected_eta = collect(m2_expanded, eta)
        print(f"\n  m² (collected in η) = {m2_collected_eta}")
    except Exception:
        pass

    # Numerical scan: check sign of m² over ranges of r, V, η
    print(f"\n  --- Numerical parameter scan ---")
    print(f"  L=1, κ=1, r ∈ [0.1, 5], V ∈ [0.1, 10], η ∈ [-5, 5]")
    print(f"  α=0 (baseline)")

    import numpy as np
    from sympy import lambdify

    alpha = params.get('alpha', S.Zero)

    # Symbol list for m²
    sym_list = [r, V, eta, L, kappa]
    # Add alpha if it's a symbol
    if hasattr(alpha, 'is_Symbol') and alpha.is_Symbol:
        sym_list.append(alpha)

    try:
        m2_func = lambdify(sym_list, m2, modules='numpy')

        # Scan
        r_vals = np.linspace(0.1, 5.0, 20)
        V_vals = np.linspace(0.1, 10.0, 20)
        eta_vals = np.linspace(0.1, 5.0, 20)

        L_val = 1.0
        kappa_val = 1.0
        alpha_val = 0.0

        min_m2 = float('inf')
        min_params_found = None
        negative_count = 0
        total_count = 0

        for r_v in r_vals:
            for V_v in V_vals:
                for eta_v in eta_vals:
                    total_count += 1
                    try:
                        if hasattr(alpha, 'is_Symbol') and alpha.is_Symbol:
                            m2_val = float(m2_func(r_v, V_v, eta_v, L_val, kappa_val, alpha_val))
                        else:
                            m2_val = float(m2_func(r_v, V_v, eta_v, L_val, kappa_val))
                        if m2_val < min_m2:
                            min_m2 = m2_val
                            min_params_found = (r_v, V_v, eta_v)
                        if m2_val < 0:
                            negative_count += 1
                    except Exception:
                        pass

        print(f"\n  Scan results ({total_count} points):")
        print(f"    Minimum m²: {min_m2:.4f}")
        print(f"    Parameters at minimum: r={min_params_found[0]:.2f}, V={min_params_found[1]:.2f}, η={min_params_found[2]:.2f}")
        print(f"    Points with m² < 0: {negative_count} / {total_count}")

        if negative_count > 0:
            print(f"\n  → m² < 0 exists (spontaneous twist generation possible)")
            # More detailed scan (wider V range)
            print(f"\n  Detailed scan in high-V region (V ∈ [1, 50]):")
            V_vals_wide = np.linspace(1.0, 50.0, 30)
            r_vals_fine = np.linspace(0.1, 2.0, 15)
            for r_v in r_vals_fine[:5]:  # representative values only
                for V_v in V_vals_wide[:10]:
                    for eta_v in [0.1, 1.0, 2.0]:
                        try:
                            if hasattr(alpha, 'is_Symbol') and alpha.is_Symbol:
                                m2_val = float(m2_func(r_v, V_v, eta_v, L_val, kappa_val, alpha_val))
                            else:
                                m2_val = float(m2_func(r_v, V_v, eta_v, L_val, kappa_val))
                            if m2_val < 0:
                                print(f"    m² = {m2_val:.4f} at r={r_v:.2f}, V={V_v:.2f}, η={eta_v:.1f}")
                        except Exception:
                            pass
        else:
            print(f"\n  → m² > 0 in scan range (stable in this range)")
            print(f"     Minimum m² {min_m2:.4f} > 0")

            # Additional scan for α < 0
            if hasattr(alpha, 'is_Symbol') and alpha.is_Symbol:
                print(f"\n  Additional scan for α < 0:")
                alpha_vals = [-0.1, -0.5, -1.0, -5.0]
                for alpha_v in alpha_vals:
                    m2_min_alpha = float('inf')
                    for r_v in r_vals:
                        for V_v in V_vals:
                            for eta_v in eta_vals:
                                try:
                                    m2_val = float(m2_func(r_v, V_v, eta_v, L_val, kappa_val, alpha_v))
                                    if m2_val < m2_min_alpha:
                                        m2_min_alpha = m2_val
                                except Exception:
                                    pass
                    print(f"    α={alpha_v}: minimum m² = {m2_min_alpha:.4f}")

    except Exception as e:
        print(f"  Numerical scan error: {e}")
        print(f"  Sign analysis of m² with symbolic expression:")
        print(f"  m²(V=0) = {simplify(m2.subs(V, S.Zero))}")

    return m2


# ===========================================================================
# Main
# ===========================================================================

def main():
    print("=" * 68)
    print("  MX/VT/AX mode NyVariant analysis: vacuum dynamics (App B.7)")
    print("=" * 68)

    t_total = time.time()

    # Task 1
    eng_s3, P_s3, P_s3_zero, eng_nil3, P_nil3, P_nil3_zero = task1_p_at_omega_zero()

    # Task 2
    V_s, V0, dV_dI2, i2_factored = task2_veff_cross_terms(eng_s3)

    # Task 3
    m2 = task3_hessian_and_stability(eng_s3)

    # Final summary
    elapsed = time.time() - t_total
    print(f"\n{'=' * 68}")
    print(f"  Final summary")
    print(f"{'=' * 68}")
    print(f"\n  Task 1: Check P(ω=0)")
    print(f"    S³:  P(ω=0) = {P_s3_zero}")
    print(f"    Nil³: P(ω=0) = {P_nil3_zero}")

    if V_s is not None:
        print(f"\n  Task 2: V_eff structure")
        print(f"    Fully factorable by I₂: {'yes' if i2_factored else 'no (cross terms present)'}")
        print(f"    Symbols of ∂V_eff/∂I₂(ω=0): {dV_dI2.free_symbols if dV_dI2 else 'not obtained'}")

    if m2 is not None:
        print(f"\n  Task 3: m²(r, V, η)")
        print(f"    m² = {m2}")

    print(f"\nComputation time: {elapsed:.1f}s")
    print()
    print("*Script: scripts/paper03/torsion_mode_analysis.py")


if __name__ == '__main__':
    main()
