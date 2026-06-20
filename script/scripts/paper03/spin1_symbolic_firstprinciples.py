#!/usr/bin/env python
"""First-principles symbolic derivation of the spin-1 corrected stability for paper03 App E.

This script uses the paper03 dppu engine to derive:
  - G (= -1/2 H_velocity) and H (= -H_static) in the X/Y basis from first principles
    (no hardcoded witness values as inputs)
  - G_XX, G_XY, G_YY, H_XX, H_XY, H_YY symbolically
  - Schur-reduced eigenvalue lambda_phys = H_eff / G_XX = -56*L**2/(3*r**2)
    derived from the actual expression and verified by sp.simplify/sp.factor
  - Degenerate pencil det(H - lambda*G) = 0 finite root, confirmed to match lambda_phys
  - G evaluated at r=3, L=1 and verified against the paper03 witness constants
    (G_oo=236.870435, G_dd=532.958162, G_od=-355.305547) to ensure consistency

Normalization note:
  - This script uses the Ricci-scalar / unnormalized-u normalization
    (u = (-2, 3) unnormalized, so G_XX = (u^T G u) at r=3,L=1 is ~10007.77)
  - App E.3 uses the sqrt(13)-normalized X_k = (-2*omega + 3*delta)/sqrt(13)
    so G_XX(normalized) ~ 769.83 = G_XX(unnormalized) / 13
  - The ratio lambda_phys = H_eff / G_XX is normalization-independent;
    only the individual matrix entries differ by the factor 13.

All PASS/FAIL outcomes are determined by computed conditions, never by hardcoded assertions.

Engine: dppu.action.spin1_gradient (plain sympy, no python-flint/fastsym required).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import sympy as sp


SCRIPT_DIR = Path(__file__).resolve().parent
PAPER_ROOT = SCRIPT_DIR.parents[2]

# engine root
ENGINE_DIR = PAPER_ROOT / "script"
SCRIPTS_PAPER03_DIR = ENGINE_DIR / "scripts" / "paper03"

# output
OUTPUT_DIR = PAPER_ROOT / "data"
SUMMARY_PATH = OUTPUT_DIR / "spin1_symbolic_firstprinciples_summary.json"

# witness values from spin1_auxiliary_exclusion.py (for cross-check only)
WITNESS_G_OO = sp.Rational("236.870435")
WITNESS_G_DD = sp.Rational("532.958162")
WITNESS_G_OD = sp.Rational("-355.305547")

# tolerance for floating-point cross-check
WITNESS_REL_TOL = sp.Rational("1e-4")


def sstr(expr: Any) -> str:
    if isinstance(expr, str):
        return expr
    return sp.sstr(sp.factor(sp.simplify(expr)))


def as_float(expr: Any) -> float:
    return float(sp.N(expr, 30))


def check_zero(expr: sp.Expr) -> bool:
    """Return True if sp.simplify(expr) == 0."""
    return sp.simplify(expr) == 0


def main() -> int:
    import sys

    # Make the paper03 dppu engine importable
    sys.path.insert(0, str(ENGINE_DIR))

    from dppu.action.spin1_gradient import derive_raw_spin1_hessians, derive_xy_coefficients

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("spin1_symbolic_firstprinciples: starting first-principles symbolic derivation")

    # ------------------------------------------------------------------
    # Step 1: derive raw Hessians from the S3 engine (no witness inputs)
    # ------------------------------------------------------------------
    raw = derive_raw_spin1_hessians(axis=0)
    r = raw.r
    L = raw.L

    # Euclidean potential convention: G = -(1/2)*H_velocity, H = -H_static
    G_raw = sp.simplify(-raw.H_velocity / 2)
    H_raw = sp.simplify(-raw.H_static)

    print(f"G_raw =\n{G_raw}")
    print(f"H_raw =\n{H_raw}")

    # ------------------------------------------------------------------
    # Step 2: rotate into X/Y basis via derive_xy_coefficients
    # ------------------------------------------------------------------
    coeffs = derive_xy_coefficients(raw)
    G_xy = sp.simplify(-coeffs.H_velocity_XY / 2)
    H_xy_mat = sp.simplify(-coeffs.H_static_XY)

    g_xx = sp.factor(sp.simplify(G_xy[0, 0]))
    g_xy = sp.factor(sp.simplify(G_xy[0, 1]))
    g_yy = sp.factor(sp.simplify(G_xy[1, 1]))
    h_xx = sp.factor(sp.simplify(H_xy_mat[0, 0]))
    h_xy = sp.factor(sp.simplify(H_xy_mat[0, 1]))
    h_yy = sp.factor(sp.simplify(H_xy_mat[1, 1]))

    print(f"G_XX = {sstr(g_xx)}")
    print(f"G_XY = {sstr(g_xy)}")
    print(f"G_YY = {sstr(g_yy)}")
    print(f"H_XX = {sstr(h_xx)}")
    print(f"H_XY = {sstr(h_xy)}")
    print(f"H_YY = {sstr(h_yy)}")

    # ------------------------------------------------------------------
    # Step 3: Schur elimination and lambda_phys symbolic derivation
    # ------------------------------------------------------------------
    schur_hessian = sp.factor(sp.simplify(h_xx - h_xy**2 / h_yy))
    lambda_phys = sp.factor(sp.simplify(schur_hessian / g_xx))

    print(f"H_eff (Schur) = {sstr(schur_hessian)}")
    print(f"lambda_phys = {sstr(lambda_phys)}")

    # The expected symbolic form is -56*L**2/(3*r**2).
    # We derive this check FROM the computed expression — not from a hardcoded constant.
    expected_form = -56 * L**2 / (3 * r**2)
    lambda_matches_formula = check_zero(sp.simplify(lambda_phys - expected_form))
    print(f"lambda_phys matches -56*L**2/(3*r**2): {lambda_matches_formula}")

    # Negativity check: substitute into computed expression
    lambda_at_r3_l1 = sp.factor(lambda_phys.subs({r: 3, L: 1}))
    lambda_negative_at_sample = bool(lambda_at_r3_l1 < 0)
    print(f"lambda_phys at r=3, L=1: {sstr(lambda_at_r3_l1)} < 0: {lambda_negative_at_sample}")

    # ------------------------------------------------------------------
    # Step 4: generalized pencil det(H - lambda*G) = 0 finite root
    # ------------------------------------------------------------------
    lam = sp.symbols("lambda")
    pencil_det = sp.factor(sp.simplify((H_xy_mat - lam * G_xy).det()))
    finite_roots = sp.solve(sp.Eq(pencil_det, 0), lam)
    finite_roots_factored = [sp.factor(sp.simplify(root)) for root in finite_roots]
    # pick the finite root (smallest absolute value at r=3,L=1)
    def abs_at_sample(expr: sp.Expr) -> float:
        return abs(float(sp.N(expr.subs({r: 3, L: 1}), 30)))

    finite_roots_factored_sorted = sorted(finite_roots_factored, key=abs_at_sample)
    finite_root = finite_roots_factored_sorted[0]
    print(f"finite root of det(H-lambda*G)=0: {sstr(finite_root)}")

    pencil_matches_schur = check_zero(sp.simplify(finite_root - lambda_phys))
    print(f"finite root matches lambda_phys: {pencil_matches_schur}")

    # ------------------------------------------------------------------
    # Step 5: G first-principles vs witness cross-check (r=3, L=1)
    #
    # Normalization difference:
    #   G_Ricci = G from the Ricci-scalar normalization (this script, kappa=1)
    #   G_witness = G from spin1_auxiliary_exclusion.py (per-axis action integral)
    #
    # The witness is computed with the full action integral including the S^3 x S^1
    # volume factor Vol = 2*pi^2*L*r^3.  The action Lagrangian density is R/(2*kappa^2).
    # With kappa=1, the kinetic block accumulates a factor of (Vol/2) = pi^2*L*r^3:
    #
    #   G_witness = pi^2 * L * r^3 * G_Ricci
    #
    # We verify this by computing G_Ricci * (pi^2 * L * r^3) at r=3, L=1
    # and comparing to the witness values.  The small residuals are due to
    # rounding in the witness constants (6 decimal places).
    # ------------------------------------------------------------------
    r_val, L_val = 3, 1
    G_raw_at_sample = G_raw.subs({r: r_val, L: L_val})
    vol_factor = sp.pi**2 * L_val * r_val**3  # pi^2 * L * r^3

    g_oo_scaled = float(sp.N(vol_factor * G_raw_at_sample[0, 0], 30))
    g_dd_scaled = float(sp.N(vol_factor * G_raw_at_sample[1, 1], 30))
    g_od_scaled = float(sp.N(vol_factor * G_raw_at_sample[0, 1], 30))

    # relative error check
    def rel_error(derived: float, witness: float) -> float:
        if witness == 0:
            return abs(derived)
        return abs((derived - witness) / witness)

    g_oo_relerr = rel_error(g_oo_scaled, float(WITNESS_G_OO))
    g_dd_relerr = rel_error(g_dd_scaled, float(WITNESS_G_DD))
    g_od_relerr = rel_error(g_od_scaled, float(WITNESS_G_OD))

    g_oo_match = g_oo_relerr < float(WITNESS_REL_TOL)
    g_dd_match = g_dd_relerr < float(WITNESS_REL_TOL)
    g_od_match = g_od_relerr < float(WITNESS_REL_TOL)

    print(f"G[0,0] * vol_factor = {g_oo_scaled:.6f}, witness = {float(WITNESS_G_OO):.6f}, relerr = {g_oo_relerr:.2e}, match = {g_oo_match}")
    print(f"G[1,1] * vol_factor = {g_dd_scaled:.6f}, witness = {float(WITNESS_G_DD):.6f}, relerr = {g_dd_relerr:.2e}, match = {g_dd_match}")
    print(f"G[0,1] * vol_factor = {g_od_scaled:.6f}, witness = {float(WITNESS_G_OD):.6f}, relerr = {g_od_relerr:.2e}, match = {g_od_match}")
    print(f"Normalization: witness_G = pi^2 * L * r^3 * G_Ricci (volume factor from action integral)")

    # ------------------------------------------------------------------
    # Step 6: Normalization note for G_XX
    # The raw X/Y transform used here is the unit-normalized X_k = (-2*omega + 3*delta)/sqrt(13).
    # So G_XX from derive_xy_coefficients is the normalized Ricci-scalar version.
    # The unnormalized version (u=(-2,3) vector) gives G_XX * 13.
    # The action-integral (witness) version further scales by vol_factor = pi^2 * L * r^3.
    # lambda_phys is independent of all these choices.
    # ------------------------------------------------------------------
    g_xx_at_sample = sp.N(g_xx.subs({r: r_val, L: L_val}), 30)
    g_xx_unnorm_at_sample = 13 * g_xx_at_sample  # unnormalized u=(-2,3) convention (Ricci)
    g_xx_action_norm = float(sp.N(vol_factor * g_xx_at_sample, 30))
    g_xx_action_unnorm = float(sp.N(vol_factor * g_xx_unnorm_at_sample, 30))
    print(f"G_XX (Ricci, sqrt(13)-normalized, r=3,L=1) = {float(g_xx_at_sample):.6f}")
    print(f"G_XX (Ricci, unnormalized u=(-2,3), r=3,L=1) = {float(g_xx_unnorm_at_sample):.6f}")
    print(f"G_XX (action-integral, sqrt(13)-normalized, r=3,L=1) = {g_xx_action_norm:.4f}  [App E.3: ~769.83]")
    print(f"G_XX (action-integral, unnormalized u=(-2,3), r=3,L=1) = {g_xx_action_unnorm:.4f}  [= 13 * 769.83]")
    print("Note: lambda_phys = H_eff/G_XX is normalization-independent; all normalizations give -56/27 at r=3,L=1.")

    # ------------------------------------------------------------------
    # Step 7: G_YY = 0 and H_YY > 0 structural checks
    # ------------------------------------------------------------------
    g_yy_is_zero = check_zero(g_yy)
    h_yy_positive_at_sample = bool(sp.N(h_yy.subs({r: 3, L: 1}), 30) > 0)
    h_xy_nonzero_at_sample = bool(abs(float(sp.N(h_xy.subs({r: 3, L: 1}), 30))) > 1e-6)
    g_xx_positive_at_sample = bool(sp.N(g_xx.subs({r: 3, L: 1}), 30) > 0)

    print(f"G_YY = 0 (exact): {g_yy_is_zero}")
    print(f"H_YY > 0 at r=3,L=1: {h_yy_positive_at_sample}")
    print(f"H_XY != 0 at r=3,L=1: {h_xy_nonzero_at_sample}")
    print(f"G_XX > 0 at r=3,L=1: {g_xx_positive_at_sample}")

    # ------------------------------------------------------------------
    # Aggregate checks
    # ------------------------------------------------------------------
    checks = {
        "g_yy_exactly_zero": g_yy_is_zero,
        "h_yy_positive_at_sample": h_yy_positive_at_sample,
        "h_xy_nonzero_at_sample": h_xy_nonzero_at_sample,
        "g_xx_positive_at_sample": g_xx_positive_at_sample,
        "lambda_phys_matches_formula": lambda_matches_formula,
        "lambda_phys_negative_at_sample": lambda_negative_at_sample,
        "pencil_root_matches_lambda": pencil_matches_schur,
        # G cross-check: Ricci G * (pi^2 * L * r^3) should match witness G within 1e-4 rel
        "g_oo_firstprinciples_matches_witness": g_oo_match,
        "g_dd_firstprinciples_matches_witness": g_dd_match,
        "g_od_firstprinciples_matches_witness": g_od_match,
    }

    all_pass = all(checks.values())
    status_str = "PASS" if all_pass else "FAIL"
    exit_label = (
        "PAPER03_SPIN1_SYMBOLIC_FIRSTPRINCIPLES_PASS"
        if all_pass
        else "PAPER03_SPIN1_SYMBOLIC_FIRSTPRINCIPLES_FAIL"
    )

    summary: dict[str, Any] = {
        "script": str(Path(__file__).resolve()),
        "status": status_str,
        "exit_label": exit_label,
        "derivation_method": "first_principles_via_dppu_action_spin1_gradient",
        "hardcoded_witness_inputs_used": False,
        "normalization": (
            "Ricci-scalar normalization for symbolic derivation. "
            "Witness G is in action-integral normalization: witness_G = pi^2*L*r^3 * G_Ricci. "
            "G_XX (action, sqrt(13)-normalized) ~ 769.83 at r=3,L=1. "
            "G_XX (action, unnormalized u=(-2,3)) ~ 10008 at r=3,L=1. "
            "lambda_phys = H_eff/G_XX is normalization-independent."
        ),
        "G_XX_symbolic": sstr(g_xx),
        "G_XY_symbolic": sstr(g_xy),
        "G_YY_symbolic": sstr(g_yy),
        "H_XX_symbolic": sstr(h_xx),
        "H_XY_symbolic": sstr(h_xy),
        "H_YY_symbolic": sstr(h_yy),
        "H_eff_schur_symbolic": sstr(schur_hessian),
        "lambda_phys_symbolic_derived": sstr(lambda_phys),
        "lambda_phys_expected_form": "-56*L**2/(3*r**2)",
        "lambda_phys_at_r3_l1": sstr(lambda_at_r3_l1),
        "finite_pencil_root_symbolic": sstr(finite_root),
        "vol_normalization_factor": "pi^2 * L * r^3",
        "G_oo_scaled_at_r3_l1": f"{g_oo_scaled:.6f}",
        "G_dd_scaled_at_r3_l1": f"{g_dd_scaled:.6f}",
        "G_od_scaled_at_r3_l1": f"{g_od_scaled:.6f}",
        "witness_G_oo": str(float(WITNESS_G_OO)),
        "witness_G_dd": str(float(WITNESS_G_DD)),
        "witness_G_od": str(float(WITNESS_G_OD)),
        "G_oo_rel_error": f"{g_oo_relerr:.3e}",
        "G_dd_rel_error": f"{g_dd_relerr:.3e}",
        "G_od_rel_error": f"{g_od_relerr:.3e}",
        "G_XX_action_normalized_at_r3_l1": f"{g_xx_action_norm:.4f}",
        "G_XX_action_unnormalized_at_r3_l1": f"{g_xx_action_unnorm:.4f}",
        "checks": checks,
    }

    SUMMARY_PATH.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))

    print(f"status={status_str}")
    print(f"exit_label={exit_label}")
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
