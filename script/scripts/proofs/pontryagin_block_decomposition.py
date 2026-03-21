"""
Algebraic structure analysis of P=0 (App B.2–B.3)
==========================================

Previously established:
  - R_{i3j3} != 0 (the 33 block exists)
  - Yet P equiv 0

Goal: Compute the "block-by-block breakdown" of P = Sigma R_{abcd} (*R)^{abcd}
      and algebraically track which pairs cancel.

Approach:
  1. Decompose P into contributions from SS×SS, SS×S3, SS×33, S3×S3, S3×33, 33×33
  2. Symbolically compute each contribution and verify whether it is zero
  3. Read off the reason for cancellation from index patterns

Author: Muacca
Date: 2026-03-15
"""

import sys
import os
import io
import time

# Windows cp932 workaround
if sys.stdout.encoding and sys.stdout.encoding.lower() not in ('utf-8', 'utf8'):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from dppu.topology.unified import UnifiedEngine, DOFConfig

from sympy import S, simplify, symbols, factor
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant
from dppu.engine.logger import NullLogger
from dppu.utils.levi_civita import levi_civita_4d


# ===========================================================================
# Block classification utilities
# ===========================================================================

def first_pair_type(a, b):
    """Type of the (a,b) pair: 'SS' or 'S3'"""
    n3 = (1 if a == 3 else 0) + (1 if b == 3 else 0)
    if n3 == 0:
        return 'SS'
    elif n3 == 1:
        return 'S3'
    else:
        return '33'


def block_label(a, b, c, d):
    """Block label for R_{abcd}: (type of first pair) x (type of second pair)"""
    return f"{first_pair_type(a,b)}x{first_pair_type(c,d)}"


# ===========================================================================
# Decompose contributions to P by block type
# ===========================================================================

def compute_P_by_blocks(engine_name: str, engine) -> dict:
    """
    Decompose P = (1/2) Sigma_{abcdef} epsilon_{abef} R_{abcd} R_{efcd}
    by block type and compute each contribution.

    Uses the map induced by epsilon_{abef} (dual of the (a,b) pair -> (e,f) pair)
    to track which block combinations contribute to P.
    """
    print(f"\n{'='*65}")
    print(f"  {engine_name}: Block-by-block contribution analysis of P")
    print(f"{'='*65}")

    R_sym = engine.data['riemann_abcd']
    dim = engine.data['dim']

    # P = Sigma_{a,b,c,d} R_{abcd} × (1/2 Sigma_ef epsilon^{abef} R_{efcd})
    # Since epsilon is totally antisymmetric, restrict to a<b, e<f and multiply by 2
    # P = Sigma_{a<b, c,d} R_{abcd} × Sigma_{e<f} epsilon_{abef} R_{efcd}

    # Accumulate contributions by block type
    contributions = {}   # key: "type(ab)x type(ef)" -> SymPy expr

    for a in range(dim):
        for b in range(dim):
            if a >= b:
                continue
            for c in range(dim):
                for d in range(dim):
                    R_abcd = R_sym[a, b, c, d]
                    if R_abcd == S.Zero:
                        continue

                    for e in range(dim):
                        for f in range(dim):
                            if e >= f:
                                continue
                            eps = levi_civita_4d(a, b, e, f)
                            if eps == 0:
                                continue
                            R_efcd = R_sym[e, f, c, d]
                            if R_efcd == S.Zero:
                                continue

                            # Block label for this contribution
                            label_ab = first_pair_type(a, b)
                            label_ef = first_pair_type(e, f)
                            label_cd = first_pair_type(c, d)
                            key = f"({label_ab}→{label_ef}) × ({label_cd})"

                            term = eps * R_abcd * R_efcd
                            contributions[key] = contributions.get(key, S.Zero) + term

    # Simplify and display each contribution
    print(f"\n  [Block contribution breakdown]")
    print(f"  (ab-dual->ef) x cd  :  contribution value")
    print(f"  {'-'*60}")

    P_total = S.Zero
    for key in sorted(contributions.keys()):
        val = simplify(contributions[key])
        contributions[key] = val
        is_zero = (val == S.Zero)
        P_total += val
        print(f"  {key:35s}: {'0' if is_zero else str(val)[:60]}")

    P_total_simplified = simplify(P_total)
    print(f"\n  P total (before applying factor 1/2): {P_total_simplified}")
    print(f"  P (1/2 x total): {simplify(P_total_simplified / 2)}")

    return contributions


# ===========================================================================
# Algebraic structure analysis of the 33 block
# ===========================================================================

def analyze_33_block_structure(engine_name: str, engine):
    """
    Detailed algebraic structure analysis of the 33 block (R_{i3j3}).

    Verifies the hypothesis for S³: R_{i3j3} = (4L²ω²/r⁴)(3δ_{ij} - 1).
    """
    print(f"\n{'='*65}")
    print(f"  {engine_name}: Algebraic structure of the 33 block")
    print(f"{'='*65}")

    R_sym = engine.data['riemann_abcd']
    dim = engine.data['dim']

    # Retrieve R_{i3j3} (i,j in {0,1,2}) as a matrix
    print("\n  R_{i3j3} matrix (i,j in {0,1,2}):")
    R_33_matrix = []
    for i in range(3):
        row = []
        for j in range(3):
            val = R_sym[i, 3, j, 3]
            row.append(val)
        R_33_matrix.append(row)
        row_str = "  ".join(str(v)[:25] for v in row)
        print(f"    i={i}: [{row_str}]")

    # Check ratio of diagonal to off-diagonal entries
    diag_val = R_33_matrix[0][0]
    offdiag_val = R_33_matrix[0][1] if dim > 1 else None

    if offdiag_val is not None and diag_val != S.Zero:
        ratio = simplify(offdiag_val / diag_val) if diag_val != S.Zero else 'N/A'
        print(f"\n  Diagonal entry: {diag_val}")
        print(f"  Off-diagonal entry: {offdiag_val}")
        print(f"  Ratio (offdiag/diag): {ratio}")

    # Compute Sigma_i R_{i3i3} (trace)
    trace_33 = sum(R_33_matrix[i][i] for i in range(3))
    trace_33_simplified = simplify(trace_33)
    print(f"\n  Tr(R_{{i3j3}}) = Σ_i R_{{i3i3}} = {trace_33_simplified}")

    # Sigma_{i,j} R_{i3j3}^2 (squared Frobenius norm)
    frob_sq = sum(simplify(R_33_matrix[i][j]**2)
                  for i in range(3) for j in range(3))
    print(f"  Σ_{{i,j}} R_{{i3j3}}^2 = {simplify(frob_sq)}")

    # Inner product of SS x 33 (candidate for direct contribution to P)
    print(f"\n  P contribution from SS x 33: Sigma_{{i<j,k<l}} epsilon_{{ij,kl}} R_{{ijkl}} R_{{...}}")
    print("  (See output of compute_P_by_blocks for details)")


# ===========================================================================
# Track the action of the Hodge dual on explicit components of R
# ===========================================================================

def trace_hodge_action_on_33(engine_name: str, engine):
    """
    Compute the Hodge dual (*R)_{i3j3} of the 33 block
    and track which block it maps into.
    """
    print(f"\n{'='*65}")
    print(f"  {engine_name}: Computing Hodge(*R)_{{i3j3}}")
    print(f"{'='*65}")

    R_sym = engine.data['riemann_abcd']
    dim = engine.data['dim']

    # (*R)_{abcd} = (1/2) Sigma_{ef} epsilon_{abef} R_{efcd}
    # 33 block: (a,b) = (i,3) or (3,i) — here we compute the (i,3) type

    print("\n  (*R)_{i3j3} = (1/2) Sigma_{ef} epsilon_{i3,ef} R_{efj3}")
    print("  (e,f) with epsilon_{i3,ef} != 0 are spatial pairs (SS type) only")

    for i in range(3):
        for j in range(3):
            dual_val = S.Zero
            for e in range(dim):
                for f in range(dim):
                    eps = levi_civita_4d(i, 3, e, f)
                    if eps == 0:
                        continue
                    R_efj3 = R_sym[e, f, j, 3]
                    dual_val += eps * R_efj3

            dual_val_simplified = simplify(dual_val / 2)
            print(f"    (*R)_{{{i}3{j}3}} = {dual_val_simplified}")


# ===========================================================================
# Decisive test: cancellation only within SS<->S3 vs. independent cancellation of 33
# ===========================================================================

def test_independent_block_cancellation(engine_name: str, engine):
    """
    Decompose P into the following three parts and verify whether each is zero:
      P_A = Sigma_{SS first pair} R_{ab,cd} (*R)^{abcd}
      P_B = Sigma_{S3 first pair} R_{ab,cd} (*R)^{abcd}
      P_C = Sigma_{33 first pair} R_{ab,cd} (*R)^{abcd}
    Check whether each vanishes individually, or whether P_A + P_B + P_C = 0.
    """
    print(f"\n{'='*65}")
    print(f"  {engine_name}: Independent block-by-block verification of P")
    print(f"{'='*65}")

    R_sym = engine.data['riemann_abcd']
    dim = engine.data['dim']

    def compute_dual_R(a, b, c, d):
        """(*R)_{abcd} = (1/2) Sigma_{ef} epsilon_{abef} R_{efcd}"""
        val = S.Zero
        for e in range(dim):
            for f in range(dim):
                eps = levi_civita_4d(a, b, e, f)
                if eps == 0:
                    continue
                val += eps * R_sym[e, f, c, d]
        return val / 2

    P_SS = S.Zero  # contributions where the first pair is of SS type
    P_S3 = S.Zero  # contributions where the first pair is of S3 type
    P_33 = S.Zero  # contributions where the first pair is of 33 type (either a or b equals 3)

    for a in range(dim):
        for b in range(dim):
            pt = first_pair_type(a, b)
            for c in range(dim):
                for d in range(dim):
                    R_abcd = R_sym[a, b, c, d]
                    if R_abcd == S.Zero:
                        continue
                    dualR = compute_dual_R(a, b, c, d)
                    if dualR == S.Zero:
                        continue
                    term = R_abcd * dualR
                    if pt == 'SS':
                        P_SS += term
                    elif pt == 'S3':
                        P_S3 += term
                    else:
                        P_33 += term

    P_SS_s = simplify(P_SS)
    P_S3_s = simplify(P_S3)
    P_33_s = simplify(P_33)
    P_total = simplify(P_SS + P_S3 + P_33)

    print(f"\n  P_SS (first pair SS): {P_SS_s}")
    print(f"  P_S3 (first pair S3): {P_S3_s}")
    print(f"  P_33 (first pair 33): {P_33_s}")
    print(f"  P_total             : {P_total}")

    print(f"\n  Sign relations among contributions:")
    print(f"    P_SS = 0 independently?  {'YES' if P_SS_s == S.Zero else 'NO'}")
    print(f"    P_S3 = 0 independently?  {'YES' if P_S3_s == S.Zero else 'NO'}")
    print(f"    P_33 = 0 independently?  {'YES' if P_33_s == S.Zero else 'NO'}")
    print(f"    P_SS + P_S3 = 0? {'YES' if simplify(P_SS + P_S3) == S.Zero else 'NO'}")
    print(f"    P_SS + P_S3 + P_33 = 0 (total)? "
          f"{'YES' if P_total == S.Zero else 'NO'}")

    return {
        'P_SS': P_SS_s, 'P_S3': P_S3_s, 'P_33': P_33_s, 'P_total': P_total
    }


# ===========================================================================
# Main
# ===========================================================================

def main():
    print("=" * 65)
    print("Algebraic structure analysis of P equiv 0 (App B.2–B.3)")
    print("=" * 65)

    print("\n[Running engines...]")
    _cfg_eng_s3 = DOFConfig.from_engine('TwistS3S1Engine')
    _cfg_eng_s3.torsion_mode = Mode.AX
    _cfg_eng_s3.ny_variant = NyVariant.FULL
    eng_s3 = UnifiedEngine(_cfg_eng_s3, NullLogger())
    eng_s3.run()
    print("  S³ engine complete")

    _cfg_eng_nil3 = DOFConfig.from_engine('TwistNil3S1Engine')
    _cfg_eng_nil3.torsion_mode = Mode.AX
    _cfg_eng_nil3.ny_variant = NyVariant.FULL
    eng_nil3 = UnifiedEngine(_cfg_eng_nil3, NullLogger())
    eng_nil3.run()
    print("  Nil³ engine complete")

    # --- 33 block algebraic structure ---
    analyze_33_block_structure("S³", eng_s3)
    analyze_33_block_structure("Nil³", eng_nil3)

    # --- Action of Hodge dual on the 33 block ---
    trace_hodge_action_on_33("S³", eng_s3)
    trace_hodge_action_on_33("Nil³", eng_nil3)

    # --- Decisive test: independent vanishing by block ---
    print("\n[Decisive test] Contributions to P by first-pair block type")
    result_s3 = test_independent_block_cancellation("S³", eng_s3)
    result_nil3 = test_independent_block_cancellation("Nil³", eng_nil3)

    # --- Summary ---
    print(f"\n{'='*65}")
    print("  Summary: Structure of P equiv 0")
    print(f"{'='*65}")

    for name, result in [("S³", result_s3), ("Nil³", result_nil3)]:
        print(f"\n  [{name}]")
        ss0 = result['P_SS'] == S.Zero
        s30 = result['P_S3'] == S.Zero
        p330 = result['P_33'] == S.Zero
        tot0 = result['P_total'] == S.Zero

        if ss0 and s30 and p330:
            print("    -> All blocks independently zero (strongest result)")
        elif ss0 and s30 and not p330:
            print("    -> SS and S3 independently zero; 33 nonzero but total is zero")
        elif not ss0 and not s30 and tot0:
            from sympy import simplify as sp
            if sp(result['P_SS'] + result['P_S3']) == S.Zero:
                print("    -> P_SS + P_S3 = 0 (SS and S3 cancel each other)")
                if p330:
                    print("      P_33 = 0 holds independently")
                else:
                    print("      P_33 cancels separately via 33x33 etc.")
            else:
                print("    -> Complex cancellation (combination of three or more terms)")
        print(f"    P_SS={result['P_SS']}, "
              f"P_S3={result['P_S3']}, "
              f"P_33={result['P_33']}")


if __name__ == '__main__':
    _t_start = time.time()
    main()
    print(f"\nComputation time: {time.time() - _t_start:.1f}s")
    print()
    print("*Script: scripts/proofs/pontryagin_block_decomposition.py")
