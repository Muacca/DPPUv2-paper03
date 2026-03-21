"""
Additional Verification: Symbolic Tracking of the "Pfaffian Sum Cancellation"
Mechanism in VT Mode (Proposal 1) (App B.3)
=========================================================================

Problem:
    AX mode: each block Pf = 0 individually -> P = 0 (cancellation by decomposability)
    VT mode: each block Pf != 0, yet the sum Sigma Pf = 0 -> P = 0
    MX mode: each block Pf != 0, sum != 0 -> P != 0

Goal:
    Symbolically trace why the "fine cancellation mechanism" in VT mode holds.

    Track the sign relations of Pfaffians and the structure of V terms.
    -> Can add a lemma to paper03 as a "VT-specific separate mechanism".

PLAN:
    1. Run VectorTwistS3S1Engine (VT)
    2. Obtain the Pfaffian of each (c,d) block symbolically
    3. Factor/collect each to confirm simplified form
    4. Identify the algebraic structure that makes the sum zero (sign-flip pairs, etc.)
    5. Perform the same analysis for Nil3

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

from sympy import S, simplify, Matrix, factor, expand, collect, symbols
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant
from dppu.engine.logger import NullLogger


def get_pfaffian_blocks(R):
    """
    Obtain the Pfaffian of each (c,d) block symbolically.
    Also returns the simplified form after factor/collect.
    """
    results = {}
    for c in range(4):
        for d in range(4):
            if c < d:
                omega_mat = Matrix.zeros(4, 4)
                for a in range(4):
                    for b in range(4):
                        omega_mat[a, b] = R[a, b, c, d]

                pf_raw = (
                    omega_mat[0, 1] * omega_mat[2, 3]
                    - omega_mat[0, 2] * omega_mat[1, 3]
                    + omega_mat[0, 3] * omega_mat[1, 2]
                )
                pf = simplify(pf_raw)
                rank = omega_mat.rank()

                try:
                    pf_factored = factor(pf)
                except Exception:
                    pf_factored = pf

                results[(c, d)] = {
                    'pf': pf,
                    'pf_factored': pf_factored,
                    'rank': rank,
                    'matrix': omega_mat,
                }
    return results


def analyze_vt_cancellation(topo_name, engine_name:str):
    """
    Analyze the Pfaffian cancellation mechanism in VT mode for the specified topology.
    """
    print(f"\n  {'=' * 64}")
    print(f"  {topo_name} — VT mode cancellation mechanism analysis")
    print(f"  {'=' * 64}")

    t0 = time.time()
    _cfg_eng = DOFConfig.from_engine(engine_name)
    _cfg_eng.torsion_mode = Mode.VT
    _cfg_eng.ny_variant = NyVariant.FULL
    eng = UnifiedEngine(_cfg_eng, NullLogger())
    eng.run()
    elapsed = time.time() - t0
    print(f"  Engine completed ({elapsed:.1f}s)")

    R = eng.data['riemann_abcd']
    blocks = get_pfaffian_blocks(R)

    # List all Pfaffians
    print(f"\n  Pfaffian of each (c,d) block (VT mode, factored):")
    print(f"  {'(c,d)':>6}  {'Rank':>5}  {'Pf':>10}  {'detailed'}")
    print(f"  {'-' * 68}")

    pf_list = []
    for (c, d), info in blocks.items():
        pf = info['pf']
        pf_f = info['pf_factored']
        rank = info['rank']
        n3 = (1 if c == 3 else 0) + (1 if d == 3 else 0)
        lbl = ['SS', 'S3', '33'][min(n3, 2)]
        is_zero = (pf == S.Zero)
        pf_str = "0" if is_zero else "≠0"
        print(f"  ({c},{d})[{lbl}]: rank={rank:>2}, Pf={pf_str}")
        if not is_zero:
            print(f"    Pf_factored = {pf_f}")
        pf_list.append(((c, d), pf))

    # Verify the sum
    P_total = simplify(sum(pf for _, pf in pf_list))
    print(f"\n  Total P = Sigma Pf = {P_total}")
    is_cancelled = (P_total == S.Zero)
    print(f"  -> Cancellation confirmed: {'✓ P = 0 (cancellation holds)' if is_cancelled else '✗ P != 0'}")

    if not is_cancelled:
        return

    # Analysis of the cancellation structure
    print(f"\n  --- Analysis of cancellation structure ---")

    nonzero_pfs = [(cd, pf) for cd, pf in pf_list if pf != S.Zero]
    if not nonzero_pfs:
        print("  All Pf = 0 (AX mode-type cancellation, no cancellation mechanism)")
        return

    print(f"  Number of nonzero Pfaffians: {len(nonzero_pfs)}")

    # Verify pair by pair (are there sign-flip pairs?)
    print(f"\n  Search for sign-flip pairs:")
    found_pairs = []
    checked = set()
    for i, (cd1, pf1) in enumerate(nonzero_pfs):
        for j, (cd2, pf2) in enumerate(nonzero_pfs):
            if i < j and (i, j) not in checked:
                checked.add((i, j))
                diff = simplify(pf1 + pf2)
                if diff == S.Zero:
                    print(f"    {cd1} and {cd2}: Pf1 + Pf2 = 0 -> sign-flip pair ✓")
                    found_pairs.append((cd1, cd2, pf1))

    if not found_pairs:
        print("  No simple sign-flip pairs found (more complex cancellation)")

        # Cancellation from 3 or more combinations? Examine partial sums
        print(f"\n  Search for partial sums (cumulative from the start):")
        running = S.Zero
        for k, (cd, pf) in enumerate(nonzero_pfs):
            running = simplify(running + pf)
            print(f"    After adding block {cd}: cumulative sum = {running if running != S.Zero else '0'}")

    # Decompose each Pf into individual terms and search for common factors
    print(f"\n  Expansion and common factor analysis of each Pf:")
    for cd, pf in nonzero_pfs[:6]:  # Display at most 6 blocks
        pf_expanded = expand(pf)
        print(f"    Pf{cd} = {pf_expanded}")

    # Grouping: sum of SS blocks, sum of S3 blocks
    ss_pfs = [pf for (c, d), pf in pf_list
              if c < 3 and d < 3 and (c, d) != (0, 0)]
    s3_pfs = [pf for (c, d), pf in pf_list if c < 3 and d == 3]

    ss_sum = simplify(sum(ss_pfs)) if ss_pfs else S.Zero
    s3_sum = simplify(sum(s3_pfs)) if s3_pfs else S.Zero

    print(f"\n  Partial sums by block type:")
    print(f"    Sum of SS blocks ((c<3,d<3)) = {ss_sum}")
    print(f"    Sum of S3 blocks ((c<3,d=3)) = {s3_sum}")
    print(f"    SS + S3 = {simplify(ss_sum + s3_sum)}")

    if simplify(ss_sum + s3_sum) == S.Zero:
        print(f"  -> Cancellation is a mutual annihilation of 'sum of SS blocks' and 'sum of S3 blocks' ✓")
    elif ss_sum == S.Zero and s3_sum == S.Zero:
        print(f"  -> Each type independently zero (stronger cancellation) ✓")


def analyze_mx_no_cancellation(topo_name, engine_name:str):
    """
    Confirm that no cancellation occurs in MX mode (as a contrast).
    """
    print(f"\n  --- Contrast: {topo_name} MX mode (confirming no cancellation) ---")

    t0 = time.time()
    _cfg_eng = DOFConfig.from_engine(engine_name)
    _cfg_eng.torsion_mode = Mode.MX
    _cfg_eng.ny_variant = NyVariant.FULL
    eng = UnifiedEngine(_cfg_eng, NullLogger())
    eng.run()
    elapsed = time.time() - t0

    R = eng.data['riemann_abcd']
    blocks = get_pfaffian_blocks(R)

    nonzero_count = sum(1 for info in blocks.values() if info['pf'] != S.Zero)
    P_total = simplify(sum(info['pf'] for info in blocks.values()))

    print(f"  ({elapsed:.1f}s) Non-zero Pf blocks: {nonzero_count}/6")

    if P_total == S.Zero:
        print(f"  P = 0 (unexpected result)")
    else:
        print(f"  P != 0 -> no cancellation ✓ (characteristic of MX mode)")


def main():
    print("=" * 68)
    print("  Additional Verification: Symbolic Tracking of VT Mode Cancellation Mechanism (Proposal 1) (App B.3)")
    print("  (20260301-41 additional_verification.md -> Proposal 1)")
    print("=" * 68)

    t_total = time.time()

    # VT analysis for S3xS1
    analyze_vt_cancellation("S3xS1 (VectorTwist)", 'VectorTwistS3S1Engine')

    # VT analysis for Nil3xS1
    analyze_vt_cancellation("Nil3xS1 (TwistNil3)", 'TwistNil3S1Engine')

    # Contrast: MX mode (no cancellation)
    analyze_mx_no_cancellation("S3xS1", 'VectorTwistS3S1Engine')

    elapsed = time.time() - t_total
    print(f"\nComputation time: {elapsed:.1f}s")
    print()
    print("*Script: scripts/proofs/vt_mode_cancellation.py")


if __name__ == '__main__':
    main()
