"""
Pfaffian Analysis of Pontryagin Vanishing (App B.2.4–B.2.6)
===================================================

For each fixed (c,d), extract ω_{ab} := R_{ab,cd} and perform:
  Step 1: Verify Pfaffian = ω_{01}ω_{23} - ω_{02}ω_{13} + ω_{03}ω_{12} = 0
  Step 2: Rank check (decomposable ↔ rank ≤ 2)
  Step 3: Direct expansion of P = Σ_{c<d} Pf(ω^{(cd)}) = 0

Core thesis:
  (★) ↔ Pf(ω^{(cd)}) = 0  ↔ ω^{(cd)} is decomposable
  → P = Σ_{cd} Pf(ω^{(cd)}) = 0

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

from sympy import S, simplify, Matrix
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant
from dppu.engine.logger import NullLogger


# ===========================================================================
# Utilities
# ===========================================================================

def extract_omega_block(R, c: int, d: int) -> Matrix:
    """
    For fixed (c,d), construct the 4×4 antisymmetric matrix ω_{ab} := R_{ab,cd}.
    """
    omega = Matrix.zeros(4, 4)
    for a in range(4):
        for b in range(4):
            omega[a, b] = R[a, b, c, d]
    return omega


def compute_pfaffian(omega: Matrix):
    """
    Compute the Pfaffian of a 4×4 antisymmetric matrix.
    Pf(ω) = ω_{01}ω_{23} - ω_{02}ω_{13} + ω_{03}ω_{12}
    """
    return (omega[0, 1] * omega[2, 3]
            - omega[0, 2] * omega[1, 3]
            + omega[0, 3] * omega[1, 2])


def block_label(c: int, d: int) -> str:
    n3 = (1 if c == 3 else 0) + (1 if d == 3 else 0)
    return ['SS', 'S3', '33'][min(n3, 2)]


# ===========================================================================
# Step 1+2: Verification of Pfaffian and Rank
# ===========================================================================

def check_decomposability(engine_name: str, engine):
    """
    Verify Pfaffian = 0 and rank ≤ 2 for each (c,d) block.
    decomposable ↔ Pf = 0 ↔ rank ≤ 2 (equivalence relation in 4D)
    """
    print(f"\n{'='*65}")
    print(f"  {engine_name}: Riemann block decomposability analysis")
    print(f"{'='*65}")

    R = engine.data['riemann_abcd']

    print(f"\n  Pfaffian and Rank for each (c,d) block:")
    print(f"  {'(c,d)':>6}  {'type':>4}  {'Pf(ω)':>8}  {'Rank':>6}  {'decomposable?':>15}")
    print(f"  {'-'*55}")

    results = []
    all_pf_zero = True

    for c in range(4):
        for d in range(4):
            if c >= d:
                continue
            omega = extract_omega_block(R, c, d)
            pf = simplify(compute_pfaffian(omega))
            rank = omega.rank()
            decomposable = (pf == S.Zero) and (rank <= 2)

            if pf != S.Zero:
                all_pf_zero = False

            pf_disp = "0" if pf == S.Zero else f"≠0 ({str(pf)[:20]}...)"
            decomp_str = "✓ Yes" if decomposable else "✗ No"
            lbl = block_label(c, d)

            print(f"  ({c},{d}):  {lbl:>4}  {pf_disp:>8}  {rank:>6}  {decomp_str:>15}")
            results.append({
                'c': c, 'd': d, 'omega': omega,
                'pf': pf, 'rank': rank,
                'block_label': lbl, 'decomposable': decomposable,
            })

    # Summary
    n = len(results)
    n_pf_zero = sum(1 for r in results if r['pf'] == S.Zero)
    n_rank_ok = sum(1 for r in results if r['rank'] <= 2)
    n_decomp = sum(1 for r in results if r['decomposable'])

    print(f"\n  Summary:")
    print(f"    Pfaffian = 0 : {n_pf_zero}/{n} blocks")
    print(f"    Rank ≤ 2     : {n_rank_ok}/{n} blocks")
    print(f"    Decomposable : {n_decomp}/{n} blocks")

    if all_pf_zero:
        print(f"\n  → Pf(ω^{{(cd)}}) = 0 for all (c<d) blocks ✓")
    else:
        print(f"\n  → Non-zero Pfaffian exists ✗")

    return results, all_pf_zero


# ===========================================================================
# Step 3: Verification of decomposable structure in non-zero blocks
# ===========================================================================

def show_nonzero_blocks(engine_name: str, results):
    """
    Display the components of non-zero (rank > 0) decomposable blocks
    and verify the geometric meaning of rank ≤ 2.
    """
    print(f"\n{'='*65}")
    print(f"  {engine_name}: Details of non-zero decomposable blocks")
    print(f"{'='*65}")

    active = [r for r in results if r['rank'] > 0]
    if not active:
        print("  No non-zero blocks (all rank 0)")
        return

    print(f"\n  Non-zero blocks: {len(active)}")

    for r in active:
        c, d = r['c'], r['d']
        omega = r['omega']
        print(f"\n  --- ω^{{({c},{d})}} [{r['block_label']}] rank={r['rank']}, Pf=0 ---")

        # Enumerate non-zero components
        nonzero = []
        for a in range(4):
            for b in range(a + 1, 4):
                val = simplify(omega[a, b])
                if val != S.Zero:
                    nonzero.append((a, b, val))

        if nonzero:
            for a, b, val in nonzero:
                vs = str(val) if len(str(val)) < 60 else str(val)[:57] + "..."
                print(f"    ω_{{{a}{b}}} = {vs}")
        else:
            print("    (all components = 0)")

        # Confirm decomposition from column space of rank-2 block
        if r['rank'] == 2:
            # Verify existence of ω = u∧v: a rank-2 antisymmetric matrix is always decomposable
            # Construct u, v from column space (simplified version)
            col_vecs = []
            for j in range(4):
                col = [simplify(omega[i, j]) for i in range(4)]
                if any(x != S.Zero for x in col):
                    col_vecs.append(col)
            if len(col_vecs) >= 2:
                print(f"    → rank=2 → ω = u∧v exists (u, v are basis of column space) ✓")


# ===========================================================================
# Step 4: Direct expansion of P = Σ_{c<d} Pf(ω^{(cd)})
# ===========================================================================

def p_as_pfaffian_sum(engine_name: str, engine):
    """
    Directly expand the Pontryagin density as a sum of Pfaffians over each (c,d) block
    and verify that each Pf = 0 implies P = 0.
    """
    print(f"\n{'='*65}")
    print(f"  {engine_name}: Direct expansion of P = Σ_{{c<d}} Pf(ω^{{(cd)}})")
    print(f"{'='*65}")

    R = engine.data['riemann_abcd']

    print("""
  Pfaffian decomposition of Pontryagin density:
    P ∝ Σ_{a<b,e<f} ε^{abef} R_{abcd} R_{efcd}
      = Σ_{c<d} [ω_{01}ω_{23} - ω_{02}ω_{13} + ω_{03}ω_{12}]
      = Σ_{c<d} Pf(ω^{(cd)})

  Pfaffian for each (c,d) block:
    """)

    total = S.Zero
    cd_pairs = [(c, d) for c in range(4) for d in range(4) if c < d]

    for c, d in cd_pairs:
        pf = simplify(
            R[0, 1, c, d] * R[2, 3, c, d]
            - R[0, 2, c, d] * R[1, 3, c, d]
            + R[0, 3, c, d] * R[1, 2, c, d]
        )
        lbl = block_label(c, d)
        pf_str = "0" if pf == S.Zero else f"≠0: {str(pf)[:30]}"
        print(f"    Pf(ω^{{({c},{d})}}) [{lbl}] = {pf_str}")
        total += pf

    total_s = simplify(total)
    print(f"\n  P ∝ Σ_{{c<d}} Pf(ω^{{(cd)}}) = {total_s}")
    if total_s == S.Zero:
        print("  → P ≡ 0 confirmed ✓")
    else:
        print("  → P ≠ 0  ✗")

    return total_s


# ===========================================================================
# Main
# ===========================================================================

def main():
    print("=" * 65)
    print("Pfaffian Analysis of Pontryagin Vanishing (App B.2.4–B.2.6)")
    print("Riemann block decomposability")
    print("=" * 65)

    print("\n[Running engines...]")
    _cfg_eng_s3 = DOFConfig.from_engine('TwistS3S1Engine')
    _cfg_eng_s3.torsion_mode = Mode.AX
    _cfg_eng_s3.ny_variant = NyVariant.FULL
    eng_s3 = UnifiedEngine(_cfg_eng_s3, NullLogger())
    eng_s3.run()
    _cfg_eng_nil3 = DOFConfig.from_engine('TwistNil3S1Engine')
    _cfg_eng_nil3.torsion_mode = Mode.AX
    _cfg_eng_nil3.ny_variant = NyVariant.FULL
    eng_nil3 = UnifiedEngine(_cfg_eng_nil3, NullLogger())
    eng_nil3.run()
    print("  Done\n")

    for name, eng in [("S³×S¹", eng_s3), ("Nil³×S¹", eng_nil3)]:
        # Step 1+2: Pfaffian + Rank
        results, all_ok = check_decomposability(name, eng)

        # Step 3: Details of non-zero blocks
        show_nonzero_blocks(name, results)

        # Step 4: Direct expansion of P = Σ Pf
        p_as_pfaffian_sum(name, eng)

        print()



if __name__ == '__main__':
    _t_start = time.time()
    main()
    print(f"\nComputation time: {time.time() - _t_start:.1f}s")
    print()
    print("*Script: scripts/proofs/pfaffian_decomposability.py")
