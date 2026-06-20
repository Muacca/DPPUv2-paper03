"""
Tracing the algebraic reason why P_SS = 0 (App B.2.4–B.2.6)
======================================

Previously established: P_SS = P_S3 = 0 hold independently.

Here we expand the concrete product structure of P_SS
and trace which components cancel each other.

P_SS = Σ_{(a,b)=SS, all c,d} R_{abcd} (*R)^{abcd}
     = 2 × Σ_{c,d} [R_{01cd}R_{23cd} - R_{02cd}R_{13cd} + R_{12cd}R_{03cd}]

Further decomposed by (c,d) block:
    P_SS^{SS}  : c,d ∈ {0,1,2} only
    P_SS^{S3}  : c,d of S3 type

Author: Claude (2026-03-01)
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

from dppu.topology import FiberMode, TopologyType, UnifiedEngine, make_config

from sympy import S, simplify
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant
from dppu.utils.logger import NullLogger
from dppu.utils.epsilon import epsilon_4d


def expand_P_SS(engine_name: str, engine):
    """
    Expand P_SS by (c,d) block and display in detail.
    """
    print(f"\n{'='*65}")
    print(f"  {engine_name}: Expansion of P_SS")
    print(f"{'='*65}")

    R = engine.data['riemann_abcd']
    dim = engine.data['dim']

    # P_SS = Σ_{c,d} Σ_{(a,b)=SS, (e,f)=S3} ε_{abef} R_{abcd} R_{efcd}
    # Explicit expansion: ε_{01,23}=1, ε_{02,13}=-1, ε_{12,03}=1
    # (*R)^{01cd} = R_{23cd}
    # (*R)^{02cd} = -R_{13cd}   [ε_{0213}=-1]
    # (*R)^{12cd} = R_{03cd}

    SS_pairs = [(0,1), (0,2), (1,2)]  # Independent SS pairs
    SS_duals = [(2,3), (1,3), (0,3)]  # Corresponding duals (S3 pairs)
    SS_signs = [+1, -1, +1]           # Signs of ε

    print("\n  Hodge dual correspondence:")
    for (a,b), (e,f), s in zip(SS_pairs, SS_duals, SS_signs):
        print(f"    ε_{{{a}{b}{e}{f}}} = {s:+d}  →  (*R)^{{{a}{b},cd}} = {s:+d} R_{{{e}{f},cd}}")

    # Compute contributions separated by (c,d) into SS / S3 / 33
    print(f"\n  Breakdown of P_SS by (c,d):")
    print(f"  {'(c,d)':>8}  {'type':>6}  {'contribution':>50}")
    print(f"  {'-'*70}")

    total_SS_SS = S.Zero
    total_SS_S3 = S.Zero

    for c in range(dim):
        for d in range(dim):
            # Contribution for this (c,d)
            term = S.Zero
            for (a,b), (e,f), s in zip(SS_pairs, SS_duals, SS_signs):
                # Both orientations: factor of 2 for (a,b) and (b,a)
                R_abcd = R[a, b, c, d]
                R_efcd = R[e, f, c, d]
                term += 2 * s * R_abcd * R_efcd

            term_s = simplify(term)
            if term_s == S.Zero:
                continue

            n3 = (1 if c == 3 else 0) + (1 if d == 3 else 0)
            cd_type = 'SS' if n3 == 0 else ('S3' if n3 == 1 else '33')

            print(f"  ({c},{d}): {cd_type:>6}  {str(term_s)[:55]}")

            if cd_type == 'SS':
                total_SS_SS += term_s
            elif cd_type == 'S3':
                total_SS_S3 += term_s

    print(f"\n  P_SS^{{SS}} (c,d spatial pairs): {simplify(total_SS_SS)}")
    print(f"  P_SS^{{S3}} (c,d mixed pairs):   {simplify(total_SS_S3)}")
    print(f"  P_SS total:                      {simplify(total_SS_SS + total_SS_S3)}")

    return {
        'P_SS_SS': simplify(total_SS_SS),
        'P_SS_S3': simplify(total_SS_S3),
    }


def expand_P_S3(engine_name: str, engine):
    """
    Expand P_S3 and display in detail.
    P_S3 = Σ_{(a,b)=S3, all c,d} R_{abcd} (*R)^{abcd}
    """
    print(f"\n{'='*65}")
    print(f"  {engine_name}: Expansion of P_S3")
    print(f"{'='*65}")

    R = engine.data['riemann_abcd']
    dim = engine.data['dim']

    # S3 pairs: (0,3),(1,3),(2,3) and their symmetric counterparts (3,0),(3,1),(3,2)
    # Hodge: ε_{03,12}=1, ε_{13,02}=-1, ε_{23,01}=1 (etc.)
    S3_pairs = [(0,3), (1,3), (2,3)]
    S3_duals = [(1,2), (0,2), (0,1)]  # To SS pairs
    S3_signs = [+1, -1, +1]

    print("\n  S3→SS dual correspondence:")
    for (a,b), (e,f), s in zip(S3_pairs, S3_duals, S3_signs):
        # Verify the value of ε_{abef}
        eps = epsilon_4d(a, b, e, f)
        print(f"    ε_{{{a}{b}{e}{f}}} = {eps:+d}  →  (*R)^{{{a}{b},cd}} ∝ {eps:+d} R_{{{e}{f},cd}}")

    # Breakdown of P_S3 by (c,d)
    print(f"\n  Breakdown of P_S3 by (c,d):")
    total = S.Zero

    for c in range(dim):
        for d in range(dim):
            term = S.Zero
            for (a,b), (e,f), s in zip(S3_pairs, S3_duals, S3_signs):
                eps_ab_ef = epsilon_4d(a, b, e, f)
                # Both (a,b) and (b,a)
                R_abcd = R[a, b, c, d]
                R_efcd = R[e, f, c, d]
                R_bacd = R[b, a, c, d]   # = -R[a,b,c,d]
                R_fecd = R[f, e, c, d]   # = -R[e,f,c,d]
                eps_ba_fe = epsilon_4d(b, a, f, e)  # = same sign as ab_ef

                term += eps_ab_ef * R_abcd * R_efcd
                term += eps_ba_fe * R_bacd * R_fecd   # = eps * (-R) * (-R) = eps * R * R = same (both orientations)

            term_s = simplify(term)
            if term_s == S.Zero:
                continue

            n3 = (1 if c == 3 else 0) + (1 if d == 3 else 0)
            cd_type = 'SS' if n3 == 0 else ('S3' if n3 == 1 else '33')
            print(f"  ({c},{d}): {cd_type:>6}  {str(term_s)[:55]}")
            total += term_s

    print(f"\n  P_S3 total: {simplify(total)}")


def check_ss_ss_cancellation(engine_name: str, engine):
    """
    Check how P_SS^{SS} = Σ_{c,d=SS} [R_{01cd}R_{23cd} - R_{02cd}R_{13cd} + R_{12cd}R_{03cd}]
    cancels within individual (c,d) pairs.
    """
    print(f"\n{'='*65}")
    print(f"  {engine_name}: Cancellation structure of P_SS^{{SS}}")
    print(f"{'='*65}")

    R = engine.data['riemann_abcd']

    SS_pairs_cd = [(0,1), (0,2), (1,2)]

    print("\n  Contribution for each (c,d)=SS:")
    total = S.Zero
    for (c,d) in SS_pairs_cd:
        term = 2 * (R[0,1,c,d]*R[2,3,c,d] - R[0,2,c,d]*R[1,3,c,d] + R[1,2,c,d]*R[0,3,c,d])
        # (c,d) and (d,c): R[a,b,d,c] = -R[a,b,c,d] → yields the same sign
        term += 2 * (R[0,1,d,c]*R[2,3,d,c] - R[0,2,d,c]*R[1,3,d,c] + R[1,2,d,c]*R[0,3,d,c])
        t = simplify(term)
        print(f"    ({c},{d}) + ({d},{c}): {t}")  # (c,d) and (d,c): R[a,b,d,c] = -R[a,b,c,d] -> same sign
        total += t

    print(f"\n  P_SS^{{SS}} = {simplify(total)}")

    # Also verify the symmetry of the SS block
    print(f"\n  R SS block component check (verify R_{{01}} and R_{{23}} for c,d=SS):")
    for (c,d) in [(0,1),(0,2),(1,2)]:
        r01 = R[0,1,c,d]
        r23 = R[2,3,c,d]
        prod = simplify(r01 * r23)
        print(f"    R_{{01,{c}{d}}} = {r01}, R_{{23,{c}{d}}} = {r23}")


def check_r_symmetry(engine_name: str, engine):
    """
    Verify the symmetry of R (second pair symmetry R_{abcd} = R_{cdab}).
    When this holds, the exchange symmetry P_SS = P_S3 is established.
    """
    print(f"\n{'='*65}")
    print(f"  {engine_name}: Symmetry check of R")
    print(f"{'='*65}")

    R = engine.data['riemann_abcd']
    dim = engine.data['dim']

    print("\n  [Check: R_{abcd} = R_{cdab}] (pair symmetry)")
    violations = []
    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                for d in range(dim):
                    diff = simplify(R[a,b,c,d] - R[c,d,a,b])
                    if diff != S.Zero:
                        violations.append((a,b,c,d,diff))

    if not violations:
        print("    → R_{abcd} = R_{cdab} holds for all components ✓")
    else:
        print(f"    → {len(violations)} violation(s) (first 5):")
        for a,b,c,d,diff in violations[:5]:
            print(f"      R_{{{a}{b}{c}{d}}} - R_{{{c}{d}{a}{b}}} = {diff}")

    print("\n  [Check: R_{abcd} = R_{bacd} antisymmetry]")
    viol2 = []
    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                for d in range(dim):
                    diff = simplify(R[a,b,c,d] + R[b,a,c,d])
                    if diff != S.Zero:
                        viol2.append((a,b,c,d,diff))

    if not viol2:
        print("    → R_{abcd} = -R_{bacd} holds for all components ✓")
    else:
        print(f"    → {len(viol2)} violation(s)")


# ===========================================================================
# Main
# ===========================================================================

def main():
    print("=" * 65)
    print("Tracing the algebraic reason why P_SS = 0 (App B.2.4–B.2.6)")
    print("=" * 65)

    print("\n[Running engines...]")
    _cfg_eng_s3 = make_config(
        TopologyType.S3,
        fiber_mode=FiberMode.TWIST,
        isotropic_twist=True,
        torsion_mode=Mode.AX,
        ny_variant=NyVariant.FULL,
    )
    eng_s3 = UnifiedEngine(_cfg_eng_s3, NullLogger())
    eng_s3.run()
    _cfg_eng_nil3 = make_config(
        TopologyType.NIL3,
        fiber_mode=FiberMode.TWIST,
        isotropic_twist=True,
        torsion_mode=Mode.AX,
        ny_variant=NyVariant.FULL,
    )
    eng_nil3 = UnifiedEngine(_cfg_eng_nil3, NullLogger())
    eng_nil3.run()
    print("  Done")

    for name, eng in [("S³", eng_s3), ("Nil³", eng_nil3)]:
        # Check symmetry of R
        check_r_symmetry(name, eng)

        # Expand P_SS
        res = expand_P_SS(name, eng)

        # Cancellation structure of P_SS^{SS}
        check_ss_ss_cancellation(name, eng)

        # Expand P_S3
        expand_P_S3(name, eng)



if __name__ == '__main__':
    _t_start = time.time()
    main()
    print(f"\nComputation time: {time.time() - _t_start:.1f}s")
    print()
    print("*Script: scripts/proofs/plucker_cancellation_proof.py")
