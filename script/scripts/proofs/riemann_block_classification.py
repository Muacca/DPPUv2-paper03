"""
Riemann Block Structure Analysis (App B.2)
==========================================

Preparation toward the sufficient-condition lemma for P≡0:

1. Symbolically extract non-zero R_{abcd} from TwistS3S1Engine and TwistNil3S1Engine
2. Classify by block type:
   - Type SS: R_{ijkl}  (all spatial indices ∈ {0,1,2})
   - Type S3: R_{ijk3}  (exactly one index equal to 3)
   - Type 33: R_{i3j3}  (two indices equal to 3, diagonal)
   - Type 3X: R_{3ikl} etc. (other patterns containing 3)
3. Verify R_{i3j3} = 0 (Condition 2 of the lemma)
4. Display index patterns for each block in detail

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

from sympy import S, simplify
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant
from dppu.engine.logger import NullLogger


# ===========================================================================
# Block classification utilities
# ===========================================================================

def classify_index(a: int) -> str:
    return 'S' if a < 3 else '3'


def block_type(a: int, b: int, c: int, d: int) -> str:
    """
    Return the block type from four indices.
    SS: all spatial
    S3: exactly one index equal to 3
    33: two indices equal to 3
    3X: three or more indices equal to 3
    """
    n3 = sum(1 for x in (a, b, c, d) if x == 3)
    if n3 == 0:
        return 'SS'
    elif n3 == 1:
        return 'S3'
    elif n3 == 2:
        return '33'
    else:
        return '3X'


def pattern_str(a: int, b: int, c: int, d: int) -> str:
    """Return the index pattern of R_{abcd} as a string."""
    labels = ['0', '1', '2', '3']
    return f"R_{{{labels[a]}{labels[b]}{labels[c]}{labels[d]}}}"


def is_type_i3j3(a: int, b: int, c: int, d: int) -> bool:
    """Return True if the component is of type R_{i3j3} (the block relevant to Condition 2 of the lemma)."""
    return b == 3 and d == 3 and a != 3 and c != 3


def is_type_3i3j(a: int, b: int, c: int, d: int) -> bool:
    """Return True if the component is of type R_{3i3j} (the symmetric counterpart of R_{i3j3})."""
    return a == 3 and c == 3 and b != 3 and d != 3


def analyze_riemann_blocks(engine_name: str, engine) -> dict:
    """
    Extract and classify the R_{abcd} components from the given engine.
    """
    print(f"\n{'='*60}")
    print(f"  {engine_name}")
    print(f"{'='*60}")

    R_sym = engine.data['riemann_abcd']
    dim = engine.data['dim']

    # Collect non-zero components
    nonzero = []
    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                for d in range(dim):
                    expr = R_sym[a, b, c, d]
                    if expr != S.Zero:
                        nonzero.append((a, b, c, d, expr))

    print(f"\n  Non-zero components: {len(nonzero)}")

    # Tally by block type
    blocks = {'SS': [], 'S3': [], '33': [], '3X': []}
    for a, b, c, d, expr in nonzero:
        bt = block_type(a, b, c, d)
        blocks[bt].append((a, b, c, d, expr))

    print(f"\n  Count by block type:")
    for bt, items in blocks.items():
        print(f"    {bt}: {len(items)}", end='')
        if items and bt != 'SS':  # Show details
            patterns = set(f"R_{{{a}{b}{c}{d}}}" for a, b, c, d, _ in items)
            print(f"  → {sorted(patterns)}", end='')
        print()

    # --- Check R_{i3j3} type (lemma Condition 2) ---
    print(f"\n  [Condition 2 check] R_{{i3j3}} type components:")
    i3j3_components = [(a, b, c, d, e) for a, b, c, d, e in nonzero if is_type_i3j3(a, b, c, d)]
    if not i3j3_components:
        print("    -> No non-zero components found  R_{i3j3} = 0 confirmed")
    else:
        print(f"    -> {len(i3j3_components)} components are non-zero! Condition 2 NOT satisfied")
        for a, b, c, d, expr in i3j3_components[:5]:
            print(f"      {pattern_str(a,b,c,d)} = {expr}")

    # --- Also check R_{3i3j} type ---
    type_3i3j = [(a, b, c, d, e) for a, b, c, d, e in nonzero if is_type_3i3j(a, b, c, d)]
    if type_3i3j:
        print(f"    R_{{3i3j}} type: {len(type_3i3j)}")

    # --- SS block details ---
    print(f"\n  [Type SS = R_{{ijkl}}]  (purely spatial block)")
    if blocks['SS']:
        # Show only independent components (a<b, c<d, and (a,b) <= (c,d))
        independent_SS = []
        for a, b, c, d, expr in blocks['SS']:
            if a < b and c < d and (a, b) <= (c, d):
                independent_SS.append((a, b, c, d, expr))
        print(f"    Independent components (a<b, c<d, (a,b)<=( c,d)): {len(independent_SS)}")
        for a, b, c, d, expr in independent_SS:
            print(f"      R_{{{a}{b}{c}{d}}} = {expr}")
    else:
        print("    None")

    # --- S3 block details ---
    print(f"\n  [Type S3 = R_{{ijk3}} etc.]  (mixed block)")
    if blocks['S3']:
        # Organize independent components by the position of index 3
        s3_with_d3 = [(a, b, c, d, e) for a, b, c, d, e in blocks['S3'] if d == 3 and a < b]
        s3_with_c3 = [(a, b, c, d, e) for a, b, c, d, e in blocks['S3'] if c == 3 and a < b]
        s3_with_a3 = [(a, b, c, d, e) for a, b, c, d, e in blocks['S3'] if a == 3 and c < d]
        s3_with_b3 = [(a, b, c, d, e) for a, b, c, d, e in blocks['S3'] if b == 3 and c < d]

        if s3_with_d3:
            print(f"    R_{{ijk3}} type (d=3, a<b): {len(s3_with_d3)}")
            for a, b, c, d, expr in s3_with_d3[:6]:
                print(f"      R_{{{a}{b}{c}{d}}} = {expr}")
        if s3_with_a3:
            print(f"    R_{{3ijk}} type (a=3, c<d): {len(s3_with_a3)}")
            for a, b, c, d, expr in s3_with_a3[:6]:
                print(f"      R_{{{a}{b}{c}{d}}} = {expr}")
    else:
        print("    None")

    # --- 33 block details ---
    print(f"\n  [Type 33]  (block containing two indices equal to 3)")
    if blocks['33']:
        for a, b, c, d, expr in blocks['33']:
            print(f"    R_{{{a}{b}{c}{d}}} = {expr}")
    else:
        print("    None")

    # --- Index pattern summary ---
    print(f"\n  [Index pattern tally]")
    pattern_dict = {}
    for a, b, c, d, _ in nonzero:
        # Pattern indicating which of the four indices equal 3
        pattern = (classify_index(a), classify_index(b),
                   classify_index(c), classify_index(d))
        pattern_dict[pattern] = pattern_dict.get(pattern, 0) + 1

    for pat, cnt in sorted(pattern_dict.items()):
        print(f"    ({pat[0]},{pat[1]},{pat[2]},{pat[3]}): {cnt}")

    return {
        'nonzero': nonzero,
        'blocks': blocks,
        'i3j3_zero': len(i3j3_components) == 0,
    }


# ===========================================================================
# Hodge dual block mapping check
# ===========================================================================

def analyze_hodge_block_mapping():
    """
    Verify which block type the Hodge dual maps to in four-dimensional Euclidean space.

    dual(SS) -> S3 or 33?
    dual(S3) -> SS or 33?
    """
    from dppu.utils.levi_civita import levi_civita_4d

    print(f"\n{'='*60}")
    print("  Hodge Dual block mapping check")
    print(f"{'='*60}")
    print("""
  (*R)_{abcd} = (1/2) ε_{abef} R^{ef}_{cd}
              = (1/2) Σ_{ef} ε_{abef} R_{efcd}

  Below: enumerate the (ab) -> (ef) mapping for non-zero ε_{abef} != 0
    """)

    dim = 4
    print("  Mappings where ε_{abef} != 0:")
    mapping = {}  # (block_ab) -> set of (block_ef)
    for a in range(dim):
        for b in range(dim):
            if a >= b:
                continue
            for e in range(dim):
                for f in range(dim):
                    if e >= f:
                        continue
                    eps = levi_civita_4d(a, b, e, f)
                    if eps != 0:
                        ba = block_type(a, b, 0, 0)[:2]  # First two indices only
                        be = block_type(e, f, 0, 0)[:2]
                        # Simplified labels: SS/S3/33
                        label_ab = 'S3' if (a == 3 or b == 3) else 'SS'
                        label_ef = 'S3' if (e == 3 or f == 3) else 'SS'
                        key = label_ab
                        mapping.setdefault(key, set()).add(label_ef)
                        print(f"    ε_{{{a}{b}{e}{f}}} = {eps:+d}  "
                              f"({label_ab}: ({a},{b})) → ({label_ef}: ({e},{f}))")

    print(f"\n  Mapping summary:")
    for src, dsts in sorted(mapping.items()):
        print(f"    dual({src}) -> {sorted(dsts)}")


# ===========================================================================
# Main
# ===========================================================================

def main():
    print("=" * 60)
    print("Riemann Block Structure Analysis (App B.2)")
    print("Verification of lemma preconditions for P≡0 (Conditions 1-3)")
    print("=" * 60)

    # --- S³ engine ---
    print("\nRunning [S³×S¹ with Isotropic Twist]...")
    _cfg_eng_s3 = DOFConfig.from_engine('TwistS3S1Engine')
    _cfg_eng_s3.torsion_mode = Mode.AX
    _cfg_eng_s3.ny_variant = NyVariant.FULL
    eng_s3 = UnifiedEngine(_cfg_eng_s3, NullLogger())
    eng_s3.run()
    result_s3 = analyze_riemann_blocks("S³×S¹ with Isotropic Twist", eng_s3)

    # --- Nil³ engine ---
    print("\nRunning [Nil³×S¹ with Isotropic Twist]...")
    _cfg_eng_nil3 = DOFConfig.from_engine('TwistNil3S1Engine')
    _cfg_eng_nil3.torsion_mode = Mode.AX
    _cfg_eng_nil3.ny_variant = NyVariant.FULL
    eng_nil3 = UnifiedEngine(_cfg_eng_nil3, NullLogger())
    eng_nil3.run()
    result_nil3 = analyze_riemann_blocks("Nil³×S¹ with Isotropic Twist", eng_nil3)

    # --- Hodge mapping check ---
    analyze_hodge_block_mapping()

    # --- Lemma condition summary ---
    print(f"\n{'='*60}")
    print("  Lemma condition check summary")
    print(f"{'='*60}")

    for name, result in [("S³", result_s3), ("Nil³", result_nil3)]:
        blocks = result['blocks']
        cond1 = len(blocks['SS']) > 0 or len(blocks['S3']) > 0  # Whether restricted to SS/S3 only
        cond1_clean = len(blocks['3X']) == 0  # No 3X blocks
        cond2 = result['i3j3_zero']           # R_{i3j3} = 0
        cond3_check = len(blocks['33']) == 0  # No 33 blocks (stronger condition)

        print(f"\n  [{name}]")
        print(f"    Condition 1 (non-zero only in SS or S3): "
              f"{'✓' if cond1_clean else 'x'}  "
              f"(3X={len(blocks['3X'])}, 33={len(blocks['33'])})")
        print(f"    Condition 2 (R_{{i3j3}} = 0):              "
              f"{'✓' if cond2 else 'x'}")
        print(f"    No 33 blocks (stronger condition):        "
              f"{'✓' if cond3_check else 'x'}  "
              f"({len(blocks['33'])} components)")
        print(f"    -> Lemma preconditions: {'SATISFIED' if (cond1_clean and cond2) else 'NOT SATISFIED'}")

    print()
    print("  Related analysis: Index contraction analysis for P=Σ R_{abcd} (*R)^{abcd}")


if __name__ == '__main__':
    _t_start = time.time()
    main()
    print(f"\nComputation time: {time.time() - _t_start:.1f}s")
    print()
    print("*Script: scripts/proofs/riemann_block_classification.py")
