#!/usr/bin/env python3
"""
Symbolic Block Structure Verification for All Topologies

Purpose:
    Verify that R^{ab}_{cd} has no cd=Mixed components for all three
    topologies: S^3 x S^1, T^3 x S^1, Nil^3 x S^1.

Author: Muacca
Date: 2026-02-07
"""

import sys
from pathlib import Path
import sympy as sp
from sympy import symbols, simplify, S, factor
from sympy.tensor.array import MutableDenseNDimArray
import warnings

# Setup path: scripts/proofs/ -> scripts/ -> _DPPUv2_Phase2/
PHASE2_DIR = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PHASE2_DIR))

import time as _time
_t_start = _time.time()

from dppu.torsion import Mode, NyVariant
from dppu.topology import UnifiedEngine, DOFConfig

def levi_civita_4d(i, j, k, l):
    """4D Levi-Civita symbol with epsilon_{0123} = +1."""
    indices = [i, j, k, l]
    if len(set(indices)) != 4:
        return 0
    inversions = 0
    arr = list(indices)
    for m in range(4):
        for n in range(m + 1, 4):
            if arr[m] > arr[n]:
                inversions += 1
    return 1 if inversions % 2 == 0 else -1

def cd_block(c, d):
    """Return 'S' for spatial, 'M' for mixed."""
    if c == 3 or d == 3:
        return 'M'
    return 'S'

def analyze_topology(engine_name, topo_name, mode, ny_variant):
    """Analyze block structure for one topology."""
    print(f"\n{'=' * 70}")
    print(f"# Topology: {topo_name}")
    print(f"{'=' * 70}\n")

    # Build engine
    print(f"Building {topo_name} engine ({mode.value} + {ny_variant.value} mode)...")
    cfg = DOFConfig.from_engine(engine_name)
    cfg.torsion_mode = mode
    cfg.ny_variant = ny_variant
    engine = UnifiedEngine(cfg)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        engine.run()

    # Get data
    R_a_bcd = engine.data['riemann']
    eta = engine.data['metric_frame']
    dim = engine.data['dim']

    # Build R^{ab}_{cd}
    R_ab_cd = MutableDenseNDimArray.zeros(dim, dim, dim, dim)

    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                for d in range(dim):
                    val = S.Zero
                    for e in range(dim):
                        val += eta[b, e] * R_a_bcd[a, e, c, d]
                    if val != S.Zero:
                        R_ab_cd[a, b, c, d] = simplify(val)

    # Count components by block
    blocks = {
        ('S', 'S'): 0, ('S', 'M'): 0,
        ('M', 'S'): 0, ('M', 'M'): 0,
    }
    mixed_components = []
    spatial_components = []

    for a in range(dim):
        for b in range(a+1, dim):
            for c in range(dim):
                for d in range(c+1, dim):
                    component = R_ab_cd[a, b, c, d]
                    if component != S.Zero:
                        ab_blk = 'M' if (a == 3 or b == 3) else 'S'
                        cd_blk = cd_block(c, d)
                        blocks[(ab_blk, cd_blk)] += 1

                        if cd_blk == 'M':
                            mixed_components.append((a, b, c, d, component))
                        else:
                            spatial_components.append((a, b, c, d, component))

    # Report
    print(f"\nBlock structure (cd indices):")
    print(f"  cd=Spatial: {blocks[('S','S')] + blocks[('M','S')]} components")
    print(f"  cd=Mixed:   {blocks[('S','M')] + blocks[('M','M')]} components")

    if len(mixed_components) == 0:
        print(f"\n  VERIFIED: cd=Mixed block is EMPTY")
    else:
        print(f"\n  WARNING: cd=Mixed block has {len(mixed_components)} components:")
        for a, b, c, d, comp in mixed_components:
            print(f"    R^{{{a}{b}}}_{{{c}{d}}} = {comp}")

    # Compute P symbolically
    R_star_ab_cd = MutableDenseNDimArray.zeros(dim, dim, dim, dim)
    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                for d in range(dim):
                    val = S.Zero
                    for e in range(dim):
                        for f in range(dim):
                            eps = levi_civita_4d(c, d, e, f)
                            if eps != 0:
                                val += sp.Rational(1, 2) * eps * R_ab_cd[a, b, e, f]
                    if val != S.Zero:
                        R_star_ab_cd[a, b, c, d] = simplify(val)

    P = S.Zero
    for a in range(dim):
        for b in range(a+1, dim):
            for c in range(dim):
                for d in range(c+1, dim):
                    term = R_ab_cd[a, b, c, d] * R_star_ab_cd[a, b, c, d]
                    P += term

    P_simplified = simplify(P)

    print(f"\n  P = <R, *R> = {P_simplified}")

    return {
        'topology': topo_name,
        'cd_spatial': blocks[('S','S')] + blocks[('M','S')],
        'cd_mixed': blocks[('S','M')] + blocks[('M','M')],
        'P': P_simplified,
        'verified': P_simplified == S.Zero and len(mixed_components) == 0
    }

# ============================================================
# Main
# ============================================================

print("=" * 70)
print("# Symbolic Block Structure Verification: All Topologies")
print("=" * 70)

results = []

# Test all topologies
topologies = [
    ('S3S1Engine',   "S^3 x S^1"),
    ('T3S1Engine',   "T^3 x S^1"),
    ('Nil3S1Engine', "Nil^3 x S^1"),
]

for engine_name, topo_name in topologies:
    result = analyze_topology(engine_name, topo_name, Mode.MX, NyVariant.FULL)
    results.append(result)

# Summary
print("\n")
print("=" * 70)
print("# SUMMARY")
print("=" * 70)
print()

print("| Topology      | cd=Spatial | cd=Mixed | P = <R,*R> | Verified |")
print("|---------------|------------|----------|------------|----------|")
for r in results:
    verified = "YES" if r['verified'] else "NO"
    P_str = "0" if r['P'] == S.Zero else str(r['P'])[:20]
    print(f"| {r['topology']:<13} | {r['cd_spatial']:>10} | {r['cd_mixed']:>8} | {P_str:>10} | {verified:>8} |")

print()

all_verified = all(r['verified'] for r in results)
if all_verified:
    print("=" * 70)
    print("THEOREM VERIFIED FOR ALL TOPOLOGIES:")
    print("  P = <R, *R> = 0 is an algebraic identity for M^3 x S^1 with")
    print("  minisuperspace ansatz and EC connection, regardless of M^3.")
    print("=" * 70)
else:
    print("WARNING: Not all topologies verified. Check results above.")

print()

print(f"\nComputation time: {_time.time() - _t_start:.1f}s")
print()
print("*Script: scripts/proofs/symbolic_block_all_topologies.py")
