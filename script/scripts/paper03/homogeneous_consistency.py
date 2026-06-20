"""
Homogeneous-ansatz consistency check (App F.7).

Scope:
    - Check only the 3D homogeneous frame sector.
    - Use constant frame deformations e^i = r M^i_j sigma^j.
    - Verify antisymmetry and Jacobi for transformed structure constants.
    - Do not claim that the full 4D geometry forms a new Lie group manifold.

Author: Muacca
Date: 2026-03-15
"""

from __future__ import annotations

import io
import sys
import time

import sympy as sp


def _configure_utf8_stdio() -> None:
    for name in ("stdout", "stderr"):
        stream = getattr(sys, name)
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
            continue
        except Exception:
            pass

        buffer = getattr(stream, "buffer", None)
        if buffer is not None:
            setattr(sys, name, io.TextIOWrapper(buffer, encoding="utf-8", errors="replace"))


def build_spin2_basis() -> list[sp.Matrix]:
    sqrt2 = sp.sqrt(2)
    sqrt6 = sp.sqrt(6)
    return [
        sp.diag(1, 1, -2) / sqrt6,
        sp.diag(1, -1, 0) / sqrt2,
        sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]) / sqrt2,
        sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]) / sqrt2,
        sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]) / sqrt2,
    ]


def make_s3_structure_constants() -> sp.MutableDenseNDimArray:
    c = sp.MutableDenseNDimArray([0] * 27, (3, 3, 3))
    for i, j, k, value in (
        (0, 1, 2, 1),
        (1, 2, 0, 1),
        (2, 0, 1, 1),
        (0, 2, 1, -1),
        (2, 1, 0, -1),
        (1, 0, 2, -1),
    ):
        c[i, j, k] = value
    return c


def make_nil3_structure_constants() -> sp.MutableDenseNDimArray:
    c = sp.MutableDenseNDimArray([0] * 27, (3, 3, 3))
    c[2, 0, 1] = 1
    c[2, 1, 0] = -1
    return c


def build_sample_shape_matrix(coeffs: list[sp.Expr]) -> sp.Matrix:
    shear_basis = build_spin2_basis()
    s_mat = sp.zeros(3)
    for coeff, basis in zip(coeffs, shear_basis):
        s_mat += coeff * basis
    return sp.eye(3) + s_mat


def transform_structure_constants(
    c_old: sp.MutableDenseNDimArray, m_mat: sp.Matrix
) -> sp.MutableDenseNDimArray:
    m_inv = m_mat.inv()
    c_new = sp.MutableDenseNDimArray([0] * 27, (3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                value = 0
                for a in range(3):
                    for b in range(3):
                        for c in range(3):
                            value += (
                                m_mat[i, a]
                                * m_inv[b, j]
                                * m_inv[c, k]
                                * c_old[a, b, c]
                            )
                c_new[i, j, k] = sp.simplify(sp.cancel(value))
    return c_new


def antisymmetry_residuals(c_arr: sp.MutableDenseNDimArray) -> list[tuple[int, int, int, sp.Expr]]:
    residuals: list[tuple[int, int, int, sp.Expr]] = []
    for i in range(3):
        for j in range(3):
            for k in range(3):
                value = sp.simplify(c_arr[i, j, k] + c_arr[i, k, j])
                if value != 0:
                    residuals.append((i, j, k, value))
    return residuals


def jacobi_residuals(
    c_arr: sp.MutableDenseNDimArray,
) -> list[tuple[int, int, int, int, sp.Expr]]:
    residuals: list[tuple[int, int, int, int, sp.Expr]] = []
    for d in range(3):
        for a in range(3):
            for b in range(3):
                for c in range(3):
                    value = 0
                    for e in range(3):
                        value += (
                            c_arr[d, a, e] * c_arr[e, b, c]
                            + c_arr[d, b, e] * c_arr[e, c, a]
                            + c_arr[d, c, e] * c_arr[e, a, b]
                        )
                    value = sp.simplify(sp.cancel(value))
                    if value != 0:
                        residuals.append((d, a, b, c, value))
    return residuals


def main() -> None:
    _configure_utf8_stdio()

    start = time.perf_counter()
    samples = [
        (
            "sample-1",
            [
                sp.Rational(1, 7),
                sp.Rational(2, 11),
                sp.Rational(-1, 5),
                sp.Rational(3, 13),
                sp.Rational(-2, 17),
            ],
        ),
        (
            "sample-2",
            [
                sp.Rational(-1, 9),
                sp.Rational(1, 6),
                sp.Rational(2, 15),
                sp.Rational(-1, 8),
                sp.Rational(1, 10),
            ],
        ),
    ]

    topologies = [
        ("S3", make_s3_structure_constants()),
        ("Nil3", make_nil3_structure_constants()),
    ]

    print("=" * 72)
    print("Homogeneous ansatz consistency (App F.7)")
    print("=" * 72)
    print()
    print("[1] General constant-frame formula")
    print("  e^i = r M^i_j sigma^j with constant invertible M")
    print("  C'^i_jk = M^i_a (M^-1)^b_j (M^-1)^c_k C^a_bc")
    print()
    print("[2] Exact component checks for representative full-shear samples")

    for sample_name, coeffs in samples:
        m_mat = build_sample_shape_matrix(coeffs)
        det_m = sp.simplify(m_mat.det())
        print()
        print(f"  {sample_name}: det(M) = {sp.sstr(det_m)}")
        print(f"    coeffs = {[sp.sstr(v) for v in coeffs]}")

        for topo_name, c_old in topologies:
            c_new = transform_structure_constants(c_old, m_mat)
            anti = antisymmetry_residuals(c_new)
            jac = jacobi_residuals(c_new)
            print(
                f"    {topo_name}: antisymmetry residuals = {len(anti)}, "
                f"Jacobi residuals = {len(jac)}"
            )
            assert not anti
            assert not jac

    elapsed = time.perf_counter() - start
    print()
    print(f"\nComputation time: {elapsed:.1f}s")
    print()
    print("*Script: scripts/paper03/homogeneous_consistency.py")


if __name__ == "__main__":
    import os as _os, sys as _sys
    _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..'))
    from dppu.utils.tee_logger import setup_log, teardown_log
    setup_log(__file__, log_dir='output')
    main()
    teardown_log()
