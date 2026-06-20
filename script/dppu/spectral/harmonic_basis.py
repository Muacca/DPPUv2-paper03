"""Closed / harmonic left-invariant 1-form basis engine.

Pure compute module — no I/O, no logging, no side effects on import.

Harmonic left-invariant 1-form basis engine for homogeneous Bianchi geometries.

Public entry points
-------------------
structure_constants(n)
closure_matrix_from_c(constants)
basis_matrix(nullspace)
span_equal(left, right)
co_closed_audit(constants, k_form_basis)

Additional helpers also available:
  epsilon, closure_matrix, matrix_rows, matrix_inline,
  nonzero_structure_constants, rotation_xz, rotation_xy,
  has_off_diagonal, wedge, d_form, d_sigma_forms,
  hodge_one_form, co_closed_scalar, audit_case,
  generic_co_closed_check, generic_identity_check
"""

from __future__ import annotations

from typing import Any, Dict, List, Sequence, Tuple

import sympy as sp


# ---------------------------------------------------------------------------
# Constants (needed by closure_matrix_from_c and d_sigma_forms)
# ---------------------------------------------------------------------------

PAIR_ORDER: Sequence[Tuple[int, int]] = ((0, 1), (0, 2), (1, 2))


# ---------------------------------------------------------------------------
# Pure compute functions
# ---------------------------------------------------------------------------

def epsilon(i: int, j: int, k: int) -> int:
    if len({i, j, k}) < 3:
        return 0
    perm = [i, j, k]
    inversions = sum(1 for a in range(3) for b in range(a + 1, 3) if perm[a] > perm[b])
    return -1 if inversions % 2 else 1


def structure_constants(n: sp.Matrix) -> Dict[Tuple[int, int, int], sp.Expr]:
    constants: Dict[Tuple[int, int, int], sp.Expr] = {}
    for i in range(3):
        for j in range(3):
            for k in range(3):
                value = sp.Integer(0)
                for l in range(3):
                    value += epsilon(j, k, l) * n[l, i]
                constants[(i, j, k)] = sp.simplify(value)
    return constants


def closure_matrix_from_c(constants: Dict[Tuple[int, int, int], sp.Expr]) -> sp.Matrix:
    rows = []
    for j, k in PAIR_ORDER:
        rows.append([sp.simplify(constants[(i, j, k)]) for i in range(3)])
    return sp.Matrix(rows)


def closure_matrix(n: sp.Matrix) -> sp.Matrix:
    return closure_matrix_from_c(structure_constants(n))


def basis_matrix(nullspace: Sequence[sp.Matrix]) -> sp.Matrix:
    if not nullspace:
        return sp.zeros(3, 0)
    return sp.Matrix.hstack(*nullspace)


def span_equal(left: sp.Matrix, right: sp.Matrix) -> bool:
    rank_left = left.rank()
    rank_right = right.rank()
    if rank_left != rank_right:
        return False
    if left.cols == 0 and right.cols == 0:
        return True
    combined = sp.Matrix.hstack(left, right)
    return combined.rank() == rank_left


def matrix_rows(matrix: sp.Matrix) -> List[List[str]]:
    return [[str(sp.simplify(matrix[i, j])) for j in range(matrix.cols)] for i in range(matrix.rows)]


def matrix_inline(matrix: sp.Matrix) -> str:
    if matrix.cols == 0:
        return "[]"
    rows = []
    for i in range(matrix.rows):
        rows.append("[" + ",".join(str(sp.simplify(matrix[i, j])) for j in range(matrix.cols)) + "]")
    return "[" + "; ".join(rows) + "]"


def nonzero_structure_constants(constants: Dict[Tuple[int, int, int], sp.Expr]) -> List[str]:
    rows = []
    for (i, j, k), value in sorted(constants.items()):
        if sp.simplify(value) != 0:
            rows.append(f"C^{i + 1}_{j + 1}{k + 1}={sp.simplify(value)}")
    return rows


def rotation_xz() -> sp.Matrix:
    return sp.Matrix(
        [
            [sp.Rational(3, 5), 0, sp.Rational(4, 5)],
            [0, 1, 0],
            [-sp.Rational(4, 5), 0, sp.Rational(3, 5)],
        ]
    )


def rotation_xy() -> sp.Matrix:
    return sp.Matrix(
        [
            [sp.Rational(5, 13), sp.Rational(12, 13), 0],
            [-sp.Rational(12, 13), sp.Rational(5, 13), 0],
            [0, 0, 1],
        ]
    )


def has_off_diagonal(matrix: sp.Matrix) -> bool:
    return any(sp.simplify(matrix[i, j]) != 0 for i in range(3) for j in range(3) if i != j)


def wedge(left: Dict[Tuple[int, ...], sp.Expr], right: Dict[Tuple[int, ...], sp.Expr]) -> Dict[Tuple[int, ...], sp.Expr]:
    result: Dict[Tuple[int, ...], sp.Expr] = {}
    for left_key, left_coeff in left.items():
        for right_key, right_coeff in right.items():
            combined = left_key + right_key
            if len(set(combined)) != len(combined):
                continue
            inversions = sum(1 for i in range(len(combined)) for j in range(i + 1, len(combined)) if combined[i] > combined[j])
            sign = -1 if inversions % 2 else 1
            key = tuple(sorted(combined))
            result[key] = sp.simplify(result.get(key, 0) + sign * left_coeff * right_coeff)
    return {key: sp.simplify(value) for key, value in result.items() if sp.simplify(value) != 0}


def d_form(form: Dict[Tuple[int, ...], sp.Expr], d_sigma: Sequence[Dict[Tuple[int, ...], sp.Expr]]) -> Dict[Tuple[int, ...], sp.Expr]:
    result: Dict[Tuple[int, ...], sp.Expr] = {}
    for key, coeff in form.items():
        if not key:
            continue
        for pos, basis_idx in enumerate(key):
            prefix = {key[:pos]: sp.Integer(1)}
            suffix = {key[pos + 1 :]: sp.Integer(1)}
            term = wedge(wedge(prefix, d_sigma[basis_idx]), suffix)
            sign = -1 if pos % 2 else 1
            for term_key, term_coeff in term.items():
                result[term_key] = sp.simplify(result.get(term_key, 0) + sign * coeff * term_coeff)
    return {key: sp.simplify(value) for key, value in result.items() if sp.simplify(value) != 0}


def d_sigma_forms(constants: Dict[Tuple[int, int, int], sp.Expr]) -> List[Dict[Tuple[int, int], sp.Expr]]:
    forms: List[Dict[Tuple[int, int], sp.Expr]] = []
    for i in range(3):
        two_form: Dict[Tuple[int, int], sp.Expr] = {}
        for j, k in PAIR_ORDER:
            coeff = sp.simplify(-constants[(i, j, k)])
            if coeff != 0:
                two_form[(j, k)] = coeff
        forms.append(two_form)
    return forms


def hodge_one_form(coeffs: Sequence[sp.Expr]) -> Dict[Tuple[int, int], sp.Expr]:
    # Orientation sigma1 wedge sigma2 wedge sigma3 and Euclidean spatial metric.
    return {
        (1, 2): sp.simplify(coeffs[0]),
        (0, 2): sp.simplify(-coeffs[1]),
        (0, 1): sp.simplify(coeffs[2]),
    }


def co_closed_scalar(constants: Dict[Tuple[int, int, int], sp.Expr], coeffs: Sequence[sp.Expr]) -> sp.Expr:
    d_star_alpha = d_form(hodge_one_form(coeffs), d_sigma_forms(constants))
    # delta alpha = -*d*alpha.  For a 3-form, this is the negative volume
    # coefficient.  Only vanishing matters for this audit.
    return sp.simplify(-d_star_alpha.get((0, 1, 2), sp.Integer(0)))


def co_closed_audit(constants: Dict[Tuple[int, int, int], sp.Expr], k_form_basis: sp.Matrix) -> Dict[str, Any]:
    a1, a2, a3 = sp.symbols("a1 a2 a3")
    generic_delta = co_closed_scalar(constants, [a1, a2, a3])
    basis_results = []
    for col in range(k_form_basis.cols):
        vector = k_form_basis[:, col]
        delta_value = co_closed_scalar(constants, list(vector))
        basis_results.append({"basis_index": col, "delta_alpha": str(sp.simplify(delta_value)), "is_zero": bool(sp.simplify(delta_value) == 0)})
    all_basis_zero = all(row["is_zero"] for row in basis_results)
    return {
        "generic_delta_alpha": str(sp.simplify(generic_delta)),
        "basis_delta_results": basis_results,
        "all_K_form_basis_co_closed": bool(all_basis_zero),
        "status": "PASS" if all_basis_zero else "FAIL",
    }


def audit_case(topology: str, bianchi: str, n: sp.Matrix, case_id: str) -> Dict[str, Any]:
    c = structure_constants(n)
    m_closed = closure_matrix_from_c(c)
    k_form_basis = basis_matrix(m_closed.nullspace())
    ker_n_basis = basis_matrix(n.nullspace())
    co_closed = co_closed_audit(c, k_form_basis)
    dim_k_form = k_form_basis.cols
    dim_ker_n = ker_n_basis.cols
    span_match = span_equal(k_form_basis, ker_n_basis)
    dim_match = dim_k_form == dim_ker_n
    rank_n = n.rank()
    rank_closed = m_closed.rank()
    status = "PASS" if dim_match and span_match and co_closed["status"] == "PASS" else "FAIL"
    return {
        "topology": topology,
        "bianchi_type": bianchi,
        "case_id": case_id,
        "n_matrix": matrix_rows(n),
        "closure_matrix": matrix_rows(m_closed),
        "nonzero_structure_constants": nonzero_structure_constants(c),
        "rank_n": int(rank_n),
        "rank_closed_form_operator": int(rank_closed),
        "dim_ker_n": int(dim_ker_n),
        "dim_K_form": int(dim_k_form),
        "K_form_basis": matrix_rows(k_form_basis),
        "ker_n_basis": matrix_rows(ker_n_basis),
        "co_closed_check": co_closed,
        "harmonic_basis_status": "PASS" if span_match and co_closed["status"] == "PASS" else "FAIL",
        "span_match": bool(span_match),
        "dimension_match": bool(dim_match),
        "non_diagonal_n": bool(has_off_diagonal(n)),
        "status": status,
    }


def generic_co_closed_check() -> Dict[str, Any]:
    n11, n12, n13, n22, n23, n33, a1, a2, a3 = sp.symbols("n11 n12 n13 n22 n23 n33 a1 a2 a3")
    n = sp.Matrix([[n11, n12, n13], [n12, n22, n23], [n13, n23, n33]])
    constants = structure_constants(n)
    delta_value = co_closed_scalar(constants, [a1, a2, a3])
    is_zero = sp.simplify(delta_value) == 0
    return {
        "description": "For symmetric class-A n, every left-invariant 1-form is co-closed in the orthonormal invariant frame.",
        "delta_alpha": str(sp.simplify(delta_value)),
        "is_zero": bool(is_zero),
        "status": "PASS" if is_zero else "FAIL",
    }


def generic_identity_check() -> Dict[str, Any]:
    n11, n12, n13, n22, n23, n33 = sp.symbols("n11 n12 n13 n22 n23 n33")
    n = sp.Matrix([[n11, n12, n13], [n12, n22, n23], [n13, n23, n33]])
    m_closed = closure_matrix(n)
    p = sp.Matrix([[0, 0, 1], [0, -1, 0], [1, 0, 0]])
    residual = sp.simplify(m_closed - p * n)
    residual_zero = residual == sp.zeros(3)
    det_p = sp.simplify(p.det())
    status = "PASS" if residual_zero and det_p != 0 else "FAIL"
    return {
        "description": "For pair order (12,13,23), M_closed = P*n for symmetric class-A n.",
        "pair_order": ["12", "13", "23"],
        "P_matrix": matrix_rows(p),
        "det_P": str(det_p),
        "residual_matrix": matrix_rows(residual),
        "residual_zero": bool(residual_zero),
        "status": status,
    }
