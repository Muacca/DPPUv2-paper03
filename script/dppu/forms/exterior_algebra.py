"""
Exterior Differential Form Algebra
===================================

Pure exterior algebra operations on differential forms.

Forms are represented as dictionaries mapping ordered index tuples
to SymPy coefficients::

    Form = Dict[Tuple[int, ...], SymPy_expr]

The empty tuple ``()`` denotes a scalar 0-form.
Index tuples are always kept in ascending sorted order; the
appropriate sign is absorbed into the coefficient when necessary.

Example::

    from dppu.forms.exterior_algebra import basis, wedge, form_to_str
    e0 = basis(0)   # {(0,): 1}
    e1 = basis(1)   # {(1,): 1}
    e01 = wedge(e0, e1)  # {(0, 1): 1}

Author: Muacca
"""

from __future__ import annotations

from collections import defaultdict
from typing import Callable, Dict, Iterable, Optional, Tuple

from sympy import S, factor, simplify

# ---------------------------------------------------------------------------
# Type alias
# ---------------------------------------------------------------------------

Form = Dict[Tuple[int, ...], object]


# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------

def clean_form(form: Form, normalize: Optional[Callable] = None) -> Form:
    """Remove zero-coefficient terms from *form*.

    Args:
        form:      Input form dict.
        normalize: Optional callable applied to each coefficient before
                   zero-testing.  Defaults to ``sympy.simplify`` when
                   not provided (backward-compatible).
    """
    _norm = normalize if normalize is not None else simplify
    out: Form = {}
    for key, val in form.items():
        sval = _norm(val)
        if sval != S.Zero:
            out[key] = sval
    return out


# ---------------------------------------------------------------------------
# Constructors
# ---------------------------------------------------------------------------

def basis(i: int) -> Form:
    """Return the elementary 1-form ``e^i``."""
    return {(i,): S.One}


def scalar_form(value, normalize: Optional[Callable] = None) -> Form:
    """Wrap *value* as a 0-form.  Returns ``{}`` if *value* is zero.

    Args:
        value:     Scalar coefficient.
        normalize: Optional coefficient normalizer.  Defaults to
                   ``sympy.simplify``.
    """
    _norm = normalize if normalize is not None else simplify
    value = _norm(value)
    return {} if value == S.Zero else {(): value}


# ---------------------------------------------------------------------------
# Arithmetic
# ---------------------------------------------------------------------------

def add_forms(*forms: Form, normalize: Optional[Callable] = None) -> Form:
    """Return the sum of any number of forms.

    Args:
        forms:     Forms to add.
        normalize: Optional coefficient normalizer passed to
                   :func:`clean_form`.  Defaults to ``sympy.simplify``.
    """
    out: Form = {}
    for form in forms:
        for key, val in form.items():
            out[key] = out.get(key, S.Zero) + val
    return clean_form(out, normalize=normalize)


def scale_form(value, form: Form, normalize: Optional[Callable] = None) -> Form:
    """Return *value* · *form*.

    Args:
        value:     Scalar multiplier.
        form:      Form to scale.
        normalize: Optional coefficient normalizer passed to
                   :func:`clean_form`.  Defaults to ``sympy.simplify``.
    """
    return clean_form(
        {key: value * val for key, val in form.items()},
        normalize=normalize,
    )


def wedge(
    left: Form,
    right: Form,
    normalize: Optional[Callable] = None,
) -> Form:
    """Return the exterior product ``left ∧ right``.

    Args:
        left:      Left factor form.
        right:     Right factor form.
        normalize: Optional coefficient normalizer passed to
                   :func:`clean_form`.  Defaults to ``sympy.simplify``.
    """
    out: Form = {}
    for left_idx, left_val in left.items():
        for right_idx, right_val in right.items():
            if set(left_idx).intersection(right_idx):
                continue  # Repeated index → zero
            seq = list(left_idx) + list(right_idx)
            inversions = sum(
                1
                for i in range(len(seq))
                for j in range(i + 1, len(seq))
                if seq[i] > seq[j]
            )
            sign = -S.One if inversions % 2 else S.One
            key = tuple(sorted(seq))
            out[key] = out.get(key, S.Zero) + sign * left_val * right_val
    return clean_form(out, normalize=normalize)


# ---------------------------------------------------------------------------
# Display and extraction
# ---------------------------------------------------------------------------

def form_to_str(form: Form) -> str:
    """Human-readable representation of *form*."""
    if not form:
        return "0"
    terms = []
    for key in sorted(form):
        basis_name = "1" if key == () else "e" + "".join(str(i) for i in key)
        terms.append(f"{factor(form[key])}*{basis_name}")
    return " + ".join(terms)


def coefficient(form: Form, key: Tuple[int, ...] = (0, 1, 2, 3)) -> object:
    """Extract and simplify the coefficient for *key* in *form*."""
    return simplify(form.get(key, S.Zero))


# ---------------------------------------------------------------------------
# Maurer–Cartan exterior derivative for class-A left-invariant 1-forms
# ---------------------------------------------------------------------------

def d_basis(
    index: int,
    structure_entries: Dict[Tuple[int, int, int], object],
    normalize: Optional[Callable] = None,
) -> Form:
    """Exterior derivative of the left-invariant basis 1-form ``σ^index``.

    Uses the Maurer–Cartan equation::

        dσ^up = -(1/2) C^up_{ij} σ^i ∧ σ^j

    where *structure_entries* is a dict ``{(up, i, j): coeff}`` in the
    paper class-A convention (1-based spatial indices; index 0 is the
    time direction and has a trivially zero derivative).

    Args:
        index:             Which basis 1-form to differentiate (1-based).
                           Index 0 returns an empty form (time component).
        structure_entries: Dict of nonzero structure constants
                           ``{(up, i, j): coeff}``.
        normalize:         Optional coefficient normalizer.  When ``None``
                           the coefficients are normalized with
                           ``sympy.expand`` (matching the paper07 class-A
                           implementation exactly).

    Returns:
        Form representing ``dσ^index``.

    Notes:
        - Entries with ``up != index`` are ignored.
        - If no entries exist for *index* (e.g. T3, or index == 0),
          the returned form is ``{}``, consistent with the class-A
          implementation where the index==0 special case is recovered
          naturally via absent entries.
    """
    import sympy as sp
    _expand = normalize if normalize is not None else sp.expand
    if index == 0:
        return {}
    out: Dict[Tuple[int, ...], object] = defaultdict(lambda: S.Zero)
    for (up, i, j), coeff in structure_entries.items():
        if up != index:
            continue
        merged = [i, j]
        if len(set(merged)) != len(merged):
            continue
        inversions = sum(
            1 for a in range(len(merged)) for b in range(a + 1, len(merged))
            if merged[a] > merged[b]
        )
        sign = (-1 if inversions % 2 else 1)
        key = tuple(sorted(merged))
        out[key] += -S.Rational(1, 2) * coeff * sign
    return {k: _expand(v) for k, v in out.items() if v != 0}


def d_form(
    form: Form,
    structure_entries: Dict[Tuple[int, int, int], object],
    normalize: Optional[Callable] = None,
) -> Form:
    """Exterior derivative of a left-invariant 1-form on a class-A Lie group.

    Applies :func:`d_basis` linearly to each basis component of *form*,
    then sums the results.  Only 1-form components (single-index keys) are
    differentiated; higher-degree components are ignored.

    Args:
        form:              Input 1-form (only ``len(key)==1`` entries are
                           differentiated).
        structure_entries: Structure constant dict ``{(up, i, j): coeff}``.
        normalize:         Optional coefficient normalizer; forwarded to
                           :func:`d_basis` and :func:`add_forms`.

    Returns:
        2-form ``d(form)``.
    """
    import sympy as sp
    _expand = normalize if normalize is not None else sp.expand
    result: Form = {}
    for key, coeff in form.items():
        if len(key) != 1:
            continue
        db = d_basis(key[0], structure_entries, normalize=normalize)
        # scale each coefficient
        scaled: Form = {k: _expand(coeff * v) for k, v in db.items() if v != 0}
        result = add_forms(result, scaled, normalize=normalize)
    return result
