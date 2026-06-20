"""
Spin-1 Gradient and Shift Derivation Helpers
============================================

Utilities for the Phase5B spin-1 X/Y coefficient derivation.

The homogeneous raw coefficients are obtained by running the existing S3
topology, Levi-Civita connection, Riemann, and Ricci-scalar primitives.  The
full-gradient coefficients are then built as the covariant principal-symbol
completion of the same raw rank-one Hessian in the 5B-2 vector-gradient basis.

This module intentionally performs no report/file I/O.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp
from sympy import Matrix, S

from ..curvature.ricci import compute_ricci_scalar
from ..curvature.riemann import compute_riemann_tensor
from ..engine.levi_civita import compute_christoffel_frame
from ..topology.s3 import S3Geometry
from ..topology.unified import FiberMode, TopologyType, make_config
from ..torsion.mode import Mode
from ..torsion.nieh_yan import NyVariant


DIM = 4
LORENTZ_METRIC = Matrix.diag(-1, 1, 1, 1)


@dataclass(frozen=True)
class RawSpin1Hessians:
    """Raw omega/delta Hessians extracted from the S3 Ricci scalar."""

    r: sp.Symbol
    L: sp.Symbol
    H_velocity: Matrix
    H_static: Matrix
    R_velocity: sp.Expr
    R_static: sp.Expr


@dataclass(frozen=True)
class XYSpin1Coefficients:
    """X/Y coefficients in the Ricci-scalar normalization."""

    r: sp.Symbol
    L: sp.Symbol
    A_X: sp.Expr
    C_X: sp.Expr
    C_XY: sp.Expr
    C_Y: sp.Expr
    Y_factor: sp.Expr
    reduced_static_X: sp.Expr
    H_velocity_XY: Matrix
    H_static_XY: Matrix


def _make_s3_data(*, enable_velocity: bool) -> tuple[dict[str, Any], Any]:
    cfg = make_config(
        TopologyType.S3,
        enable_squash=False,
        enable_shear=False,
        enable_offdiag_shear=False,
        fiber_mode=FiberMode.BOTH,
        isotropic_twist=False,
        torsion_mode=Mode.AX,
        ny_variant=NyVariant.FULL,
        enable_velocity=enable_velocity,
        signature="lorentzian",
    )
    engine = S3Geometry(cfg)
    engine.step_E4_1_setup()
    engine.step_E4_2_metric_and_frame()
    return engine.data["params"], engine.data["structure_constants"]


def _ricci_scalar_from_structure_constants(c_arr: Any) -> sp.Expr:
    gamma_lc = compute_christoffel_frame(
        c_arr,
        DIM,
        metric=LORENTZ_METRIC,
        metric_inv=LORENTZ_METRIC,
    )
    riemann = compute_riemann_tensor(gamma_lc, c_arr, DIM)
    return sp.cancel(compute_ricci_scalar(riemann, DIM, metric_inv=LORENTZ_METRIC))


def _zero_spin1_substitutions(params: dict[str, Any]) -> dict[Any, Any]:
    result: dict[Any, Any] = {}
    for key, value in params.items():
        if key.startswith("omega") or key.startswith("delta") or key.startswith("sd"):
            if hasattr(value, "is_Symbol") and value.is_Symbol:
                result[value] = S.Zero
        if key.startswith("cd") and value is not None:
            if hasattr(value, "is_Symbol") and value.is_Symbol:
                result[value] = S.One
    return result


def derive_raw_spin1_hessians(axis: int = 0) -> RawSpin1Hessians:
    """
    Derive raw omega/delta velocity and static Hessians from the S3 Ricci scalar.

    ``axis=0`` corresponds to ``omega1`` and ``delta0``.  The isotropic S3
    background makes the three axes equivalent; callers can verify this by
    repeating the extraction for axes 1 and 2 if desired.
    """
    if axis not in {0, 1, 2}:
        raise ValueError("axis must be 0, 1, or 2")

    params_v, c_vel = _make_s3_data(enable_velocity=True)
    r = params_v["r"]
    l_scale = params_v["L"]
    r_velocity = _ricci_scalar_from_structure_constants(c_vel)
    v_omega = params_v[f"v_omega{axis + 1}"]
    v_delta = params_v[f"v_delta{axis}"]
    h_velocity = Matrix(
        [
            [sp.diff(r_velocity, v_omega, v_omega), sp.diff(sp.diff(r_velocity, v_omega), v_delta)],
            [sp.diff(sp.diff(r_velocity, v_delta), v_omega), sp.diff(r_velocity, v_delta, v_delta)],
        ]
    )
    h_velocity = h_velocity.applyfunc(lambda item: sp.factor(sp.cancel(item)))

    params_s, c_static = _make_s3_data(enable_velocity=False)
    r_static = _ricci_scalar_from_structure_constants(c_static)
    omega = params_s[f"omega{axis + 1}"]
    sd = params_s[f"sd{axis}"]
    zero = _zero_spin1_substitutions(params_s)
    h_static = Matrix(
        [
            [
                sp.diff(r_static, omega, omega).subs(zero),
                sp.diff(sp.diff(r_static, omega), sd).subs(zero),
            ],
            [
                sp.diff(sp.diff(r_static, sd), omega).subs(zero),
                sp.diff(r_static, sd, sd).subs(zero),
            ],
        ]
    )
    h_static = h_static.applyfunc(lambda item: sp.factor(sp.cancel(item)))

    return RawSpin1Hessians(
        r=r,
        L=l_scale,
        H_velocity=h_velocity,
        H_static=h_static,
        R_velocity=r_velocity,
        R_static=r_static,
    )


def xy_transform_matrix(r: sp.Symbol, l_scale: sp.Symbol) -> Matrix:
    """
    Return the raw-from-XY transform for the Phase5B X/Y basis.

    The basis is the one fixed by Phase5A/5B-2:
      X = (-a omega + delta)/sqrt(1+a^2), Y = (omega + a delta)/sqrt(1+a^2),
      a = 2L/r.
    The returned matrix maps [X, Y]^T to [omega, delta]^T.
    """
    a_ratio = 2 * l_scale / r
    norm = sp.sqrt(1 + a_ratio**2)
    return Matrix([[-a_ratio / norm, 1 / norm], [1 / norm, a_ratio / norm]])


def derive_xy_coefficients(raw: RawSpin1Hessians) -> XYSpin1Coefficients:
    """Rotate raw Hessians into the X/Y basis and extract scalar coefficients."""
    transform = xy_transform_matrix(raw.r, raw.L)
    h_velocity_xy = (transform.T * raw.H_velocity * transform).applyfunc(
        lambda item: sp.factor(sp.cancel(item))
    )
    h_static_xy = (transform.T * raw.H_static * transform).applyfunc(
        lambda item: sp.factor(sp.cancel(item))
    )

    a_x = sp.factor(sp.cancel(h_velocity_xy[0, 0] / 2))
    c_x = sp.factor(sp.cancel(h_static_xy[0, 0] / 2))
    c_xy = sp.factor(sp.cancel(h_static_xy[0, 1]))
    c_y = sp.factor(sp.cancel(h_static_xy[1, 1] / 2))
    y_factor = sp.factor(sp.cancel(-c_xy / (2 * c_y)))
    reduced_static_x = sp.factor(sp.cancel(c_x - c_xy**2 / (4 * c_y)))

    return XYSpin1Coefficients(
        r=raw.r,
        L=raw.L,
        A_X=a_x,
        C_X=c_x,
        C_XY=c_xy,
        C_Y=c_y,
        Y_factor=y_factor,
        reduced_static_X=reduced_static_x,
        H_velocity_XY=h_velocity_xy,
        H_static_XY=h_static_xy,
    )


def spin1_operator_coefficients(coeffs: XYSpin1Coefficients) -> list[dict[str, Any]]:
    """
    Build the 5B-2 vector-gradient operator coefficients.

    The principal symbol is the covariant completion
      A_X D_i X_j D^i X^j
    of the rank-one X velocity block extracted from the dppu Ricci scalar.
    Using the 5B-2 decomposition gives:
      D_i X_j D^i X^j = S_X:S_X + theta_X^2/3 + curl_X^2/2.
    Y and X/Y mixed principal-gradient coefficients are zero in this branch.
    """
    a_x = coeffs.A_X
    zero = S.Zero
    return [
        {
            "operator_id": "GRAD-XX-DIV",
            "topology": "S3",
            "operator": "theta_X^2",
            "ricci_coefficient": sp.factor(sp.cancel(a_x / 3)),
            "curvature_commutator_part": zero,
        },
        {
            "operator_id": "GRAD-XX-CURL",
            "topology": "S3",
            "operator": "(curl X)^2",
            "ricci_coefficient": sp.factor(sp.cancel(a_x / 2)),
            "curvature_commutator_part": zero,
        },
        {
            "operator_id": "GRAD-XX-SHEAR",
            "topology": "S3",
            "operator": "S_X:S_X",
            "ricci_coefficient": a_x,
            "curvature_commutator_part": zero,
        },
        {
            "operator_id": "GRAD-YY-DIV",
            "topology": "S3",
            "operator": "theta_Y^2",
            "ricci_coefficient": zero,
            "curvature_commutator_part": zero,
        },
        {
            "operator_id": "GRAD-YY-CURL",
            "topology": "S3",
            "operator": "(curl Y)^2",
            "ricci_coefficient": zero,
            "curvature_commutator_part": zero,
        },
        {
            "operator_id": "GRAD-YY-SHEAR",
            "topology": "S3",
            "operator": "S_Y:S_Y",
            "ricci_coefficient": zero,
            "curvature_commutator_part": zero,
        },
        {
            "operator_id": "GRAD-XY-DIV",
            "topology": "S3",
            "operator": "theta_X theta_Y",
            "ricci_coefficient": zero,
            "curvature_commutator_part": zero,
        },
        {
            "operator_id": "GRAD-XY-CURL",
            "topology": "S3",
            "operator": "(curl X).(curl Y)",
            "ricci_coefficient": zero,
            "curvature_commutator_part": zero,
        },
        {
            "operator_id": "GRAD-XY-SHEAR",
            "topology": "S3",
            "operator": "S_X:S_Y",
            "ricci_coefficient": zero,
            "curvature_commutator_part": zero,
        },
    ]


def spin1_shift_source_coefficients(coeffs: XYSpin1Coefficients, n_lapse: Any | None = None) -> dict[str, Any]:
    """
    Derive the spin-1 shift velocity block from the ADM covariant square.

    The spin-1 transverse shift projection is normalized to the same component
    coordinate as X_i.  The velocity block is constructed as
      a (Xdot_i - N_i)^2,
    where ``a = N r^3 A_X`` is the density coefficient from the Ricci scalar.
    """
    n = sp.symbols("N", positive=True, real=True) if n_lapse is None else n_lapse
    a_density = sp.factor(sp.cancel(n * coeffs.r**3 * coeffs.A_X))
    b_density = sp.factor(sp.cancel(-a_density))
    c_density = a_density
    schur = sp.factor(sp.cancel(a_density - b_density**2 / c_density))
    return {
        "a_density": a_density,
        "b_density": b_density,
        "c_density": c_density,
        "schur_density": schur,
        "pi_contribution": "pi_X_i = 2*a_density*(Xdot_i - N_i)",
        "momentum_source": "H_i receives -2*a_density*(Xdot_i - N_i) in the spin-1 projection",
    }


def curvature_commutator_audit(coeffs: XYSpin1Coefficients) -> dict[str, Any]:
    """Evaluate the round-S3 commutator channel and its T3 flat limit."""
    kappa = sp.symbols("kappa", real=True)
    x_sq = sp.symbols("Xsq", real=True)
    commutator_density = sp.factor(sp.cancel(2 * kappa * coeffs.A_X * x_sq))
    return {
        "s3_kappa": S.One,
        "t3_kappa": S.Zero,
        "representative_density": commutator_density,
        "s3_value": commutator_density.subs(kappa, S.One),
        "t3_value": commutator_density.subs(kappa, S.Zero),
    }
