"""
S3Geometry — S³×S¹ Topology Implementation
============================================

EC+NY engine on S³ (SU(2) Lie group) × S¹.

Supported DOFs:
    ε (squash)         : symmetry breaking SO(3) → U(1)
    s (shear)          : shear in the 0–1 plane (further breaks U(1))
    q₃,q₄,q₅ (off-diag shear): off-diagonal shear in the T₃,T₄,T₅ directions
    TWIST (ω₁,ω₂,ω₃)  : twist of the S¹ fiber
    MIXING (δ₀,δ₁,δ₂) : rotation mixing of ê^i and ê³ (i=0,1,2)

Spatial block of structure constants (i,j,k ∈ {0,1,2}):

    [Diagonal frame q₃=q₄=q₅=0]
    C^i_{jk} = (4/r) · factor_i · ε_{ijk}

    factor_0 = (1+ε)^{2/3} (1+s)²
    factor_1 = (1+ε)^{2/3} / (1+s)²
    factor_2 = (1+ε)^{−4/3}               (independent of s)

    [With off-diagonal shear]
    F = F_diag(ε,s) × G(q₃,q₄,q₅)   where G = exp(q₃T₃)exp(q₄T₄)exp(q₅T₅)
    C^i_{jk} = (4/r) · Σ F[i,l] F⁻¹[j,m] F⁻¹[k,n] ε_{lmn}

Additional fiber terms:
    TWIST        → C³_{jk}  (_add_s3_twist_C)
    MIXING/BOTH  → apply δ₀,δ₁,δ₂ rotation M=R₀₃×R₁₃×R₂₃ to all components

Volume factor: 2π² · L · r³
"""

from dataclasses import dataclass
from sympy import S, Rational, sqrt, cosh, sinh, Matrix, cancel
from sympy.tensor.array import MutableDenseNDimArray
from typing import Any, Dict, List, Tuple

from .base_topology import TopologyEngine
from ..utils.epsilon import epsilon_3d


class S3Geometry(TopologyEngine):
    """
    Engine dedicated to the S³×S¹ topology.

    Can be instantiated directly or via the ``UnifiedEngine`` factory.

    Examples
    --------
    Direct usage::

        from dppu.topology.s3 import S3Geometry
        from dppu.topology import make_config, TopologyType, FiberMode
        from dppu.torsion.mode import Mode
        from dppu.torsion.nieh_yan import NyVariant

        cfg = make_config(
            TopologyType.S3,
            enable_squash=True,
            enable_shear=True,
            fiber_mode=FiberMode.TWIST,
            isotropic_twist=False,
            torsion_mode=Mode.AX,
        )
        engine = S3Geometry(cfg)
        engine.run()
    """

    # ── E4.1 hook ─────────────────────────────────────────────────────────────

    def _build_radial_and_deformation_params(self, params: Dict) -> None:
        """
        Add S³-specific symbols to params.

            r       : radial coordinate (positive)
            epsilon : squash parameter (Symbol only when enable_squash=True)
            s       : shear parameter  (Symbol only when enable_shear=True)
            q3,q4,q5: physical parameters for off-diagonal shear (for reference; do not appear in Veff)
            z3,z4,z5: rationalized symbols z_i = exp(q_i/√2) (when enable_offdiag_shear=True)
                      These are the actual free variables of Veff.
                      Pass z_i = np.exp(q_i/√2) during numerical evaluation.
                      q₃: (0,1)-plane shear (T₃),  q₄: (0,2)-plane (T₄),  q₅: (1,2)-plane (T₅)
                      Frame: F = F_diag(ε,s) × G(z₃,z₄,z₅),  z=1 → G=I (matches legacy)

        Speedup: Method 2 (z = e^x rationalization)
            cosh(q/√2) = (z+1/z)/2,  sinh(q/√2) = (z-1/z)/2 (z = exp(q/√2))
            SymPy processes rational functions faster than transcendental functions,
            so the entire pipeline is computed as a rational function in z.
            No effect on accuracy (algebraic identity; numerical error < 4e-8).
        """
        from sympy import symbols, Symbol
        cfg = self.config

        params['r'] = symbols('r', positive=True, real=True)

        params['epsilon'] = (
            symbols('epsilon', real=True) if cfg.enable_squash else S.Zero
        )
        params['s'] = (
            symbols('s', real=True) if cfg.enable_shear else S.Zero
        )

        enable_offdiag = getattr(cfg, 'enable_offdiag_shear', False)
        # q3/q4/q5 are retained as physical parameters (for reference in lambdify wrappers etc.)
        params['q3'] = symbols('q3', real=True) if enable_offdiag else S.Zero
        params['q4'] = symbols('q4', real=True) if enable_offdiag else S.Zero
        params['q5'] = symbols('q5', real=True) if enable_offdiag else S.Zero
        # z3/z4/z5: Method 2 rationalized symbols (the actual free variables of Veff)
        if enable_offdiag:
            params['z3'] = Symbol('z3', positive=True)
            params['z4'] = Symbol('z4', positive=True)
            params['z5'] = Symbol('z5', positive=True)
        else:
            params['z3'] = params['z4'] = params['z5'] = None

    # ── E4.2 hook ─────────────────────────────────────────────────────────────

    def _build_structure_constants(
        self, params: Dict, C: MutableDenseNDimArray
    ) -> None:
        """
        Write the S³ spatial block + fiber terms into C in-place.

        Spatial block (i,j,k ∈ {0,1,2}):
          [Diagonal frame only (q₃=q₄=q₅=0)]
            C^i_{jk} = (4/r) · f_i/(f_j·f_k) · ε_{ijk}
            f_i/(f_j·f_k) = spatial_factors[i]

          [With off-diagonal shear (q₃,q₄,q₅ ≠ 0)]
            F = F_diag(ε,s) × G(q₃,q₄,q₅)
            C^i_{jk} = (4/r) · Σ_{l,m,n} F[i,l] F⁻¹[j,m] F⁻¹[k,n] ε_{lmn}

        δ normalization for MIXING / BOTH:
          n_delta is not included in the spatial block. The rotation M(δ₀,δ₁,δ₂)
          automatically introduces cos/sin factors (= 1/n, δ/n), making them equivalent.
        """
        from .unified import FiberMode

        cfg     = self.config
        r       = params['r']
        L       = params['L']
        epsilon = params['epsilon']
        s       = params['s']
        q3      = params['q3']
        q4      = params['q4']
        q5      = params['q5']

        # ── Spatial frame diagonal factors ──────────────────────────────────────────
        #   spatial_factors[i] = f_i / (f_j · f_k)
        factor_0 = (1 + epsilon) ** Rational(2, 3)  * (1 + s) ** 2
        factor_1 = (1 + epsilon) ** Rational(2, 3)  / (1 + s) ** 2
        factor_2 = (1 + epsilon) ** Rational(-4, 3)      # independent of s
        spatial_factors = [factor_0, factor_1, factor_2]

        # ── Spatial block ─────────────────────────────────────────────────
        if q3 == S.Zero and q4 == S.Zero and q5 == S.Zero:
            # Fast path: diagonal frame (matches existing behavior; agrees exactly at q₃=q₄=q₅=0)
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        eps_val = epsilon_3d(i, j, k)
                        if eps_val != 0:
                            C[i, j, k] = eps_val * 4 * spatial_factors[i] / r
        else:
            # General path: F = F_diag × G(z₃,z₄,z₅)  [Method 2: z_i = exp(q_i/√2) rationalization]
            # G is a product of symmetric positive-definite matrices (composition of hyperbolic rotations)
            #
            # Speedup principle:
            #   cosh(q/√2) = (z + 1/z)/2,  sinh(q/√2) = (z - 1/z)/2  (z = exp(q/√2))
            #   SymPy processes rational functions faster than transcendental functions (cosh/sinh).
            #   cancel() folds the entire pipeline into rational polynomials,
            #   completing in ~10s for all 3 variables simultaneously (previously: >30 min).
            #   No effect on accuracy (algebraic identity).
            #
            # Uses params['z3'], params['z4'], params['z5'].
            # z=1 (q=0) → G=I (matches the diagonal frame).
            z3 = params['z3']
            z4 = params['z4']
            z5 = params['z5']
            c3, s3_ = (z3 + 1/z3) / 2, (z3 - 1/z3) / 2
            c4, s4_ = (z4 + 1/z4) / 2, (z4 - 1/z4) / 2
            c5, s5_ = (z5 + 1/z5) / 2, (z5 - 1/z5) / 2

            G3 = Matrix([[c3,  s3_, S.Zero],
                         [s3_, c3,  S.Zero],
                         [S.Zero, S.Zero, S.One]])
            G4 = Matrix([[c4,  S.Zero, s4_],
                         [S.Zero, S.One, S.Zero],
                         [s4_, S.Zero, c4]])
            G5 = Matrix([[S.One, S.Zero, S.Zero],
                         [S.Zero, c5,  s5_],
                         [S.Zero, s5_, c5]])

            G = G3 * G4 * G5
            F_diag = Matrix.diag(factor_0, factor_1, factor_2)
            F      = F_diag * G

            # Explicit inverse: F_inv = G5_inv × G4_inv × G3_inv × F_diag_inv
            # Each Gi is a hyperbolic rotation: Gi_inv simply negates the sign of the sinh term.
            G3_inv = Matrix([[c3,  -s3_, S.Zero],
                             [-s3_, c3,  S.Zero],
                             [S.Zero, S.Zero, S.One]])
            G4_inv = Matrix([[c4,  S.Zero, -s4_],
                             [S.Zero, S.One, S.Zero],
                             [-s4_, S.Zero, c4]])
            G5_inv = Matrix([[S.One, S.Zero, S.Zero],
                             [S.Zero, c5,  -s5_],
                             [S.Zero, -s5_, c5]])
            F_diag_inv = Matrix.diag(
                S.One / factor_0, S.One / factor_1, S.One / factor_2)
            F_inv = G5_inv * G4_inv * G3_inv * F_diag_inv

            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        v = S.Zero
                        for l in range(3):
                            for m in range(3):
                                for n in range(3):
                                    eps_val = epsilon_3d(l, m, n)
                                    if eps_val != 0:
                                        v += (F[i, l] * F_inv[j, m]
                                              * F_inv[k, n] * eps_val)
                        if v != S.Zero:
                            # Simplify as a rational function via cancel() (the core of the speedup)
                            C[i, j, k] = cancel(4 * v / r)

        # ── Fiber terms ─────────────────────────────────────────────────────
        if cfg.fiber_mode == FiberMode.TWIST:
            self._add_s3_twist_C(C, params, epsilon, s, r, L)
        elif cfg.fiber_mode in (FiberMode.MIXING, FiberMode.BOTH):
            # BOTH: add twist terms first, then apply the δ rotation
            # MIXING: no twist terms → apply δ rotation to spatial block only
            # In both cases, _apply_s3_mixing_rotation_to_C automatically generates O(δω) cross terms
            if cfg.fiber_mode == FiberMode.BOTH:
                self._add_s3_twist_C(C, params, epsilon, s, r, L)
            self._apply_s3_mixing_rotation_to_C(C, params)

        # ── Velocity terms (nonzero only when enable_velocity=True) ───────────────────
        self._compute_velocity_terms(params, C)

    def _compute_velocity_terms(
        self, params: Dict, C: MutableDenseNDimArray
    ) -> None:
        """
        Add velocity terms to the structure constants (for G metric computation; active only when enable_velocity=True).

        First-order correction to the structure constants from the τ-direction velocities
        v^ω_k = ∂_τω_k, v^δ_k = ∂_τδ_k at the isotropic point (ω=δ=ε=s=0):

            C^3_{3k} += +2 v^ω_k / r  −  v^δ_k / L    (k = 0,1,2)
            C^3_{k3}  = −C^3_{3k}                       (antisymmetric)

        Index correspondence (0-indexed k):
            k=0 → v_omega1 (ω₁ in the σ⁰ direction), v_delta0
            k=1 → v_omega2 (ω₂ in the σ¹ direction), v_delta1
            k=2 → v_omega3 (ω₃ in the σ² direction), v_delta2

        Reference: 00_plan/20260307-05_G_metric_calculation_plan.md §1.2
        """
        r = params['r']
        L = params['L']
        v_w = [params.get('v_omega1', S.Zero),
               params.get('v_omega2', S.Zero),
               params.get('v_omega3', S.Zero)]
        v_d = [params.get('v_delta0', S.Zero),
               params.get('v_delta1', S.Zero),
               params.get('v_delta2', S.Zero)]
        for k in range(3):
            vel = 2 * v_w[k] / r - v_d[k] / L
            C[3, 3, k] = C[3, 3, k] + vel
            C[3, k, 3] = C[3, k, 3] - vel

    def _compute_volume(self, params: Dict) -> Any:
        """2π² · L · r³ (volume of S³ × S¹)"""
        from sympy import pi
        return 2 * pi ** 2 * params['L'] * params['r'] ** 3

    # ── Off-diagonal shear numerical conversion helper ──────────────────────────────

    @staticmethod
    def offdiag_zi_from_qi(q3=0.0, q4=0.0, q5=0.0):
        """
        Return the numerical values of rationalized symbols z_i from physical parameters q_i.

        This conversion is required during numerical evaluation because Veff is expressed in z3,z4,z5.

        Parameters
        ----------
        q3, q4, q5 : float
            Physical parameters of off-diagonal shear (q=0 → G=I)

        Returns
        -------
        z3, z4, z5 : float
            z_i = exp(q_i / √2)

        Examples
        --------
        >>> z3, z4, z5 = S3Geometry.offdiag_zi_from_qi(q4=0.1)
        >>> # Veff_func(r, eps, z3, z4, z5) = Veff(r, eps, q4=0.1)
        """
        import math
        sq2 = math.sqrt(2.0)
        return math.exp(q3/sq2), math.exp(q4/sq2), math.exp(q5/sq2)


# ── Module-level helpers ──────────────────────────────────────────────────────

def s3_lie_frame_structure_constants(
    q: Any,
    n_lapse: Any,
    qdot: Any,
    lam: Any,
    dim: int = 4,
) -> "MutableDenseNDimArray":
    """Build S3×time Lie-frame structure constants for a cosmological ansatz.

    Returns C^a_{bc} (4D) for the coframe ``e^0 = N dt``, ``e^i = q σ^i``
    (i = 1,2,3) with explicit lapse *n_lapse*, scale factor *q*, its time
    derivative *qdot*, and spatial curvature normalisation *lam*.

    Non-zero components:

    * Time-spatial block (Hubble terms):
        C^i_{0i} = h,  C^i_{i0} = −h   where  h = qdot / (n_lapse · q)

    * Spatial block (Lie-algebra of SU(2)):
        C^i_{jk} = lam · ε_{i-1, j-1, k-1} / q   for  i, j, k ∈ {1,2,3}

    Args:
        q:        Scale-factor symbol (positive).
        n_lapse:  Lapse function symbol N.
        qdot:     Time derivative of q (symbol or expression).
        lam:      SU(2) Lie-algebra normalisation factor.
                  ``lam = 2`` gives the unit-round S3 conventions used in
                  Phase5A; ``lam = 4/r`` is the LZ-native convention.
        dim:      Space-time dimension (default 4).

    Returns:
        ``MutableDenseNDimArray`` of shape ``(dim, dim, dim)``.
    """
    from sympy.tensor.array import MutableDenseNDimArray

    c_arr = MutableDenseNDimArray.zeros(dim, dim, dim)
    h = qdot / (n_lapse * q)

    for i in range(1, dim):
        c_arr[i, 0, i] = h
        c_arr[i, i, 0] = -h

    for i in range(1, dim):
        for j in range(1, dim):
            for k in range(1, dim):
                eps_val = epsilon_3d(i - 1, j - 1, k - 1)
                if eps_val:
                    c_arr[i, j, k] = lam * eps_val / q
    return c_arr


# ── S3 SVT linearised decomposition symbols ───────────────────────────────────

@dataclass
class S3SVTSymbols:
    """S3 background coframe + linearised branch-A SVT decomposition symbols.

    Encodes the coframe perturbation δe^a_j around the S3 background
    (N dt, q σ^i) into SVT sectors:

    * scalar    : ``s``  (trace part, coefficient 1/3)
    * shift     : ``N_1, N_2, N_3``  (×3)
    * spin-2/qA : ``qA_00, qA_11, qA_01, qA_02, qA_12``  traceless symmetric  (×5)
    * spin-1/X  : ``X_01, X_02, X_12``  antisymmetric  (×3)

    Also provides component dicts ``qA_components``, ``X_components``,
    the linearised coframe perturbation dict ``de[(a, j)]``,
    and the ordered coordinate list ``coords`` (12 symbols, canonical order).
    """

    N: Any                             # lapse symbol
    q: Any                             # spatial scale-factor symbol
    s: Any                             # scalar perturbation
    N_shift: List                      # [N_1, N_2, N_3]
    qA00: Any
    qA11: Any
    qA01: Any
    qA02: Any
    qA12: Any
    qA22: Any                          # derived: -(qA00 + qA11)
    X01: Any
    X02: Any
    X12: Any
    qA_components: Dict                # {(i,j): expr} for i,j in {0,1,2}
    X_components:  Dict                # {(i,j): expr}
    de: Dict                           # {(a, j): expr} for a in {0..3}, j in {1..3}
    coords: List                       # 12 SVT symbols in canonical order


def s3_svt_symbols(
    N: Any = None,
    q: Any = None,
) -> S3SVTSymbols:
    """Build S3 background + linearised SVT decomposition symbols.

    Parameters
    ----------
    N : sympy.Symbol, optional
        Lapse symbol.  Defaults to ``sp.Symbol('N')``.
    q : sympy.Symbol, optional
        Spatial scale-factor symbol.  Defaults to ``sp.Symbol('q')``.

    Returns
    -------
    S3SVTSymbols
        Dataclass with all symbols and component dicts for stage2/stage3
        S3 linearised Dirac computations.
    """
    import sympy as sp

    if N is None:
        N = sp.Symbol("N")
    if q is None:
        q = sp.Symbol("q")

    s: Any = sp.Symbol("s")
    N_shift: List = [sp.Symbol(f"N_{i}") for i in range(1, 4)]

    qA00, qA11, qA01, qA02, qA12 = sp.symbols("qA_00 qA_11 qA_01 qA_02 qA_12")
    qA22 = -(qA00 + qA11)
    qA_components: Dict[Tuple[int, int], Any] = {
        (0, 0): qA00, (1, 1): qA11, (2, 2): qA22,
        (0, 1): qA01, (1, 0): qA01,
        (0, 2): qA02, (2, 0): qA02,
        (1, 2): qA12, (2, 1): qA12,
    }

    X01, X02, X12 = sp.symbols("X_01 X_02 X_12")
    X_components: Dict[Tuple[int, int], Any] = {
        (0, 1): X01,  (1, 0): -X01,
        (0, 2): X02,  (2, 0): -X02,
        (1, 2): X12,  (2, 1): -X12,
        (0, 0): sp.S.Zero, (1, 1): sp.S.Zero, (2, 2): sp.S.Zero,
    }

    de: Dict[Tuple[int, int], Any] = {}
    for a in range(4):
        for j in range(1, 4):
            if a == 0:
                de[(a, j)] = N_shift[j - 1]
            else:
                ia, ij = a - 1, j - 1
                val: Any = qA_components[(ia, ij)] + X_components[(ia, ij)]
                if ia == ij:
                    val = val + s / 3
                de[(a, j)] = val

    coords: List = [s] + N_shift + [qA00, qA11, qA01, qA02, qA12, X01, X02, X12]

    return S3SVTSymbols(
        N=N, q=q, s=s,
        N_shift=N_shift,
        qA00=qA00, qA11=qA11, qA01=qA01, qA02=qA02, qA12=qA12,
        qA22=qA22,
        X01=X01, X02=X02, X12=X12,
        qA_components=qA_components,
        X_components=X_components,
        de=de,
        coords=coords,
    )
