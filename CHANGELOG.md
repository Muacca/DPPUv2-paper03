# Changelog

## v3 — 2026-11-01

Changes are relative to v2 (2026-06-20).
The paper is retitled *Spin-0, Spin-1 and Spin-2 Deformations of S³×S¹ in
Einstein–Cartan–Nieh–Yan Theory*. The three spatial representation sectors
organize results for the Euclidean alpha=0 model and its stated variation families.

- Corrected the complete NY curvature coefficient from 1/4 to 1/2, the
  product volume to pi² L r³/4, and the circle convention to period one.
- Distinguished form-index Pontryagin P from internal-dual B. Standard P and
  complete NY vanish on the constant unimodular family, including MX;
  the former nonzero-MX Pontryagin interpretation is removed.
- Corrected inverse indices, mixing rotations, physical diagonal coframe
  factors and the off-diagonal shape transformation of twist.
- Retained spin-0 as a classification of scalar deformations and controls.
  Re-derived the common-coordinate Hessian and radius–torsion couplings;
  r=3 is a radius minimum at fixed eta=V=1, while the AX eta=1 radius
  potential is linear, not a flat zero mode.
- Derived the exact joint periodic circle action with the anisotropic shift,
  velocity cross terms and velocity-linear coupling. The scalar and vector
  circle actions follow as restrictions; torsion and coframe-rotation
  derivatives are retained as periodic boundary currents.
- Added the zero-twist shape derivative action and the full velocity rank.
  Quadratic sector separation at the round, zero-twist reference survives,
  while nonlinear derivative terms couple the sectors.
- Constructed a smooth coupled shape-twist family whose Euclidean action is
  unbounded below even at fixed scalar controls. Distinguished this functional
  from the static zero-twist shape potential and its round global minimum.
- Identified the twist source produced by a rotating anisotropic shape.
  The zero-twist stationary statement retains its fixed-scalar, round-shape scope.
- Derived the exact finite-amplitude twist circle action and its quartic
  direction term. The velocity matrix has rank two away from zero twist and
  vanishes at the origin. In the fixed-scalar, round-shape real Euclidean family, the
  only freely stationary twist is zero. The spin-1 triplet mass
  interpretations of v1 (stable) and v2 (tachyonic), with the
  kinetic-normalization bridge, are not retained.
- Re-derived the spin-2 curvatures with normalized exponential paths.
  The round coefficient is 24 pi² at r=3, L=kappa=1. Distinguished path
  curvature zeros from stationary transitions.
- Solved the full doublet minimization exactly: continuous bifurcation at
  epsilon=1/sqrt(2)-1, a square-root amplitude law, a C1 minimized potential
  with a second-derivative jump, and a residual rotation orbit of minima.
  Proved the round global minimum over all volume-preserving shapes at
  zero twist and fixed scalars; minimizing static twist jointly also selects zero.
- Included the cubic tensor invariant and explicit tensor–twist couplings.
  Replaced universal Landau and chirality-response claims with the stated
  geometric and variational results.
- Made the common static sign and coordinate dictionaries explicit relative to
  Papers01/02, including the normalization of the internal curvature contraction.
- Kept the curvature-squared coupling at zero. Ricci-tracefree RC and
  irreducible EC-Weyl are distinct projections; a nonzero-coupling atlas
  is outside this paper.
- Updated Japanese/English text, PDF, ten figure pairs, saved datasets,
  notebook and figure tools together. Rendering supplied data and explicit
  recomputation are separate operations with new output directories.
- The public code supplement consists of paper_results, the figure notebook
  and paper03_saved_data_viewer.py alongside the notebook. The general-purpose symbolic library and unrelated
  model utilities are no longer included; dependencies match the reproduction paths.
- Documented the numerical APIs and conventions in matching English and Japanese
  guides. The joint bulk function implements the stated spatially homogeneous
  periodic family; the coframe curvature calculation uses constant parameters.
- Streamlined the exposition around the geometric results and their variation
  domains, with code and data availability collected in Appendix H.

## v2 — 2026-06-20

The title was unchanged from v1. The spin-1 analysis eliminated the
auxiliary mixing combination before forming the generalized eigenvalue.
The reported spin-1 triplet eigenvalue changed from +4.1867 (locally stable)
to about -0.4183 (non-ghost, tachyonic), and a new appendix described the
reduction-order artifact. The nonzero MX Pontryagin interpretation and the
Landau EFT framing of v1 were retained; v3 withdraws them, as listed above.

## v1 — 2026-03-21

First release under the title *Unified Geometric Landau EFT of Homogeneous
S³×S¹ Minisuperspace in Einstein–Cartan + Nieh–Yan Theory*.

Figures and numerical conclusions of v1 and v2 use pre-v3 definitions
and should not be combined with the v3 datasets.
