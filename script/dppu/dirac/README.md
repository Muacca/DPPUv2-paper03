# Dirac Layer

-> [Japanese](README_ja.md)

Dirac constraint-closure engines used by the paper07 witnesses.  This package recomputes Branch-A local first-order Palatini models, scalar-axis projections, DBR closure, and Routh-before-Dirac route comparisons without relying on external analysis logs or intermediate CSV files.

## Overview

`dppu.dirac` is the computational core for the DBR/RBD route checks used in the paper.  The main callers are `scripts/paper07/rank_defect_activation_witness.py`, `scripts/paper07/harmonic_kernel_correspondence.py`, and `scripts/paper07/route_comparison_s3_branch_a_witness.py`.

## Modules

### `constraint_closure.py`

Builds Dirac constraint closure from the Branch-A local first-order Lagrangian and harmonic projections.

**Key functions:**

- `build_branchA_local_L2()`: construct the Branch-A local quadratic Lagrangian.
- `project_lagrangian(local, sample_l, l_sector, projection_variant="x")`: project a local model to a harmonic sector.
- `build_dirac_for_projection(proj, local)`: compute primary/secondary constraints, Poisson rank, and first/second-class split.
- `build_dirac_for_second_order_lagrangian(L_red, q_red, v_red)`: compute closure for a second-order reduced Lagrangian.
- `sample_geometry(sample_l, variant="x")`: representative harmonic geometry for route comparison.

### `local_models.py`

Class-A local-model builders for `T3`, `Nil3`, `Sol3`, and `S3`.  The rank-defect activation and harmonic-kernel witnesses build topology-specific local models from here.

**Key functions:**

- `build_branchA_local_L2_with_lc(topology)`: Branch-A local model with topology-specific LC structure.
- `build_t3_local()`: `T3` local model.
- `build_sol3_local()`: `Sol3` local model.
- `build_full_local_projection(local)`: full local projection helper.

### `scalar_axis.py`

Scalar-axis DBR closure engine.  It is used by witnesses checking `rank(A_DBR)=rank(n)` and the reduced-Hessian contact for active/inactive scalar axes.

**Key functions:**

- `build_scalar_axis_projection(local, axis)`: build an `x/y/z` scalar-axis projection.
- `dirac_standard(proj, local, rank_subs)`: run DBR standard closure.
- `reduced_hamiltonian_status(result)`: summarize the reduced Hamiltonian physical-pair status.

### `route_variants.py`

Negative-control helpers for Routh-before-Dirac routes.  The RBD variants are recomputed from the same Branch-A local model used by DBR.

**Key functions:**

- `routh_eliminate(local, also_K0)`: spatial-omega or full-K Routh elimination.
- `project_reduced_dirac(local, routh, sample_l)`: reduced Dirac closure after Routh elimination.

## Dependencies

- SymPy: symbolic Poisson brackets, constraint ranks, Hamiltonian/Hessian operations.
- [action](../action/README.md): Palatini/EH density helpers.
- [topology](../topology/README.md): class-A structure data.
- [spectral](../spectral/README.md): reduced-Hessian and harmonic-basis audits.

