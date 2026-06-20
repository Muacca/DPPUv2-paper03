# Spectral Layer

-> [Japanese](README_ja.md)

Spectral computation engines used by the paper07 witnesses: harmonic bases, constraint-surface reduced Hessian audits, and the `Sol3` compact-edge spectrum.

## Overview

`dppu.spectral` provides the closed 1-form / `ker n` correspondence, reduced-Hessian kernel contact after DBR closure, and the `Sol3` compact-edge A-dual-orbit modified-Mathieu positivity witness.  The main callers are `scripts/paper07/harmonic_kernel_correspondence.py` and `scripts/paper07/spectral_completion_witness.py`.

## Modules

### `harmonic_basis.py`

Builds closure matrices and harmonic/co-closed audits for left-invariant 1-forms from class-A structure constants.

**Key functions:**

- `structure_constants(n)`: build structure constants from a class-A structure matrix `n`.
- `closure_matrix(n)`: return the closed 1-form condition matrix.
- `basis_matrix(nullspace)`: format a nullspace basis as a matrix.
- `span_equal(left, right)`: check equality of two spans.
- `co_closed_audit(constants, k_form_basis)`: audit co-closed/harmonic status.
- `audit_case(topology, bianchi, n, case_id)`: topology-specific harmonic-basis audit.
- `generic_identity_check()`: notation check for the formal structure-map identity.

### `reduced_hessian.py`

Audits the Hessian kernel on the constraint-surface reduced phase space after DBR closure.

**Key functions:**

- `reduced_hessian_pair_audit(result)`: check the relation between reduced symplectic pairs and the Hamiltonian Hessian kernel.
- `scalar_axis_hessian_audits(topology, local, active_axes)`: run scalar-axis DBR Hessian audits.

### `sol3_compact_edge.py`

Computes the `Sol3` compact-edge A-dual-orbit modified-Mathieu finite-difference spectrum.  It does not read external CSV files or analysis logs; it locally recomputes hyperbolic `SL(2,Z)` matrix rigidity, dual orbits, grid/window refinement, and Q-level/discriminant robustness.

**Key functions:**

- `compute_tables()`: compute rigidity, lambda, convergence, compact-vs-noncompact, and robustness rows.
- `gate_rows(rigidity_rows, lambda_rows, convergence_rows, compare_rows, robustness_rows)`: build witness gate rows.

## Dependencies

- SymPy: exact symbolic closure/Hessian audit.
- NumPy/SciPy: modified-Mathieu finite-difference eigenvalue computation.
- [dirac](../dirac/README.md): scalar-axis DBR closure results.
- [topology](../topology/README.md): class-A structure data.

