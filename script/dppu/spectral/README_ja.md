# Spectral Layer

-> [English](README.md)

paper07 witness が使う harmonic basis、constraint-surface reduced Hessian、Sol3 compact-edge spectrum の計算 engine 群です。

## 概要

`dppu.spectral` は、閉 1-form 空間と `ker n` の対応、DBR closure 後の reduced Hessian kernel 接触、ならびに `Sol3` compact edge の A-dual-orbit modified-Mathieu 正性 witness を提供します。主な利用元は `scripts/paper07/harmonic_kernel_correspondence.py` と `scripts/paper07/spectral_completion_witness.py` です。

## モジュール

### `harmonic_basis.py`

class-A structure constants から left-invariant 1-form の closure matrix と harmonic/co-closed audit を構成します。

**主な関数:**

- `structure_constants(n)`: class-A structure matrix `n` から構造定数を構成します。
- `closure_matrix(n)`: closed 1-form 条件の行列を返します。
- `basis_matrix(nullspace)`: nullspace basis を matrix form に整形します。
- `span_equal(left, right)`: 2 つの span が一致するか確認します。
- `co_closed_audit(constants, k_form_basis)`: co-closed/harmonic status を監査します。
- `audit_case(topology, bianchi, n, case_id)`: topology ごとの harmonic-basis audit。
- `generic_identity_check()`: formal structure-map identity の記法チェック。

### `reduced_hessian.py`

DBR closure 後の constraint-surface reduced phase space 上で Hessian kernel を監査します。

**主な関数:**

- `reduced_hessian_pair_audit(result)`: reduced symplectic pair と Hamiltonian Hessian kernel の対応を確認します。
- `scalar_axis_hessian_audits(topology, local, active_axes)`: scalar-axis ごとの DBR Hessian audit。

### `sol3_compact_edge.py`

`Sol3` compact edge の A-dual-orbit modified-Mathieu finite-difference spectrum を計算します。外部 CSV や解析ログは読み込まず、SL(2,Z) hyperbolic matrix、dual orbit、grid/window refinement、Q-level/discriminant robustness をローカルに再計算します。

**主な関数:**

- `compute_tables()`: rigidity、lambda、convergence、compact-vs-noncompact、robustness rows を計算します。
- `gate_rows(rigidity_rows, lambda_rows, convergence_rows, compare_rows, robustness_rows)`: witness gate rows を構成します。

## 依存関係

- SymPy: exact symbolic closure/Hessian audit。
- NumPy/SciPy: modified-Mathieu finite-difference eigenvalue computation。
- [dirac](../dirac/README_ja.md): scalar-axis DBR closure results。
- [topology](../topology/README_ja.md): class-A structure data。

