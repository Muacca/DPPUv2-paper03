# Dirac Layer

-> [English](README.md)

paper07 witness が使う Dirac constraint-closure engine 群です。Branch-A local first-order Palatini model、scalar-axis projection、DBR closure、Routh-before-Dirac route comparison を、外部解析ログや中間 CSV に依存せず再計算します。

## 概要

`dppu.dirac` は論文本文で使う DBR/RBD route 判定の計算核です。主な利用元は `scripts/paper07/rank_defect_activation_witness.py`、`scripts/paper07/harmonic_kernel_correspondence.py`、`scripts/paper07/route_comparison_s3_branch_a_witness.py` です。

## モジュール

### `constraint_closure.py`

Branch-A local first-order Lagrangian と harmonic projection から Dirac constraint closure を構成します。

**主な関数:**

- `build_branchA_local_L2()`: Branch-A local quadratic Lagrangian を構成します。
- `project_lagrangian(local, sample_l, l_sector, projection_variant="x")`: local model を harmonic sector に射影します。
- `build_dirac_for_projection(proj, local)`: primary/secondary constraints、Poisson rank、first/second class split を計算します。
- `build_dirac_for_second_order_lagrangian(L_red, q_red, v_red)`: second-order reduced Lagrangian の Dirac closure を計算します。
- `sample_geometry(sample_l, variant="x")`: route comparison 用の代表 harmonic geometry を返します。

### `local_models.py`

`T3`, `Nil3`, `Sol3`, `S3` の class-A local model builder です。rank-defect activation と harmonic-kernel witness がここから topology-specific local model を作ります。

**主な関数:**

- `build_branchA_local_L2_with_lc(topology)`: topology ごとの LC structure を入れた Branch-A local model。
- `build_t3_local()`: `T3` local model。
- `build_sol3_local()`: `Sol3` local model。
- `build_full_local_projection(local)`: full local projection helper。

### `scalar_axis.py`

scalar-axis DBR closure engine です。`rank(A_DBR)=rank(n)` と active/inactive scalar-axis の reduced Hessian 接触を確認する witness から利用します。

**主な関数:**

- `build_scalar_axis_projection(local, axis)`: `x/y/z` scalar-axis projection。
- `dirac_standard(proj, local, rank_subs)`: DBR standard closure。
- `reduced_hamiltonian_status(result)`: reduced Hamiltonian の物理 pair status。

### `route_variants.py`

Routh-before-Dirac route の negative-control helper です。DBR と同じ Branch-A local model から RBD variants を再計算します。

**主な関数:**

- `routh_eliminate(local, also_K0)`: spatial-omega または full-K Routh elimination。
- `project_reduced_dirac(local, routh, sample_l)`: Routh 後の reduced Dirac closure。

## 依存関係

- SymPy: symbolic Poisson bracket、constraint rank、Hamiltonian/Hessian operations。
- [action](../action/README_ja.md): Palatini/EH density helpers。
- [topology](../topology/README_ja.md): class-A structure data。
- [spectral](../spectral/README_ja.md): reduced-Hessian and harmonic-basis audits。

