# Curvature Layer

⇒ [English](README.md)

曲率テンソルの計算と解析を担当するモジュール群。

## 概要

Riemannテンソル、Ricciスカラー、Hodge双対、自己双対性診断、Pontryagin内積、Weylテンソルを提供。

## モジュール

### riemann.py

Riemannテンソルの計算と反対称性検証。

**主要クラス・関数:**

- `compute_riemann_tensor(connection)`: Riemannテンソルを計算
- `verify_antisymmetry(R)`: 反対称性を3レベルで検証
- `RiemannAntisymmetryError`: 反対称性違反時の例外

**反対称性検証（3レベル）:**

1. **Level 1**: SymPyでのシンボリック証明
2. **Level 2**: 高精度数値での反例探索
3. **Level 3**: デバッグモード（開発用）

**Riemann曲率の定義:**

```
R^a_{bcd} = ∂_c Γ^a_{bd} - ∂_d Γ^a_{bc} + Γ^a_{ec}Γ^e_{bd} - Γ^a_{ed}Γ^e_{bc}
```

### ricci.py

Ricciテンソルとスカラーの計算。

**主要関数:**

- `compute_ricci_tensor(R)`: Ricci縮約
- `compute_ricci_scalar(Ric, metric)`: Ricciスカラー

**定義:**

```
R_{ab} = R^c_{acb}
R = η^{ab} R_{ab}
```

### hodge.py

Hodge双対演算子の実装。

**主要関数:**

- `compute_hodge_dual(R)`: R^{ab}_{cd} の Hodge双対を計算
- `cd_block(a, b)`: 6成分ブロック分類

**Hodge双対の定義:**

```
(*R)^{ab}_{cd} = (1/2) ε_{cdef} R^{ab,ef}
```

### self_duality.py

自己双対性（SD/ASD）の診断。

**主要クラス:**

- `SDExtensionMixin`: エンジンにSD機能を追加
- `CurvatureSDDiagnostics`: SD残差の計算と評価

**SD/ASD条件:**

```
R = *R   (自己双対: SD)
R = -*R  (反自己双対: ASD)
```

**主要メトリクス:**

- `sd_residual`: ||R - *R||
- `asd_residual`: ||R + *R||
- `is_sd`: sd_residual < ε
- `is_asd`: asd_residual < ε

### pontryagin.py

Pontryagin内積の計算と診断。

**主要クラス・関数:**

- `compute_pontryagin_inner_product(R)`: P = ⟨R, *R⟩ を計算
- `evaluate_sd_status(params)`: SD状態の完全診断

**Pontryagin 保護定理（paper 03）:**

| トーションモード | 結果 | 証明手法 |
|---|---|---|
| AX（軸性） | P ≡ 0（厳密） | Plücker 型消去 |
| VT（ベクトル-トレース） | P ≡ 0（厳密） | Pfaffian 分解可能性 |
| MX（混合） | 一般に P ≠ 0 | — |

AX・VT モードでは、全トポロジー（S³, T³, Nil³）・全パラメータ値で消去が成立する。MX モードでは:

```
P₀ = 2Vη(V²r² + 9η² − 36) / (9r³)
```

P ≠ 0 の源はトーション混合（Vη 結合）であり、twist ではない。

P = 0 の場合:
```
SD_residual / ||R|| = √2
ASD_residual / ||R|| = √2
```

### weyl.py

Weylテンソル（共形曲率）の計算。

**主要関数:**

- `compute_weyl_tensor(R_abcd, Ricci, R_scalar, metric, dim)`: Weylテンソル C_{abcd} を計算
- `compute_weyl_scalar(C_abcd, metric)`: Weylスカラー C² を計算

**定義:**

```
C_{abcd} = R_{abcd}
         - (1/2)(g_{ac}R_{bd} - g_{ad}R_{bc} - g_{bc}R_{ad} + g_{bd}R_{ac})
         + (R/6)(g_{ac}g_{bd} - g_{ad}g_{bc})
```

**性質:**

- トレースフリー: C^a_{bad} = 0
- Riemannテンソルと同じ対称性

## 使用例

```python
from dppu.curvature import SDExtensionMixin, CurvatureSDDiagnostics

# エンジンにSD拡張を付加
SDExtensionMixin.attach_to(engine)

# SD診断
diag = CurvatureSDDiagnostics(engine)
result = diag.evaluate_sd_status({
    'r': 1.0, 'L': 1.0, 'eta': -1.0, 'V': 2.0,
    'kappa': 1.0, 'theta_NY': 0.0
})

print(f"P = {result['P_RstarR']:.2e}")
print(f"SD/||R|| = {result['sd_residual']/result['curvature_norm']:.4f}")
```

## 依存関係

- [connection](../connection/README_ja.md): EC接続
- [utils](../utils/README_ja.md): Levi-Civita記号
