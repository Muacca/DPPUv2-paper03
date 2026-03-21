# Connection Layer

⇒ [English](README.md)

接続（connection）の計算を担当するモジュール群。

## 概要

Levi-Civita接続、contortion、Einstein-Cartan接続を計算する。

## モジュール

### levi_civita.py

トーションフリーのLevi-Civita接続を計算。

**主要関数:**

- `compute_christoffel_symbols(structure_constants, metric)`: 一般Koszul公式による計算

**Koszul公式:**

```
Γ^a_{bc} = (1/2)(C^a_{bc} + η^{ad}η_{be}C^e_{dc} - η^{ad}η_{ce}C^e_{bd})
```

bi-invariant計量の場合は簡略化される:
```
Γ^a_{bc} = (1/2)C^a_{bc}
```

### contortion.py

トーションから生じるcontortionテンソルを計算。

**主要関数:**

- `compute_contortion(torsion, metric)`: contortionテンソルを計算

**定義:**

```
K_{abc} = (1/2)(T_{abc} + T_{bca} - T_{cab})
```

### ec_connection.py

Einstein-Cartan接続（LC + contortion）を統合。

**主要関数:**

- `compute_ec_connection(christoffel, contortion)`: EC接続を計算

**定義:**

```
Γ^{EC,a}_{bc} = Γ^{LC,a}_{bc} + K^a_{bc}
```

## 使用例

```python
from dppu.connection import compute_christoffel_symbols, compute_contortion

# Levi-Civita接続
christoffel = compute_christoffel_symbols(C, eta)

# Contortion
K = compute_contortion(T, eta)

# EC接続
Gamma_EC = christoffel + K
```

## 計量整合性

EC接続は計量整合性を満たす:
```
∇_a η_{bc} = 0
```

これは自動的に検証される。

## 依存関係

- [geometry](../geometry/README_ja.md): 計量テンソル
- [torsion](../torsion/README_ja.md): トーションテンソル

## 関連モジュール

- [curvature](../curvature/README_ja.md): 接続から曲率を計算
