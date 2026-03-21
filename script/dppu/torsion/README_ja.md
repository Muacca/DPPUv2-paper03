# Torsion Layer

⇒ [English](README.md)

トーション（ねじれ）構造を定義するモジュール群。

## 概要

トーションansatz、Nieh-Yan密度の計算を提供。

## モジュール

### mode.py

トーションモードの定義。

**Mode Enum:**

| モード | 説明 | パラメータ |
|--------|------|-----------|
| `AX` | 軸性トーションのみ | η |
| `VT` | ベクトルトレースのみ | V |
| `MX` | 混合モード | η, V |

```python
from dppu.torsion import Mode

mode = Mode.MX  # 混合モード
```

### nieh_yan.py

Nieh-Yanバリアントの定義とNY密度計算。

**NyVariant Enum:**

| バリアント | 定義 | 物理的意味 |
|-----------|------|-----------|
| `TT` | T^a ∧ T_a | トーション-トーション項 |
| `REE` | e^a ∧ e^b ∧ R_ab | 曲率由来項 |
| `FULL` | TT - REE | 完全なNieh-Yan密度 |

**主要関数:**

- `compute_nieh_yan_TT(T)`: TT項を計算
- `compute_nieh_yan_REE(R)`: REE項を計算
- `compute_nieh_yan_full(T, R)`: 完全なNY密度

### ansatz.py

トーションテンソルのansatz構成。

**T1: 完全反対称（軸性）**

```
T^{ijk} = (2η/r) ε^{ijk}  (i,j,k ∈ {0,1,2})
```

空間方向のみに存在し、軸性ベクトルS^μを定義:
```
S^μ = (1/6) ε^{μνρσ} T_{νρσ}
```

**T2: ベクトルトレース**

```
T^{abc} = (1/3)(η^{ac}T^b - η^{ab}T^c)
T^μ = (0, 0, 0, V)
```

S¹方向のみに存在。

**主要関数:**

- `build_torsion_tensor(mode, eta, V, r)`: ansatzからT^{abc}を構成

### scalar.py

トーションスカラーの計算。

**トーションスカラー T:**

```
T = T_{abc} T^{abc}
```

## 使用例

```python
from dppu.torsion import Mode, NyVariant, build_torsion_tensor

# トーションテンソル構成
T = build_torsion_tensor(Mode.MX, eta=-1.0, V=2.0, r=1.0)

# Nieh-Yan密度
from dppu.torsion import compute_nieh_yan_full
N = compute_nieh_yan_full(T, R)
```

## 物理的意味

- **軸性トーション (AX)**: スピン-スピン相互作用に関連
- **ベクトルトーション (VT)**: 質量生成機構に関連
- **混合モード (MX)**: 両方の効果を含む一般的設定

## 依存関係

- [utils](../utils/README_ja.md): Levi-Civita記号

## 関連モジュール

- [connection](../connection/README_ja.md): トーションからcontortionを計算
- [curvature](../curvature/README_ja.md): Nieh-Yan項のREE部分
