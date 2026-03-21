# Action Layer

⇒ [English](README.md)

作用（action）と安定性解析を担当するモジュール群。

## 概要

ラグランジアン、有効ポテンシャル、安定性分類を提供。

## モジュール

### lagrangian.py

Einstein-Cartan + Nieh-Yanラグランジアンの構成。

**作用:**

```
S = ∫ d⁴x √|g| L
L = R/(2κ²) + θ_NY × N
```

**主要関数:**

- `build_lagrangian(R, N, kappa, theta_NY)`: ラグランジアンを構成
- `integrate_angular(L, topology)`: 角度方向を積分

### potential.py

有効ポテンシャルの計算。

**有効ポテンシャル:**

角度積分後、r依存の有効ポテンシャルV(r)を抽出:

```
S = ∫ dr × V_eff(r)
```

**主要関数:**

- `extract_effective_potential(S_integrated)`: 有効ポテンシャルを抽出
- `find_extrema(V, r_range)`: 極値を探索
- `compute_barrier_height(V, r_min)`: バリア高さを計算

### stability.py

安定性解析と分類。

**安定性タイプ:**

| タイプ | 条件 | 物理的解釈 |
|--------|------|-----------|
| type-I | V''(r*) > 0, V(r*) > 0 | 準安定（バリアあり） |
| type-II | V''(r*) > 0, V(r*) < 0 | 真の極小（真空より低い） |
| type-III | 極小なし or V''(r*) ≤ 0 | 不安定 |

**主要関数:**

- `classify_stability(V, r_star)`: 安定性を分類
- `StabilityResult`: 分類結果のデータクラス

```python
@dataclass
class StabilityResult:
    stability_type: str  # 'type-I', 'type-II', 'type-III'
    r_star: float        # 安定点の位置
    V_star: float        # 安定点でのポテンシャル値
    V_second: float      # 二次微分
    barrier_height: float  # バリア高さ (type-Iのみ)
```

## 使用例

```python
from dppu.action import build_lagrangian, classify_stability

# ラグランジアン構成
L = build_lagrangian(R, N, kappa=1.0, theta_NY=0.5)

# 安定性分類
result = classify_stability(V_eff, r_star=1.5)
print(f"Type: {result.stability_type}")
print(f"r* = {result.r_star:.3f}")
```

## Phase1との統合

Phase1スキャン結果はCSVに保存され、Phase2でのSD解析に使用される:

```csv
topology,ny_variant,V,eta,theta,r0,delta_V,stability_type
S3,FULL,2.0,-1.0,0.0,1.234,0.056,type-I
```

## 依存関係

- [curvature](../curvature/README_ja.md): Ricciスカラー
- [torsion](../torsion/README_ja.md): Nieh-Yan密度

## 関連モジュール

- [scanning](../scanning/README_ja.md): パラメータ空間のスキャン
