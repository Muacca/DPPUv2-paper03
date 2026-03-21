# Geometry Layer

⇒ [English](README.md)

幾何学的基盤を提供するモジュール群。

## 概要

計量テンソル、体積形式など、微分幾何学の基本構造を定義する。

## モジュール

### metric.py

計量テンソルの定義と操作。

**主要クラス・関数:**

- `create_frame_metric(topology)`: トポロジーに応じたフレーム計量を生成
- `get_volume_factor(topology)`: 体積因子を計算

**計量の形式:**

```
η_ab = diag(-1, +1, +1, +1)  (Lorentzian signature)
```

ただし、コンパクト化後は有効的にEuclidean signatureとなる。

## トポロジー別計量

| トポロジー | 計量構造 | 体積 |
|-----------|----------|------|
| S³×S¹ | bi-invariant on SU(2) | 2π²Lr³ |
| T³×S¹ | flat | (2π)⁴LR₁R₂R₃ |
| Nil³×S¹ | left-invariant on Heisenberg | (2π)⁴LR³ |

## 使用例

```python
from dppu.geometry import create_frame_metric

# S³用のフレーム計量
metric = create_frame_metric('S3')
```

## 依存関係

- SymPy (シンボリック計算)
- NumPy (数値評価)

## 関連モジュール

- [connection](../connection/README_ja.md): 計量から接続を計算
- [topology](../topology/README_ja.md): トポロジー固有の計量定義
