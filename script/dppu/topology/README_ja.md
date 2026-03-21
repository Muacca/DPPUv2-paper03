# Topology Layer

⇒ [English](README.md)

トポロジー固有の計算エンジンを提供するモジュール群。

## 概要

Topology Layer は、EC+NY 15ステップ計算パイプライン内で $M^3 \times S^1$ 多様体（$M^3$ = S³, T³, Nil³）の幾何を実装する。Paper 03 では全トポロジーエンジンを統一設定インターフェース（`UnifiedEngine` + `DOFConfig`）に統合し、paper 02 の個別クラス（`S3S1Engine`, `T3S1Engine`, `Nil3S1Engine`）を置き換えた。

## クラス階層

```
BaseFrameEngine          (dppu/engine/pipeline.py)
    └── TopologyEngine   (base_topology.py)   ← abstract
            ├── S3Geometry    (s3.py)
            ├── T3Geometry    (t3.py)
            └── Nil3Geometry  (nil3.py)
```

`UnifiedEngine`（`unified.py`）は `DOFConfig` インスタンスに基づき適切なサブクラスを選択するファクトリ。

---

## モジュール

### unified.py — UnifiedEngine と DOFConfig

全 paper 03 計算のエントリポイント。

**主要クラス:**

- `UnifiedEngine`: `DOFConfig` から適切なトポロジーエンジンを生成するファクトリ
- `DOFConfig`: 有効な幾何 DOF を指定する設定データクラス
- `TopologyType`: Enum — `S3` / `T3` / `NIL3`
- `FiberMode`: Enum — `NONE` / `TWIST` / `MIXING` / `BOTH`

**DOFConfig パラメータ:**

| パラメータ | 型 | 意味 |
|---|---|---|
| `topology` | `TopologyType` | 基底トポロジー |
| `enable_squash` | `bool` | ε 変形を有効化 |
| `enable_shear` | `bool` | s (T₂) shear を有効化 |
| `fiber_mode` | `FiberMode` | ファイバー DOF 選択 |
| `isotropic_twist` | `bool` | 全方向共通の単一 ω を使用 |
| `torsion_mode` | `Mode` | `AX` / `VT` / `MX` |
| `ny_variant` | `NyVariant` | Nieh-Yan バリアント |
| `enable_velocity` | `bool` | G 計量計算用の速度シンボルを有効化 |

**使用例:**

```python
from dppu.topology.unified import UnifiedEngine, DOFConfig, TopologyType, FiberMode
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant

cfg = DOFConfig(
    topology=TopologyType.S3,
    enable_squash=False,
    fiber_mode=FiberMode.BOTH,
    isotropic_twist=False,
    torsion_mode=Mode.MX,
    ny_variant=NyVariant.FULL,
)
engine = UnifiedEngine(cfg)
engine.run()

Veff = engine.data['potential']
fp   = engine.get_free_params()
```

**プリセット名（後方互換性）:**

```python
cfg = DOFConfig.from_engine('S3S1Engine')    # S³, squash=True, fiber=NONE
cfg = DOFConfig.from_engine('T3S1Engine')    # T³, fiber=NONE
cfg = DOFConfig.from_engine('Nil3S1Engine')  # Nil³, squash=True, fiber=NONE
```

レガシーエンジン 13 種の完全等価テーブルは `unified.py` の docstring を参照。

---

### base_topology.py — TopologyEngine 抽象基底クラス

ステップ E4.1/E4.2 のテンプレートメソッドパターンを実装する抽象基底クラス。

**抽象メソッド（サブクラスが実装すること）:**

- `_build_radial_and_deformation_params(params)` — r/R, ε, s シンボルを定義
- `_build_structure_constants(params, C)` — 構造定数配列 C を充填
- `_compute_volume(params)` — 体積因子式を返す

**共通コンクリートヘルパー:**

- `_build_fiber_params(params)` — ω / δ / cd / sd シンボルを追加（TWIST/MIXING/BOTH モード）
- `_add_s3_twist_C(C, ...)` — twist 項 $C^3_{jk}$ を埋め込む
- `_apply_s3_mixing_rotation_to_C(C, params)` — mixing 回転 $M(\delta_0,\delta_1,\delta_2)$ を適用
- `get_free_params()` — 有効な SymPy Symbol の dict を返す
- `get_riemann_lambdified()` — 数値スキャン用の lambdified $R_{abcd}$ を返す

---

### s3.py — S³×S¹ (S3Geometry)

**数学的構造:**

- Lie 群: SU(2)
- 構造定数: $C^i_{jk} = (4/r)\,\lambda_i\,\varepsilon_{ijk}$（ε, s でスケーリング）
- 計量: bi-invariant
- 背景曲率: $R_{\rm LC} = 24/r^2$

**サポート DOF:** ε (squash), s (shear), q₃,q₄,q₅ (off-diagonal shear), TWIST (ω₁,ω₂,ω₃), MIXING (δ₀,δ₁,δ₂)

**体積:** $V = 2\pi^2 L r^3$

**使用例:**

```python
from dppu.topology.s3 import S3Geometry
from dppu.topology.unified import DOFConfig, TopologyType

cfg = DOFConfig.from_engine('S3S1Engine')
engine = S3Geometry(cfg)
engine.run()
```

---

### t3.py — T³×S¹ (T3Geometry)

**数学的構造:**

- Lie 群: U(1)³ (Abelian)
- 構造定数: 全てゼロ
- 計量: flat
- 背景曲率: $R_{\rm LC} = 0$

**体積:** $V = (2\pi)^4 L r^3$（等方的スケーリング $R_1 = R_2 = R_3 = r$）

---

### nil3.py — Nil³×S¹ (Nil3Geometry)

**数学的構造:**

- Lie 群: Heisenberg 群
- 構造定数: $[E_0, E_1] = (1/R)E_2$
- 計量: left-invariant（**NOT bi-invariant**）
- 背景曲率: $R_{\rm LC} = -1/(2R^2)$

**重要:** Nil³ は bi-invariant 計量を持たないため、一般 Koszul 公式を使用（CONVENTIONS §5 参照）。

**体積:** $V = (2\pi)^4 L R^3$

---

## トポロジー比較

| プロパティ | S³×S¹ | T³×S¹ | Nil³×S¹ |
|----------|-------|-------|---------|
| 構造定数 | $\varepsilon_{ijk}$ | 0 | $[E_0,E_1]=E_2$ |
| 背景曲率 | +24/r² | 0 | −1/(2R²) |
| bi-invariant | Yes | Yes | **No** |
| Koszul 公式 | 簡略（bi-invariant） | 自明 | 一般 |
| squash サポート | Yes | No | Yes |
| fiber サポート | TWIST/MIXING/BOTH | TWIST | TWIST |

## 依存関係

- [engine](../engine/README_ja.md): `BaseFrameEngine`
- [geometry](../geometry/README_ja.md): 計量定義
- [connection](../connection/README_ja.md): 接続計算
- [curvature](../curvature/README_ja.md): 曲率計算
- [torsion](../torsion/README_ja.md): トーション構成
- [action](../action/README_ja.md): 作用と安定性
