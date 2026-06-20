# DPPUv2 計算エンジン v5 — スクリプトディレクトリ

⇒ [English](README.md)

**論文**: "Unified Geometric Landau EFT of Homogeneous S³×S¹ Minisuperspace in Einstein-Cartan + Nieh-Yan Theory"（paper03）

Einstein-Cartan + Nieh-Yan 理論における一様 S³×S¹ minisuperspace の数値・記号計算のための Python パッケージ群と実行スクリプト。

---

## ディレクトリ構成

```
script/
├── dppu/                      # メイン Python パッケージ（DPPUv2 Engine v5）
│   ├── geometry/              # 計量・体積形式・構造定数
│   ├── connection/            # Levi-Civita 接続・Contortion・EC 接続
│   ├── curvature/             # Riemann・Ricci・Hodge 双対・自己双対性・Pontryagin・Weyl
│   ├── torsion/               # トーションモード・Ansatz・Nieh-Yan 密度
│   ├── action/                # ラグランジアン・有効ポテンシャル・安定性分類
│   ├── topology/              # 統一エンジン（S³×S¹・T³×S¹・Nil³×S¹）
│   ├── engine/                # 計算パイプライン・ロギング・チェックポイント
│   ├── scanning/              # パラメータスキャン・Squashing 走査・Phase1 統合
│   └── utils/                 # 共通ユーティリティ（Levi-Civita 記号・記号計算・可視化）
│
└── scripts/                   # 実行スクリプト
    ├── paper03/               # paper03 固有の解析スクリプト
    ├── proofs/                # 解析的・記号的証明スクリプト
    └── visualize/             # インタラクティブビューア
```

### `docs/` - ドキュメント

技術ドキュメントと規約：
- [DPPUv2 Engine CONVENTIONS](docs/CONVENTIONS_ja.md) — エンジンコアの規約と仕様
- [DPPUv2 SymPy guideline](docs/SymPy_guideline_ja.md) — SymPy 使用ガイドラインとベストプラクティス

---

## パッケージ概要（dppu/）

| モジュール | 役割 | 主要クラス・関数 |
|-----------|------|----------------|
| [`geometry`](dppu/geometry/README_ja.md) | 計量・フレーム場定義 | `build_metric`, `frame_field` |
| [`connection`](dppu/connection/README_ja.md) | EC 接続の構築 | `levi_civita`, `contortion`, `ec_connection` |
| [`curvature`](dppu/curvature/README_ja.md) | 曲率テンソル群・SD 診断・Pontryagin | `RiemannTensor`, `SDExtensionMixin`, `compute_pontryagin_inner_product` |
| [`torsion`](dppu/torsion/README_ja.md) | トーション構造 | `Mode`, `NyVariant`, `build_torsion_tensor` |
| [`action`](dppu/action/README_ja.md) | 作用・安定性解析 | `build_lagrangian`, `classify_stability` |
| [`topology`](dppu/topology/README_ja.md) | 全トポロジー・全 DOF 対応統一エンジン | `UnifiedEngine`, `DOFConfig`, `TopologyType`, `FiberMode` |
| [`engine`](dppu/engine/README_ja.md) | 15 ステップ計算パイプライン | `BaseFrameEngine`, `ComputationLogger`, `CheckpointManager` |
| [`scanning`](dppu/scanning/README_ja.md) | パラメータ空間走査 | `run_scan`, `Phase1ResultsLoader`, `SDDiagnosticsWithPhase1` |
| [`utils`](dppu/utils/README_ja.md) | 共通ユーティリティ | `LeviCivita`, `symbolic_simplify`, `visualize_potential` |

---

## 実行スクリプト概要（scripts/）

### paper03/ — paper03 固有解析

#### Pontryagin 保護

| スクリプト | 説明 |
|-----------|------|
| `pontryagin_scan_s3.py` | S³ での Pontryagin スキャン（トーションモード・パラメータ空間） |
| `pontryagin_scan_nil3.py` | Nil³ での Pontryagin スキャン |
| `pontryagin_anisotropic_twist.py` | 異方性 twist 下での Pontryagin 検証 |
| `torsion_mode_analysis.py` | AX / VT / MX モードの Pontryagin 分類 |

#### Spin-2 セクター

| スクリプト | 説明 |
|-----------|------|
| `spin2_schur_degeneracy.py` | 等方点での Schur 5 重縮退 μ²_q = 48π² の検証 |
| `spin2_spectrum_splitting.py` | SO(3)→SO(2) 分裂と相転移 ε_c± |
| `spin2_diagonal_hessians.py` | 有限差分による H₁₁^q(ε)・H₂₂^q(ε) |
| `spin2_offdiag_hessian.py` | H₄₄^q（T₄・T₅ off-diagonal shear） |

#### Spin-1 セクター

| スクリプト | 説明 |
|-----------|------|
| `spin1_kinetic_metric.py` | velocity-mode engine による field-space metric G（約 300–600 秒） |
| `spin1_null_mode_analysis.py` | null mode 分類（Case B: non-dynamical auxiliary） |
| `spin1_cross_hessian.py` | cross-Hessian：spin-1 と spin-2 / spin-0 の decoupling 検証 |
| `spin1_auxiliary_exclusion.py` | Case C 排除の検証 |

#### Spin-0 セクター

| スクリプト | 説明 |
|-----------|------|
| `spin0_rLetaV_hessian.py` | Spin-0 Hessian（r, η, V ブロック） |
| `spin0_reta_block_eigenvalues.py` | Spin-0（r, η）ブロック固有値 |
| `spin0_torsion_cross_hessian.py` | cross-Hessian：spin-0 torsion セクター |

#### EFT

| スクリプト | 説明 |
|-----------|------|
| `eft_invariant_generation.py` | EFT 不変量基底の生成（Landau EFT 選択則） |
| `homogeneous_consistency.py` | 一様幾何整合性チェック |

### proofs/ — 解析的・記号的証明

| スクリプト | 証明内容 |
|-----------|---------|
| `analytic_proof_P_equals_zero.py` | AX/VT モードで P = ⟨R, *R⟩ ≡ 0 の代数的証明 |
| `symbolic_block_all_topologies.py` | 全トポロジー（S³・T³・Nil³）記号的ブロック構造検証 |
| `plucker_cancellation_proof.py` | **Pontryagin 保護（AX）**：Plücker 型消去証明 |
| `pfaffian_decomposability.py` | **Pontryagin 保護（VT）**：Pfaffian 分解可能性証明 |
| `pontryagin_block_decomposition.py` | Pontryagin ブロック分解（P_SS・P_S3・P_33 の独立消去） |
| `riemann_block_classification.py` | Riemann ブロック分類（R_{i3j3} ≠ 0 の確認） |
| `vt_mode_cancellation.py` | VT モード消去機構の証明 |
| `analytical_proof.py` | Lemma 1：C²(ε=0)=0・∂C²/∂ε|₀=0・∂²C²/∂ε²|₀>0（SymPy 記号的） |

### visualize/ — インタラクティブビューア

| ファイル | 内容 |
|---------|------|
| `DPPUv2_interactive_viewer_p03.py` | paper03 用インタラクティブポテンシャルビューア（Jupyter、`%matplotlib widget` が必要） |

---

## クイックスタート

### 依存関係のインストール

```bash
pip install sympy numpy scipy pandas tqdm matplotlib
```

### UnifiedEngine でエンジンを実行

```python
from dppu.topology.unified import UnifiedEngine, DOFConfig, TopologyType, FiberMode
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant

# S³×S¹ squash（paper02 の S3S1Engine と等価）：
cfg = DOFConfig(
    topology=TopologyType.S3,
    enable_squash=True,
    torsion_mode=Mode.AX,
    ny_variant=NyVariant.FULL,
)
engine = UnifiedEngine(cfg)
engine.run()

Veff = engine.data['potential']
fp   = engine.get_free_params()   # 有効な SymPy Symbol の dict
```

### レガシーエンジンのプリセット名を使用

```python
cfg = DOFConfig.from_engine('S3S1Engine')    # S³×S¹ squash
cfg = DOFConfig.from_engine('T3S1Engine')    # T³×S¹
cfg = DOFConfig.from_engine('Nil3S1Engine')  # Nil³×S¹ squash
cfg.torsion_mode = Mode.MX
engine = UnifiedEngine(cfg)
engine.run()
```

### Pontryagin 保護の証明スクリプト

```bash
# AX/VT モードで P = 0 の代数的証明
python scripts/proofs/analytic_proof_P_equals_zero.py

# Plücker 型消去（AX モード保護機構）
python scripts/proofs/plucker_cancellation_proof.py

# Pfaffian 分解可能性（VT モード保護機構）
python scripts/proofs/pfaffian_decomposability.py

# 全トポロジーの記号的ブロック構造検証
python scripts/proofs/symbolic_block_all_topologies.py
```

### Spin-2 解析

```bash
# 等方点での Schur 縮退
python scripts/paper03/spin2_schur_degeneracy.py

# SO(3)→SO(2) 分裂と相転移 ε_c±
python scripts/paper03/spin2_spectrum_splitting.py

# 有限差分による H₁₁^q(ε)・H₂₂^q(ε)
python scripts/paper03/spin2_diagonal_hessians.py
```

### Spin-1 解析

```bash
# field-space metric G（velocity モード、約 300–600 秒）
python scripts/paper03/spin1_kinetic_metric.py

# null mode 分類
python scripts/paper03/spin1_null_mode_analysis.py

# cross-Hessian decoupling 検証
python scripts/paper03/spin1_cross_hessian.py
```

---

## paper03 の主要結果と対応スクリプト

| 結果 | 論文箇所 | 内容（証明種別） | 対応スクリプト |
|-----|---------|----------------|--------------|
| **Pontryagin 保護定理** | §3 | AX/VT で P ≡ 0 — **解析的証明**（Plücker/Pfaffian 消去） | `proofs/plucker_cancellation_proof.py`、`proofs/pfaffian_decomposability.py` |
| **MX モードで P ≠ 0** | §3 | P₀ = 2Vη(V²r²+9η²−36)/(9r³)；源は捩れ混合 Vη、twist ではない | `paper03/torsion_mode_analysis.py` |
| **Schur 5 重縮退** | §4 | H_AB = μ²_q δ_AB、μ²_q = 48π² — **記号的証明** | `paper03/spin2_schur_degeneracy.py` |
| **相転移 ε_c±** | §4 | ε_c+ ≈ 0.483（squash）、ε_c− ≈ −0.293（shear）— **数値解析** | `paper03/spin2_spectrum_splitting.py`、`paper03/spin2_diagonal_hessians.py` |
| **Spin-1 physical/auxiliary 分解** | §5 | X_k 物理的・Y_k 補助的（Case B 確定）— **解析的 + 数値** | `paper03/spin1_kinetic_metric.py`、`paper03/spin1_null_mode_analysis.py` |
| **λ_phys ≈ −0.4183 < 0** | §5 | Schur 消去後の有限一般化固有値：タキオン的（非ゴースト） | `paper03/spin1_kinetic_metric.py` |
| **二次セクター decoupling** | §6 | 等方点で cross-Hessian がゼロ — **数値検証** | `paper03/spin1_cross_hessian.py`、`paper03/spin0_torsion_cross_hessian.py` |
| **EFT 不変量基底** | §6 | Landau EFT 選択則と禁制項 — **群論的** | `paper03/eft_invariant_generation.py` |

---

## トポロジー一覧

| トポロジー | Lie 群 | 背景曲率 | Koszul 公式 |
|-----------|-------|---------|-----------|
| S³×S¹ | SU(2) | +24/r² | 簡略（bi-invariant） |
| T³×S¹ | U(1)³ | 0（平坦） | 自明 |
| Nil³×S¹ | Heisenberg 群 | −1/(2R²) | 一般（bi-invariant 非成立） |

---

## チェックポイント機能

長時間計算の中断・再開が可能：

```bash
# チェックポイント付きで実行
python scripts/paper03/spin1_kinetic_metric.py \
    --output-dir output/ --checkpoint-dir ./checkpoints

# 中断後の再開（チェックポイントディレクトリを再指定するだけ）
python scripts/paper03/spin1_kinetic_metric.py \
    --output-dir output/ --checkpoint-dir ./checkpoints
```

---

## ライセンス

リポジトリルートの LICENSE ファイルを参照してください。

---

**Author**: Muacca
**Version**: DPPUv2 Engine v5
**Date**: 2026-03
