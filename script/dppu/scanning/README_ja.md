# Scanning Layer

⇒ [English](README.md)

パラメータ空間の走査とPhase1統合を担当するモジュール群。

## 概要

(V, η, θ)パラメータ空間のスキャン、Phase1結果の読み込み・補間、SD診断との統合を提供。

## モジュール

### parameter_scan.py

パラメータ空間のグリッドスキャン。

**主要関数:**

- `run_scan(...)`: パラメータ空間をスキャン

**使用例:**

```python
from dppu.scanning import run_scan

results = run_scan(
    V_points=50, eta_points=100, theta_points=20,
    V_min=0.0, V_max=5.0,
    eta_min=-5.0, eta_max=5.0,
    theta_min=0.0, theta_max=5.0,
    topologies=['S3', 'T3', 'Nil3'],
    ny_variants=['FULL', 'TT', 'REE'],
    output_dir='output/',
    n_workers=8
)
```

**出力CSV形式:**

| カラム | 説明 |
|--------|------|
| topology | S3, T3, Nil3 |
| ny_variant | FULL, TT, REE |
| V, eta, theta | パラメータ値 |
| r0 | 安定点の位置 (None if type-III) |
| delta_V | バリア高さまたは井戸深さ |
| stability_type | type-I, type-II, type-III |

### potentials.py

トポロジー別の有効ポテンシャル関数。

**POTENTIAL_FUNCTIONS辞書:**

```python
from dppu.scanning import POTENTIAL_FUNCTIONS

# トポロジー別ポテンシャル取得
V_func = POTENTIAL_FUNCTIONS['S3_FULL']
V_value = V_func(r=1.0, eta=-1.0, V_param=2.0, theta=0.0)
```

### phase1_loader.py

Phase1スキャン結果の読み込みと補間。

**Phase1ResultsLoader:**

```python
from dppu.scanning import Phase1ResultsLoader

# CSVから読み込み
loader = Phase1ResultsLoader.from_csv(
    'output/dppu_scan_S3_FULL.csv',
    theta_fixed=0.0
)

# サマリー表示
print(loader.summary())
# {'total_points': 5000, 'stable_points': 3200, ...}

# r*の補間取得
r_star = loader.get_r_star(V=2.0, eta=-1.0)

# 安定性タイプ取得
phase = loader.get_phase_type(V=2.0, eta=-1.0)
# 'I', 'II', or 'III'
```

**グリッドタイプ:**

- `regular`: 規則グリッド → RegularGridInterpolator使用
- `irregular`: 不規則グリッド → LinearNDInterpolator使用

### sd_audit.py

SD診断とPhase1の統合。

**SDDiagnosticsWithPhase1:**

```python
from dppu.scanning import SDDiagnosticsWithPhase1, Phase1ResultsLoader

loader = Phase1ResultsLoader.from_csv('scan.csv', theta_fixed=0.0)
sd_diag = SDDiagnosticsWithPhase1(engine, phase1_loader=loader)

# r*でのSD評価
result = sd_diag.evaluate_at_rstar(V=2.0, eta=-1.0, theta_NY=0.0)
print(f"r* = {result.r_star:.3f}")
print(f"SD residual = {result.sd_residual:.4f}")

# (η, V)平面のスキャン
results = sd_diag.scan_with_phase1_rstar(
    eta_range=(-5, 3, 100),
    V_range=(0.5, 5, 50),
    theta_NY=0.0
)

print(f"SD curve points: {len(results['sd_curve'])}")
print(f"Type I ∩ SD: {len(results['type_I_sd_intersection'])}")
```

**SDScanResult:**

```python
@dataclass
class SDScanResult:
    eta: float
    V: float
    r_star: float
    phase_type: str
    sd_residual: float
    asd_residual: float
    curvature_norm: float
    is_nontrivial_sd: bool
    is_nontrivial_asd: bool
```

## ワークフロー

### Phase1 → Phase2 統合

1. **Phase1**: パラメータスキャンで安定点r*(V, η)を特定
2. **Phase2**: r*でのSD診断を実行
3. **分析**: SD曲線と安定領域の交差を調査

```bash
# Phase1: パラメータスキャン
python scripts/run_parameter_scan.py --output-dir output/

# Phase2: SD解析
python scripts/run_sd_analysis.py \
    --phase1-csv output/dppu_scan_S3_FULL.csv \
    --theta 0.0
```

## 依存関係

- [topology](../topology/README_ja.md): 計算エンジン
- [curvature](../curvature/README_ja.md): SD診断
- pandas (CSV読み込み)
- scipy (補間)

## 関連モジュール

- [action](../action/README_ja.md): 安定性分類
