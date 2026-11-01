# 計算コードと再現方法

[English](README.md) / [規約](docs/CONVENTIONS_ja.md) / [実装ガイド](docs/IMPLEMENTATION_ja.md)

Paper03 v3（2026-11-01）用。Python 3.10以上。
ECはEinstein–Cartan、NYはNieh–Yanを表す。空間的一様Euclideanモデルを扱い、
完全NYの係数は一定、alpha=0である。[LaTeX原稿](../LaTeX/main.tex) を添付している。

## 図とデータの再現

`script` ディレクトリから実行する。

```powershell
python -m pip install -r requirements.txt
python -B -m paper_results render --output-dir ../rendered
python -B -m paper_results recompute --output-dir ../recomputed-data
python -B -m paper_results render --data-dir ../recomputed-data --output-dir ../recomputed-figures
```

出力先は各々未作成のディレクトリを指定する。`render` は [data/revised](../data/revised) を読み、
[manifest.json](../data/revised/manifest.json) の11データファイルのハッシュを確認して
PNG/PDF各10図を描く。幾何計算は呼ばない。`recompute` は固定格子で再計算し、
11データファイルとmanifestを出力する。掲載図は [LaTeX/figures/v3](../LaTeX/figures/v3) にある。

[依存指定](requirements.txt) は、計算・描画用のNumPy・Matplotlibと、
notebook表示用のIPython・ipywidgetsを含む。
添付コードは `paper_results`、図notebook、保存データviewerで構成する。

## 計算インターフェース

| 入口 | 出力 | 論文との対応 |
|---|---|---|
| [model.cartan](paper_results/model.py) | coframe・接続・完全EC曲率・NY・P・B | 第2節、付録A–B |
| [model.potential](paper_results/model.py) | 静的ポテンシャルUの閉形式 | 式(M3) |
| [model.scalar_data / scalar_bulk](paper_results/model.py) | スカラー勾配・Hessian、roundスカラー円方向密度 | 第3節 |
| [model.twist_bulk](paper_results/model.py) | 静的twist cost・振幅四次の微分項・速度Hessian | 第4節 |
| [model.shear_curvatures / shape_minimum](paper_results/model.py) | 経路曲率と制限doublet極小 | 第5節 |
| [circle.joint_bulk](paper_results/circle.py) | 一般bulk密度、静的・速度二次・速度一次の寄与 | 式(M4)–(M6)、付録J |
| [circle.diagonal_cross_coefficients](paper_results/circle.py) | 五つの速度交叉係数 | 式(X2) |
| [circle.coupled_periodic_family](paper_results/circle.py) | 周期shape–twist系列の固定求積 | 式(X4) |
| [data](paper_results/data.py) / [figures](paper_results/figures.py) | 保存配列・manifest・図生成 | 付録H |

`paper_results` の `eta,V` は幾何学的torsion係数である。`joint_bulk` は
r・L・M・omegaの値と一階微分を受け取る。Mは行列式1、
`tr(inv(M) @ M_prime)=0` を満たす。出力密度は周期全微分を除いたbulk密度である。
静的最小化と周期profileの変分には、論文でそれぞれ指定した定義域を用いる。

## Notebookとviewer

[図notebook](DPPUv2_Paper03_Figures_v3.ipynb) をJupyterで開き、
作業ディレクトリをリポジトリ直下または `script` とする。
既定の `RECOMPUTE=False` は保存データを描画し、`True` を明示すると再計算する。
実行ごとに `script` と同じ階層へ日時付きの新しいフォルダ（`notebook-figures-YYYYMMDD-HHMMSS`、再計算時は `notebook-data-...` も）を作る。
notebookは同じディレクトリの [保存データviewer](paper03_saved_data_viewer.py) も読み込む。
`show(data_dir=None)` はdropdownと出力widgetを返し、静的断面・速度交叉係数・
周期作用系列など8種類の保存データ表示を切り替える。新たな幾何計算や格子変更は行わない。
Jupyter本体は利用環境に必要であり、widgetはIPythonとipywidgetsを使う。

規約は [CONVENTIONS](docs/CONVENTIONS_ja.md) にまとめた。APIと戻り値は [実装ガイド](docs/IMPLEMENTATION_ja.md) に示す。
