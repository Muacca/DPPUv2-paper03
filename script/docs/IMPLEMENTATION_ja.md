# 数値計算の実装ガイド

[English](IMPLEMENTATION.md) / [規約](CONVENTIONS_ja.md) / [再現方法](../README_ja.md)

## 論文の計算

コマンドライン入口と同じく `script` からPythonを実行する。
`paper_results` はNumPy配列と閉形式を用いる。
`cartan` はC・幾何学的T・LC/EC接続・一定成分の曲率を組み立ててから縮約する。
`potential` は閉形式(M3)を評価する。

```python
import numpy as np
from paper_results.model import cartan, potential
from paper_results.circle import joint_bulk

geometry = cartan(r=3., L=1., eta=1., V=1.)
U = potential(3., eta=1., V=1.)
assert np.isclose(-geometry["volume"] * geometry["R_scalar"] / 2, U)
assert geometry["torsion_residual"] < 1e-12
bulk = joint_bulk(r=3., M=np.eye(3), M_prime=np.zeros((3, 3)))
assert np.isclose(bulk["bulk"], U)
```

配列は浮動小数点数なので、数値残差には許容誤差を用いる。厳密な恒等式は論文に示す。
`cartan` は静的計算用である。周期profileには、全微分を積分除去した式(M6)を実装する
`circle.joint_bulk` を使い、`bulk` をyの区間[0,1]で積分する。
このbulk関数の点ごとの値は、全微分を除く前の局所曲率密度と一致するとは限らない。

| 関数 | 戻り値 |
|---|---|
| `cartan` | coframe, C, Gamma, T, R, R_LC, R_scalar, TT, REE, NY, P, B, volume, torsion_residualなどの辞書 |
| `scalar_data` | (r,L,eta,V)順の勾配と4×4 Hessian |
| `shape_minimum` | 非負のdoublet振幅と最小shape関数 |
| `twist_bulk` | 静的twist寄与、微分寄与、3×3速度Hessian |
| `scalar_bulk` | roundスカラーbulk密度 |
| `joint_bulk` | static, velocity_quadratic, velocity_linear, bulk, gamma, lapse, shift |
| `diagonal_cross_coefficients` | G_rw, G_sw, G_rs, G_rL, G_sL |
| `coupled_periodic_family` | 周波数系列、静的・微分・全作用の配列、C_A、round参照値 |

主要な作用関数は正の尺度を確認する。実有限配列と行列式1のshapeを使い、
`joint_bulk` ではさらに配列次元とM_primeの接条件を確認する。
行列指数に使うshear基底は対称tracelessである。

## データと描画

`data.recompute(output)` は新しい出力先に固定格子のデータを書く。
`data.load_data(source)` はschemaと保存データのハッシュを確認し、CSV配列とスカラーJSONを返す。
`figures.render(source,output)` はこの配列から10組の図を描く。
manifestはパラメータ・格子・model.py/circle.py・データのハッシュを記録する。
読み込み時に照合するのはデータのハッシュで、現在のsourceとの自動照合は行わない。
明示的な再計算によって、新しいsourceとデータの対応が記録される。
