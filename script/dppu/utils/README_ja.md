# Utils Layer

⇒ [English](README.md)

共通ユーティリティを提供するモジュール群。

## 概要

Levi-Civita記号、シンボリック計算ヘルパー、可視化ユーティリティなど、複数モジュールで使用される共通機能。

## モジュール

### levi_civita.py

Levi-Civita記号（完全反対称テンソル）。

**主要関数:**

- `epsilon_symbol(i, j, k)`: 3次元Levi-Civita記号
- `levi_civita_4d(a, b, c, d)`: 4次元Levi-Civita記号

**3次元:**

```python
from dppu.utils import epsilon_symbol

eps = epsilon_symbol(0, 1, 2)  # = +1
eps = epsilon_symbol(1, 0, 2)  # = -1
eps = epsilon_symbol(0, 0, 2)  # = 0
```

**4次元:**

```python
from dppu.utils import levi_civita_4d

eps = levi_civita_4d(0, 1, 2, 3)  # = +1
eps = levi_civita_4d(1, 0, 2, 3)  # = -1
```

**性質:**

- 完全反対称: 任意の2インデックス交換で符号反転
- 重複インデックスでゼロ
- 正規化: ε₀₁₂₃ = +1（または ε₀₁₂ = +1）

### symbolic.py

シンボリック計算のヘルパー関数。

**主要関数:**

- `prove_zero(expr)`: 式がゼロであることを証明
- `find_nonzero_witness(expr, params)`: 非ゼロの反例を探索
- `generate_test_points(param_ranges, n_points)`: テスト点を生成

**ゼロ証明:**

```python
from dppu.utils import prove_zero
import sympy as sp

x = sp.Symbol('x')
expr = x**2 - x**2

is_zero, proof_type = prove_zero(expr)
# is_zero=True, proof_type='symbolic'
```

**反例探索:**

```python
from dppu.utils import find_nonzero_witness

expr = some_complex_expression
witness = find_nonzero_witness(expr, {'r': (0.1, 5.0), 'eta': (-5, 5)})

if witness:
    print(f"Non-zero at: {witness['params']}")
    print(f"Value: {witness['value']}")
else:
    print("No non-zero witness found")
```

### visualization.py

論文図表用の標準化プロットスタイルとヘルパー関数。

**主要関数:**

- `set_style()`: matplotlib グローバルスタイルを設定（seaborn whitegrid, serif フォント, 300 dpi）
- `save_plot(filename, output_dir="output")`: 標準設定でプロットを保存

**使用例:**

```python
from dppu.utils.visualization import set_style, save_plot
import matplotlib.pyplot as plt

set_style()
plt.plot(x, y)
plt.xlabel(r'$\varepsilon$')
plt.ylabel(r'$\mu^2$')
save_plot("spin2_spectrum.pdf", output_dir="output")
```

## 使用例

### トーション構成でのLevi-Civita

```python
from dppu.utils import epsilon_symbol

# 軸性トーション T^{ijk} = (2η/r) ε^{ijk}
def build_axial_torsion(eta, r):
    T = {}
    for i in range(3):
        for j in range(3):
            for k in range(3):
                T[(i,j,k)] = (2*eta/r) * epsilon_symbol(i, j, k)
    return T
```

### Hodge dualでの4次元Levi-Civita

```python
from dppu.utils import levi_civita_4d

# (*R)^{ab}_{cd} = (1/2) ε_{cdef} R^{ab,ef}
def hodge_dual(R):
    R_star = np.zeros((4,4,4,4))
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    for e in range(4):
                        for f in range(4):
                            R_star[a,b,c,d] += 0.5 * levi_civita_4d(c,d,e,f) * R[a,b,e,f]
    return R_star
```

## 数学的背景

### Levi-Civita記号の定義

**3次元:**
```
ε_{ijk} = { +1  (i,j,k)が(0,1,2)の偶置換
          { -1  (i,j,k)が(0,1,2)の奇置換
          {  0  いずれかのインデックスが重複
```

**4次元:**
```
ε_{abcd} = { +1  (a,b,c,d)が(0,1,2,3)の偶置換
           { -1  (a,b,c,d)が(0,1,2,3)の奇置換
           {  0  いずれかのインデックスが重複
```

### 恒等式

```
ε_{ijk} ε_{ilm} = δ_{jl}δ_{km} - δ_{jm}δ_{kl}
ε_{abcd} ε_{abef} = 2(δ_{ce}δ_{df} - δ_{cf}δ_{de})
```

## 依存関係

- SymPy（シンボリック計算）
- NumPy（数値計算）
- matplotlib + seaborn（可視化）
- mpmath（高精度演算、オプション）

## 関連モジュール

- [curvature](../curvature/README_ja.md): Hodge dual計算
- [torsion](../torsion/README_ja.md): トーション構成
