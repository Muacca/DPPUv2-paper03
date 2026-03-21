## Appendix E. Spin-1 Sector: $H$, $G$, and Generalized Eigenvalue Analysis

本付録では、spin-1 sector の raw 二次形式と final basis への変換をまとめる。

### E.1 Raw quadratic blocks

等方点 $r_0=3,\ L=1$ では

$$
H_2=
\begin{pmatrix}
52.6379 & -631.6544\\
-631.6544 & 3789.9295
\end{pmatrix},
\qquad
G_2=
\begin{pmatrix}
236.870435 & -355.305547\\
-355.305547 & 532.958162
\end{pmatrix}.
$$

raw Hessian のみでは physical instability を判定できず、 $(H,G)$ の組として扱う必要がある。

### E.2 Null mode and physical mode

null direction は

$$
Y_k\propto(\omega_k,\delta_k)=\left(1,\frac23\right),
\qquad
k=0,1,2
$$

であり、physical direction は

$$
X_k=\frac{-2\omega_k+3\delta_k}{\sqrt{13}}
$$

で与えられる。

この基底では

$$
G\,Y_k=0,
\qquad
G_{\rm phys}=769.8286\,I_3,
\qquad
H_{\rm phys}=3223.06\,I_3
$$

となり、

$$
\lambda_{\rm phys}=4.1867>0
$$

が従う。

### E.3 Auxiliary interpretation summary

$Y_k$ は gauge でも variable redundancy でもない。実際、

$$
\delta g_{\tau,\sigma^0}\neq0,
\qquad
H_{\rm null}=894.845>0,
\qquad
\lambda_{\rm eff}\to\infty
$$

であるため、 $Y_k$ は kinetic term を持たない non-dynamical auxiliary と解釈される。variable redundancy 排除の詳細は E.4 にまとめる。

### E.4 variable redundancy 排除の詳細

§5.3 での $Y_k$ の variable redundancy 排除を詳述する。

**variable redundancy の必要条件**は次の 3 つである。

1. $\delta g_{\mu\nu}=0$（metric 変分がゼロ）
2. $H_{\rm null}=0$（auxiliary 方向 Hessian がゼロ）
3. 有効固有値 $\lambda_{\rm eff}$ が有限

以下、各条件が満たされないことを示す。

**条件 1 の排除（ $\delta g_{\tau,\sigma^0}\neq0$ ）**

軸 $k=0$ に対して、 $Y_0$ 方向の変分は twist 成分と mixing 成分それぞれ

$$
\delta e^3_{\rm twist}=+\frac{2L}{r_0}e^0,
\qquad
\delta e^3_{\rm mixing}=-\frac{2L}{r_0}e^0
$$

を与え、両者の合計では $\delta e^3=0$ となる。しかし coframe $e^0$ の変分は

$$
\delta e^0=+\frac{2L}{r_0}e^3
$$

として残り、これにより metric 変分

$$
\delta g_{\tau,\sigma^0}=(\delta e^0)_\tau\,(e^0)_{\sigma^0}
=\frac{2L}{r_0}\,(e^3)_\tau\,(e^0)_{\sigma^0}
\neq0
$$

が生じる。 $r_0=3,\ L=1$ では $2L/r_0=2/3$ であり、metric 変分は非ゼロである。したがって必要条件 1 は排除される。 $k=1,2$ 軸でも同様の結論が成り立つ。

**条件 2 の排除（ $H_{\rm null}>0$ ）**

null direction $Y_k=(1, 2/3)^T$（各軸成分）に沿う Hessian 成分は

$$
Y=\begin{pmatrix}1\\
2/3\end{pmatrix},\qquad
H_2=
\begin{pmatrix}
52.6379 & -631.6544\\
-631.6544 & 3789.9295
\end{pmatrix}
$$

$$
H_{\rm null}=Y^T H_2 Y=894.845>0
$$

である。したがって必要条件 2 も排除される。

**条件 3 の排除（ $\lambda_{\rm eff}\to\infty$ ）**

一般化固有値の Rayleigh quotient は

$$
\lambda_{\rm eff}(v)=\frac{v^T H_2\, v}{v^T G_2\, v}.
$$

$v=Y_k$ のとき $G_2\,Y_k=0$ であるから分母がゼロとなる。一方 $H_{\rm null}=Y^T H_2 Y=894.845>0$ であるから分子は正に留まり、

$$
\lambda_{\rm eff}(Y_k)\to+\infty.
$$

したがって必要条件 3 も排除される。

以上より $Y_k$ は variable redundancy ではなく、**non-dynamical auxiliary** であることが確定する。

### E.5 Final statement

spin-1 sector の最終結論は次の 2 点である。

1. propagating triplet は $X_k$ のみである。
2. $Y_k$ は constraint-like auxiliary direction であり、physical mass spectrum に含めない。

したがって、raw basis の負固有値をそのまま instability と読むことはできない。

### E.6 Computational verification

本付録の結果は次のスクリプトで確認されている。

- **$H_2$, $G_2$ の数値計算と physical/null 基底**:
  - `paper03/spin1_null_mode_analysis.py`
- **$G$ metric（field-space metric）の strict 計算**:
  - `paper03/spin1_kinetic_metric.py`
- **variable redundancy 排除の数値確認**（ $\delta g\neq0$, $H_{\rm null}>0$, $\lambda_{\rm eff}\to\infty$ ）:
  - `paper03/spin1_auxiliary_exclusion.py`
