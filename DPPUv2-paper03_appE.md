## Appendix E. Spin-1 Sector: Schur Reduction and Finite Generalized Eigenvalue

本付録では、spin-1 sector の raw 二次形式、rank-deficient kinetic metric、non-dynamical auxiliary の同定、および Schur 消去後の finite generalized eigenvalue をまとめる。

### E.1 Raw quadratic blocks

等方点 $r_0=3,\ L=1$ では、各軸 $k=0,1,2$ ごとに

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
\end{pmatrix}
$$

を得る。ここで $H_2$ は potential Hessian、 $G_2$ は field-space kinetic metric である。raw Hessian 単独の固有値は physical mass 判定には不十分であり、rank-deficient kinetic metric と auxiliary constraint を含む二次形式として扱う。

### E.2 Null direction and auxiliary character

kinetic metric には null direction

$$
Y_k\propto(\omega_k,\delta_k)=\left(1,\frac23\right)
$$

が存在し、

$$
G_2Y_k\simeq0,
\qquad
H_{YY}=Y_k^TH_2Y_k\simeq894.845>0
$$

を満たす。したがって $Y_k$ は kinetic term を持たないが potential curvature は持つ second-class auxiliary である。

$Y_k$ は gauge でも variable redundancy でもない。実際、軸 $k=0$ に対して $Y_0$ 方向の変分は

$$
\delta e^3_{\rm twist}=+\frac{2L}{r_0}e^0,
\qquad
\delta e^3_{\rm mixing}=-\frac{2L}{r_0}e^0
$$

を与え、両者の合計では $\delta e^3=0$ となる。しかし coframe $e^0$ の変分は

$$
\delta e^0=+\frac{2L}{r_0}e^3
$$

として残るため、

$$
\delta g_{\tau,\sigma^0}\neq0
$$

が生じる。さらに $H_{YY}>0$ であり、 $G_2Y_k=0$ のため $Y_k$ 方向の有効 Rayleigh quotient は $\lambda_{\rm eff}\to\infty$ となる。よって $Y_k$ は variable redundancy ではなく、non-dynamical auxiliary と解釈される。

### E.3 Non-null coordinate and Schur data

null direction に対する補空間座標として

$$
X_k=\frac{-2\omega_k+3\delta_k}{\sqrt{13}}
$$

を用いる。この $X_k$ は kinetic metric の非零方向を張る便利な座標であり、 $X$-block の係数は

$$
G_{XX}\simeq769.8286,
\qquad
H_{XX}\simeq3223.06
$$

である。ただし $Y_k$ は auxiliary として同じ二次形式内に残るため、finite physical eigenvalue は $X$-block だけからではなく、 $X$-$Y$ system の Schur 消去から定まる。

必要な残りの係数は

$$
H_{YY}\simeq894.845,
\qquad
H_{XY}
=\frac{24081847}{3750}
\simeq6421.8\neq0
$$

である。非ゼロの $H_{XY}$ は、kinetic-null 方向を消去する前に落としてはならない auxiliary coupling である。

### E.4 Schur elimination

spin-1 sector の二次 Lagrangian は、各軸について

$$
\mathcal{L}_{X,Y}
=\frac12\dot X^TG_{XX}\dot X
-\frac12X^TH_{XX}X
-X^TH_{XY}Y
-\frac12Y^TH_{YY}Y
$$

と書ける。 $Y$ は kinetic term を持たないため、代数 EOM

$$
\frac{\partial V}{\partial Y}=0
\qquad\Rightarrow\qquad
Y=-\frac{H_{XY}}{H_{YY}}X
$$

で消去される。これにより physical Hessian は

$$
H_{\rm eff}=H_{XX}-\frac{H_{XY}^2}{H_{YY}}
$$

となり、finite generalized eigenvalue は

$$
\lambda_{\rm phys}
=\frac{H_{\rm eff}}{G_{XX}}
=\frac{H_{XX}-H_{XY}^2/H_{YY}}{G_{XX}}
\simeq -0.4183<0
$$

である。symbolic Ricci 正規化では同じ Schur 構造が

$$
\lambda_{\rm phys}^{\rm symbolic}
=-\frac{56L^2}{3r^2}<0
$$

として得られる。両者の magnitude は正規化が異なるため直接一致させないが、符号と物理分類は一致する。

### E.5 Corrected mode and ghost test

finite generalized eigenvalue に対応する corrected mode は、Euclidean unit normalization で概ね

$$
v_{\rm corr}\simeq(0.9816,\ 0.1909)
$$

である。この mode の kinetic norm は

$$
v_{\rm corr}^TG_2v_{\rm corr}\simeq114.5>0
$$

である。したがって corrected spin-1 mode は ghost ではなく、負の mass squared を持つ tachyonic mode である。

### E.6 Final statement

spin-1 sector の結論は次の通りである。

1. $Y_k$ は gauge でも redundancy でもなく、non-dynamical auxiliary である。
2. $X$ - $Y$ cross term $H_{XY}$ は非ゼロであり、finite physical eigenvalue の評価には Schur 消去が必要である。
3. $Y_k$ を Schur 消去した finite generalized eigenvalue は $\lambda_{\rm phys}\simeq-0.4183<0$ である。
4. corrected physical mode は kinetic norm 正の非ゴースト・タキオンである。

### E.7 Computational verification

本付録の計算では、再現性スクリプトにより次の値を確認する。

- $\lambda_{\rm phys}\simeq-0.4183$
- $\lambda_{\rm phys}^{\rm symbolic}=-56L^2/(3r^2)$
- $H_{XY}=24081847/3750$
- $v_{\rm corr}\simeq(0.9816,0.1909)$
- $v_{\rm corr}^TGv_{\rm corr}\simeq114.5>0$

#### スクリプト

- **第一原理 symbolic 導出: `dppu.action.spin1_gradient` エンジンを用いて $G$, $H$ を記号で導出し、 $\lambda_{\rm phys}^{\rm symbolic}=-56L^2/(3r^2)$ を式から導き、退化ペンシルの有限根との一致を確認、さらに $G$ を volume factor $\pi^2 L r^3$ でスケールして witness 値と突合する**:
  - `paper03/spin1_symbolic_firstprinciples.py`
- **raw $G,H$ witness から $G_{XX},H_{XX},H_{XY},H_{YY}$ を組み立て、Schur 消去後の $\lambda_{\rm phys}\simeq-0.4183<0$ と kinetic norm 正を確認する（symbolic 導出も内包）**:
  - `paper03/spin1_kinetic_metric.py`
- **$Y_k$ が kinetic-null かつ non-dynamical auxiliary であることを確認し、Schur 消去後の spin-1 classification と接続する**:
  - `paper03/spin1_null_mode_analysis.py`
- **$Y_k$ が gauge / variable redundancy ではないことを、coframe 変分、 $H_{YY}>0$、 $\lambda_{\rm eff}\to\infty$ で確認する**:
  - `paper03/spin1_auxiliary_exclusion.py`

