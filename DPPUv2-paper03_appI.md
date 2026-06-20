## Appendix I. Reduction-order Artifact in the Spin-1 Sector

本付録では、rank-deficient field-space metric を持つ sector で、auxiliary 消去より先に補空間へ射影すると finite generalized eigenvalue ではない量が physical eigenvalue のように現れる具体例を整理する。目的は、spin-1 sector の Schur reduction を補足し、Appendix E の計算でなぜ auxiliary coupling を保持する必要があるかを明確にすることである。

### I.1 Projected quotient before auxiliary elimination

spin-1 の各軸では、raw variables $(\omega_k,\delta_k)$ に対する witness 行列

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

を用いる。kinetic metric には null direction

$$
n=\left(1,\frac23\right)^T
$$

があり、これに対する単純な補空間方向として

$$
u=(-2,3)^T
$$

を取ることができる。この方向に沿う Rayleigh quotient は

$$
\lambda_{\rm proj}
=\frac{u^TH_2u}{u^TG_2u}
=\frac{41899.77}{10007.77}
\simeq4.1867>0
$$

である。この quotient は auxiliary 消去前の projected coefficient の比であり、rank-deficient pencil の finite generalized eigenvalue そのものではない。

### I.2 Test against the generalized eigenvalue equation

$G_2$ は rank deficient であり、

$$
G_2n\simeq0
$$

を満たす。このとき、単純な射影 Rayleigh quotient は退化ペンシル

$$
\det(H_2-\lambda G_2)=0
$$

の有限一般化固有値と一致するとは限らない。実際、projection value を代入すると

$$
\det(H_2-\lambda_{\rm proj}G_2)\simeq -2.2\times10^6\neq0
$$

であり、さらに

$$
(H_2-\lambda_{\rm proj}G_2)u\neq0
$$

である。したがって $u=(-2,3)$ は一般化固有ベクトルではなく、 $\lambda_{\rm proj}=4.1867$ は finite generalized eigenvalue ではない。

### I.3 Missing auxiliary coupling

projection quotient と Schur-reduced eigenvalue が異なる原因は、non-dynamical auxiliary $Y_k$ との Hessian cross term にある。実際、

$$
H_{YY}=n^TH_2n\simeq894.845>0,
\qquad
H_{XY}=n^TH_2u
=\frac{24081847}{3750}
\simeq6421.8\neq0.
$$

したがって $n$ 方向は gauge/redundancy ではなく second-class auxiliary であり、physical mode を得る前に代数 EOM によって消去しなければならない。 $H_{XY}$ を保持しない reduction では、auxiliary の backreaction が欠落する。

### I.4 Comparison with the Schur-reduced result

正しい reduction order は次の通りである。

1. rank-deficient kinetic metric の null direction $Y_k$ を同定する。
2. $Y_k$ が second-class auxiliary であることを確認する。
3. $X$ - $Y$ cross term を含めた quadratic form から、 $Y$ の代数 EOM を解く。
4. Schur 補

$$
H_{\rm eff}=H_{XX}-\frac{H_{XY}^2}{H_{YY}}
$$

を用いて physical finite eigenvalue を評価する。

この手順に従うと

$$
\lambda_{\rm phys}
=\frac{H_{XX}-H_{XY}^2/H_{YY}}{G_{XX}}
\simeq -0.4183<0
$$

となる。symbolic Ricci 正規化でも同じ Schur-reduced structure は

$$
\lambda_{\rm phys}^{\rm symbolic}
=-\frac{56L^2}{3r^2}<0
$$

を与える。したがって、auxiliary direction を Schur 消去しない場合の positive projection value と、Schur 消去後の finite generalized eigenvalue は符号も物理分類も異なる。正しい分類は、kinetic norm 正の non-ghost tachyon である。

### I.5 Basis-independence lesson

補空間方向 $u=(-2,3)$ を使うこと自体が本質的な問題ではない。問題は、rank-deficient metric を持つ sector で auxiliary direction を Schur 消去せずに、消去前の projected quotient を physical generalized eigenvalue と同一視することである。実際、同じ補空間方向から出発しても、 $Y$ を正しく消去すれば finite eigenvalue は負になる。

したがって本例から得られる実用的な教訓は次の通りである。

> rank-deficient field-space metric を持つ sector では、naive 直交射影の Rayleigh quotient を physical generalized eigenvalue と同一視してはならない。second-class auxiliary が存在する場合、Schur 消去または退化ペンシルの finite generalized eigenvalue を用いる必要がある。

### I.6 Relation to reproducibility checks

本稿v1 ではSchur 消去前の projected quotient を physical generalized eigenvalue と同一視していたが、v2 ではSchur 消去後の finite generalized eigenvalue を採用した。そのため、以下の本文・付録について訂正した。

| verification target | section / appendix |
|---|---|
| v1 で採用していた $4.1867>0$ は projected quotient であり、finite generalized eigenvalue ではないこと。 | Appendix I |
| Schur 消去後の有限一般化固有値が witness 正規化で $\lambda_{\rm phys}\simeq -0.4183<0$、symbolic 正規化で $-56L^2/(3r^2)<0$ となること。 | §5.2, Appendix E |
| corrected spin-1 mode が kinetic norm 正の非ゴースト・タキオンであること。 | §5.2, §5.5, Appendix E |
| spin-1 不安定が curvature-scale の tachyonic non-ghost mode であり、frame-degeneracy endpoint の発散や radial decompactification とは別系統であること。 | §5.5, §7 |

### I.7 Computational verification

本付録の数値は以下のスクリプトで確認されている。

- **Schur 消去前の projected quotient $\lambda_{\rm proj}\simeq4.1867>0$ を再現し、それが finite generalized eigenvalue ではないことを確認する**:
  - `paper03/spin1_appI_projection_artifact.py`

