## Appendix B. Detailed Proof of the Pontryagin Protection Theorem

本付録では、§3 の Pontryagin 保護定理を代数的に導出する。

### B.1 Statement

主定理は次の通りである。

> **Theorem P0.**
> AX mode および VT mode では、homogeneous minisuperspace の対象トポロジーに対して
$P=\langle R,\star R\rangle\equiv0$
が成立する。

### B.2 AX sector: algebraic proof of Plucker/Pfaffian cancellation

#### B.2.1 Setup: EC connection for S³×S¹ AX mode

##### コフレームの定義

本節全体で用いる $S^3\times S^1$ コフレームを次のように定める。

$$
e^i = r\,M^i{}_j\,\sigma^j \quad (i=0,1,2),
\qquad
e^3 = L\!\left(d\tau + \omega_k\,\sigma^k\right)
$$

ここで $\sigma^i$ は $S^3$ の左不変 1-形式、 $M^i{}_j=(\exp h)^i{}_j$（ $h=\sum_A q_A T_A$: shear tensor）、 $\omega_k$（ $k=0,1,2$）は vector twist 成分である。等方 AX mode では $h=0$（ $M=I$ ）、 $\omega_k=\omega$（全成分共通の scalar）とする。spin-1 mixing $\delta_k$（spatial leg $e^k$ への $d\tau$ 成分付加）は独立な自由度であり、App E で扱う。

---

以下、空間指標 $i,j,k\in\{0,1,2\}$、fiber 指標 $a=3$（ $\tau$ 方向）とする。等方 twist $\omega$ に対する構造定数は Maurer-Cartan 方程式

$$
de^k = \frac{1}{r}\varepsilon_{ijk}\,e^i\wedge e^j,
\qquad
de^3 = \frac{2L\omega}{r^2}\varepsilon_{ij}\,e^i\wedge e^j
$$

から

$$
C^k_{ij}=\frac{2}{r}\varepsilon_{ijk},
\qquad
C^3_{ij}=\frac{4L\omega}{r^2}\varepsilon_{ij}
$$

と決まる。ここで $\varepsilon_{ij}$ は $\{0,1,2\}$ 上の巡回符号（ $\varepsilon_{01}=\varepsilon_{12}=\varepsilon_{20}=1$ ）である。 
なお、 $\tau$ 方向は空間 $S^3$ に構造定数を生成しないことから、 $C^i_{3j}=0$ となる。

**Levi-Civita 接続** 
( Koszul 公式 $\Gamma^{a}_{bc}=\tfrac12(C^{a}\_{bc}+C^{c}\_{ba}-C^{b}\_{ac})$ ):

$$
\Gamma^k_{ij}=\frac{1}{r}\varepsilon_{ijk},
\qquad
\Gamma^3_{ij}=\frac{2L\omega}{r^2}\varepsilon_{ij},
\qquad
\Gamma^i_{3j}=-\frac{2L\omega}{r^2}\varepsilon_{ij}.
$$

これは

$$
\Gamma^2_{01}=\tfrac12(C^2_{01}+C^1_{20}-C^0_{12})=\tfrac{1}{r}(\varepsilon_{201}+\varepsilon_{120}-\varepsilon_{012})=\tfrac{1}{r}(1+1-1)=\tfrac{1}{r}
$$

として確認できる。

**AX torsion の contorsion**:
$T_{ijk}=\tfrac{2\eta}{r}\varepsilon_{ijk}$（空間成分のみ非零）に対し contorsion 公式

$$
K_{abc}=\frac{1}{2}(T_{abc}+T_{bca}-T_{cab})
$$

を適用すると、完全反対称テンソルでは $K_{ijk}=\tfrac{\eta}{r}\varepsilon_{ijk}$、他はゼロ。

**EC 接続**
$\omega^a_{bc}=\Gamma^a_{bc}+K^a_{bc}$: 二つのパラメータを定義する。

$$
c \equiv \frac{1+\eta}{r},
\qquad
f \equiv \frac{2L\omega}{r^2}.
$$

非零成分:

$$
\omega^k_{ij}=c\,\varepsilon_{ijk},
\qquad
\omega^3_{ij}=f\,\varepsilon_{ij},
\qquad
\omega^i_{3j}=-f\,\varepsilon_{ij},
\qquad
\omega^i_{j3}=+f\,\varepsilon_{ij}.
$$

#### B.2.2 Riemann tensor for homogeneous spaces

等質空間では frame 微分がゼロになるため、Riemann テンソルは代数式

$$
R^a{}_{bcd}=\omega^a{}_{ec}\,\omega^e{}_{bd}-\omega^a{}_{ed}\,\omega^e{}_{bc}+\omega^a{}_{be}\,C^e{}_{cd}
$$

で与えられる（構造定数項が偏微分項を置き換える）。

#### B.2.3 Representative Riemann components for (c,d) = (0,1)

$(c,d)=(0,1)$ SS-type ペアに対して Plücker 条件に必要な成分を順次計算する。

**$R_{0101}$** は

$$
R^0{}_{101}
= -\bigl(\omega^0{}_{21}\,\omega^2{}_{10}+\omega^0{}_{31}\,\omega^3{}_{10}\bigr)
+\omega^0{}_{12}\,C^2{}_{01}+\omega^0{}_{13}\,C^3{}_{01}.
$$

ここに以下を代入すると

$\omega^0_{21}=-c, \qquad \omega^2_{10}=-c, \qquad \omega^0_{31}=-f, \qquad \omega^3_{10}=-f \qquad \omega^0_{12}=c, \qquad \omega^0_{13}=f, \\ 
C^2_{01}=2/r, \qquad C^3_{01}=2f$

**$R_{0101}$** は以下のように整理される。

$$
R^0{}_{101}=-(c^2+f^2)+\frac{2c}{r}+2f^2
=-c^2+f^2+\frac{2c}{r}.
$$

さらに、 $c=(1+\eta)/r$ を代入し、以下を得る。

$$
\boxed{R_{0101}=\frac{1-\eta^2}{r^2}+f^2.}
$$

**$R_{2301}=R^2{}_{301}$** は

$$
R^2{}_{301}=\omega^2{}_{e0}\,\omega^e{}_{31}-\omega^2{}_{e1}\,\omega^e{}_{30}+\omega^2{}_{3e}\,C^e{}_{01}
$$

$\omega^{2}{}\_{e0}$ は空間指標の組で $\omega^{2}\_{i0}=c\varepsilon_{2i0}$ だが $\varepsilon_{2i0}$ は全て 0（ $i=0$ ）または $\varepsilon_{210}=-1$（ $i=1$ ）。
$\omega^2_{10}=-c$ と $\omega^1_{31}=0$（反対称）の組が積になるため第一項は $(-c)\cdot 0=0$ のみ。
第二項も $\omega^2_{01}=c$, $\omega^1_{30}=f$, $\omega^2_{31}=f$, $\omega^2_{30}=-f$ で計算するとゼロ。structure constant 項は $\omega^2_{32}=0$ により消える。

$$
\boxed{R_{2301}=0.}
$$

**$R_{0201}=R^0{}_{201}$** は

$$
R^0{}_{201}=0-\omega^0{}_{e1}\,\omega^e{}_{20}+\text{struct.}
=0-(-f^2)+(-2f^2)=-f^2.
$$

よって

$$
\boxed{R_{0201}=-f^2.}
$$

**$R_{1301}=R^1{}_{301}$** は

$$
R^1{}_{301}=cf-0-\frac{2f}{r}=f\!\left(c-\frac{2}{r}\right)=\frac{f(\eta-1)}{r}.
$$

よって

$$
\boxed{R_{1301}=\frac{f(\eta-1)}{r}.}
$$

**$R_{0301}=R^0{}_{301}$** も同様の計算により

$$
\boxed{R_{0301}=\frac{f(1-\eta)}{r}=-R_{1301}.}
$$

**$R_{1201}=R^1{}_{201}$** は

$$
R^1{}_{201}=-f^2+2f^2=f^2.
$$

よって

$$
\boxed{R_{1201}=f^2.}
$$

#### B.2.4 Verification of Pfaffian = 0 for (c,d) = (0,1)

$$
\mathrm{Pf}\,\omega^{(01)}
=R_{0101}\cdot R_{2301}-R_{0201}\cdot R_{1301}+R_{0301}\cdot R_{1201}
$$

$$
=\left(\frac{1-\eta^2}{r^2}+f^2\right)\cdot 0
-(-f^2)\cdot\frac{f(\eta-1)}{r}
+\frac{f(1-\eta)}{r}\cdot f^2
$$

$$
=0+\frac{f^3(\eta-1)}{r}+\frac{f^3(1-\eta)}{r}
=\frac{f^3}{r}\bigl[(\eta-1)+(1-\eta)\bigr]=0.
\qquad\checkmark
$$

#### B.2.5 Verification for (c,d) = (0,3) S3-type

$(c,d)=(0,3)$ に対して、 $C^e_{03}=0$（ $\tau$ 方向からの structure constant はゼロ）を利用すると:

$$
R_{0103}=-fc,\quad
R_{2303}=f^2,\quad
R_{0203}=-fc,\quad
R_{1303}=f^2,\quad
R_{0303}=-f^2,\quad
R_{1203}=0.
$$

したがって

$$
\mathrm{Pf}\,\omega^{(03)}
=(-fc)\cdot f^2 -(-fc)\cdot f^2+(-f^2)\cdot 0
=-f^3c+f^3c=0.\qquad\checkmark
$$

#### B.2.6 Generalization to all (c,d) pairs

$(c,d)$ ペアは SS 型（ $\{01,02,12\}$ ）と S3 型（ $\{03,13,23\}$ ）の 6 種類である。

- **SS 型**: $\{0,1,2\}$ の巡回置換の対称性（S³ の SO(3) 対称性）により $(0,1)$ と同じ構造を持つ。それぞれの Plücker 条件は B.2.4 と同一の代数的理由で成立する。
- **S3 型**: $(0,3)$ の計算で見たように、 $C^e_{i3}=0$（ $\tau$ 方向は構造定数を生成しない）が key となり、すべての S3 型ペアで同様に成立する。
- **P₃₃ 項**: $R_{3\cdots 3\cdots}$ は反対称性より自明にゼロ。

以上から、Pontryagin 密度は

$$
P=2\sum_{(c,d)}\mathrm{Pf}\,\omega^{(cd)}=0
$$

と代数的に確立される。

**注**: 上記の証明は等方 twist の場合のものである。一般ベクトル twist $(\omega_0,\omega_1,\omega_2)$ についても同じ構造定数の形（ $C^3_{ij}$ が各軸方向に独立に寄与する）が保たれ、 $V_{\rm eff}=f(I_2)$ の恒等的成立と合わせて $P\equiv0$ が成立することが確認されている（ $Nil^3 \times S^1$ を含む）。

### B.3 VT sector: Pfaffian cancellation

VT mode では消去機構が異なる。各 $(c,d)$ ブロックの Pfaffian は個別には非ゼロ（rank 4）であるが、SS 群と S3 群が符号反転で相殺する。 $S^3\times S^1$ での結果:

$$
\sum_{\rm SS}\mathrm{Pf}\,\omega^{(cd)}=-\frac{16L^2V\,I_2}{3r^5},
\qquad
\sum_{\rm S3}\mathrm{Pf}\,\omega^{(cd)}=+\frac{16L^2V\,I_2}{3r^5},
$$

ゆえに

$$
P=\sum_{\rm SS}\mathrm{Pf}+\sum_{\rm S3}\mathrm{Pf}=0.
$$

この相殺は「各ブロックが decomposable」という AX mode の機構とは異なる。VT mode では V=0 でなく $\eta=0$ が成立し、EC 接続の構造が変わるため各ブロックは rank 4 になる。しかしながら、SS 群の Pfaffian 総和と S3 群の Pfaffian 総和は

$$
\sum_{\rm SS}\mathrm{Pf}=\frac{16L^2V}{3r^5}\sum_{(i,j)\in\rm SS}R_{ij,ij},
\qquad
\sum_{\rm S3}\mathrm{Pf}=-\sum_{\rm SS}\mathrm{Pf}
$$

として等大逆符号になる。この符号反転は、 $V$ が pure vector-trace torsion として $\tau$ 方向にのみ入ることで、SS と S3 の dual 構造が相補的に働くことに由来する。Nil³×S¹ でも同様の相殺が成立することが代数的に確認されている。

### B.4 MX sector: the unique exception

MX mode では $\eta\neq0$ かつ $V\neq0$ が同時に成立し、保護は破れる。 $S^3\times S^1$ での Pontryagin 密度:

$$
P(\boldsymbol{\omega})=
\frac{2V\eta\,(V^2r^2+9\eta^2-36)}{9r^3}+
\frac{40V\eta L^2}{3r^5}\,I_2.
$$

等方点 ($I_2=0$) でも

$$
P_0=\frac{2V\eta\,(V^2r^2+9\eta^2-36)}{9r^3}\neq0
$$

が残り、 $P\neq0$ の起源は torsion mixing $V\eta$ であって twist $\omega$ 自身ではない。

### B.5 NyVariant independence

Riemann テンソル（したがって $P$）は NyVariant（FULL/TT/REE）の選択に依存しない。これはパイプライン計算順序（Riemann は E4.7、Nieh-Yan 密度は E4.10）から理論的に保証される。モード別の結果:

| mode | FULL | TT | REE |
|---|---|---|---|
| AX | 0 | 0 | 0 |
| VT | 0 | 0 | 0 |
| MX | nonzero | nonzero | nonzero |

### B.6 Interpretation

この付録の要点は次の通りである。

1. AX と VT は異なる代数機構で保護される。AX は各 $(c,d)$ ブロックの個別消去（Plücker 型）、VT は SS と S3 群 Pfaffian の符号相殺。
2. MX が唯一の例外であり、twist 自身は生成起源ではない。
3. Pontryagin 構造は、最終基底の整理と両立する形で sector ごとに明確に切り分けられる。

### B.7 Computational verification

本付録の代数的計算は以下のスクリプトで記号的に確認されている。

- **AX mode: Riemann ブロック分類・各 $(c,d)$ Pfaffian ゼロ確認**:
  - `paper03/pontryagin_anisotropic_twist.py`
- **VT / MX mode の NyVariant 横断確認（全 18 ケース）**:
  - `paper03/torsion_mode_analysis.py`
- **等方 twist $P\equiv0$ スキャン（ $r$, $\eta$ 依存性）**:
  - `paper03/pontryagin_scan_s3.py`
  - `paper03/pontryagin_scan_nil3.py`
