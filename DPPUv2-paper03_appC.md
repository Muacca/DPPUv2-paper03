## Appendix C. Full Shear Tensor Basis and Embedding of $(\varepsilon,s,\phi)$

本付録では、spin-2 sector の最終基底と、raw 幾何変数の $(\varepsilon,s,\phi)$ 記述との関係をまとめる。

### C.1 Traceless symmetric basis

full shear tensor は

$$
h_{ij}=\sum_{A=1}^5 q_A (T_A)_{ij}
$$

と展開する。基底 $T_1,\dots,T_5$ の明示的な行列表現は §2.3 に掲載する。ここでは以下の性質が exact に成立することを確認する。

$$
{\rm Tr}(T_AT_B)=\delta_{AB},
\qquad
{\rm tr}(h)=0,
\qquad
{\rm tr}(h^2)=\sum_{A=1}^5 q_A^2.
$$

### C.2 Embedding of the old variables

raw 幾何変数 $(\varepsilon,s,\phi)$ は、spin-2 quintet の中の特定 3 次元 slice を与える記述として理解される。本論ではこの slice を主役にせず、最終基底 $q_A$ を一次的な変数として採用する。したがって、

- $\varepsilon$ は主として $m=0$ singlet 方向の軸対称変形
- $s$ と $\phi$ は負側 condensed phase を追うための補助パラメータ

として用い、最終的な representation-theoretic statement は常に $q_A$ で書き直す。

### C.3 Invariants in the final basis

最終基底での低次不変量は

$$
I_2={\rm tr}(h^2)=\sum_{A=1}^5 q_A^2,
\qquad
I_3=h_{ij}\omega_i\omega_j,
\qquad
I_4={\rm tr}(h^3)
$$

で整理されることが exact に確認された。

### C.4 Constant-shape homogeneous consistency

一定 frame 変形

$$
e^i=r\,M^i{}_j\,\sigma^j
$$

に対し、構造定数が

$$
C'^i{}_{jk}=
M^i{}_a (M^{-1})^b{}_j (M^{-1})^c{}_k C^a{}_{bc}
$$

と変換することを用いて homogeneous consistency を検証した。

all-five-nonzero の representative full-shear sample を 2 例とり、 $S^3$ と $Nil^3$ の両方について

- antisymmetry residual = 0
- Jacobi residual = 0

が exact に成立した。したがって final spatial ansatz は homogeneous minisuperspace frame として self-consistent である。

### C.5 Role of this appendix

この付録の役割は、raw 幾何変数から最終基底への視点の切り替えを明確にすることである。本論では、spin-2 についての最終主張は $q_A$ とその既約分解で書き、 $(\varepsilon,s,\phi)$ は相図説明の補助変数として用いる。

### C.6 Jacobian: old variables to final basis

raw 幾何変数 $(\varepsilon,s)$ と最終基底座標 $(q_1,q_2)$ の対応関係を明示する。

$\varepsilon$ は $T_1$ 方向の変形を指数的に埋め込んだもので、

$$
\varepsilon=\exp\!\left(\frac{\sqrt{6}}{2}\,q_1\right)-1,
\qquad
q_1=\frac{2}{\sqrt{6}}\ln(1+\varepsilon).
$$

等方点（$\varepsilon=0$）近傍での Jacobian は

$$
\frac{\partial q_1}{\partial\varepsilon}\bigg|_{\varepsilon=0}=\frac{2}{\sqrt{6}}=\sqrt{\frac{2}{3}}.
$$

同様に $s$ は $T_2$ 方向への変形で

$$
s=\exp\!\left(\frac{q_2}{\sqrt{2}}\right)-1,
\qquad
q_2=\sqrt{2}\ln(1+s),
\qquad
\frac{\partial q_2}{\partial s}\bigg|_{s=0}=\sqrt{2}.
$$

$s\neq0$ 相で現れる角度 $\phi$ は $(q_2,q_3)$ 平面内の SO(2) 回転として

$$
q_2=s_q\cos\phi,\qquad q_3=s_q\sin\phi,\qquad s_q=\sqrt{q_2^2+q_3^2}
$$

と対応する。 $T_4,T_5$ 方向（off-diagonal 成分）への同様の埋め込みは対称性から同形である。

等方点での合成 Jacobian の対角成分は

$$
J=\mathrm{diag}\!\left(\sqrt{\tfrac{2}{3}},\,\sqrt{2}\right)
$$

であり、Hessian の基底変換は

$$
H_q=J^{-T}\,H_{\varepsilon s}\,J^{-1}
$$

で与えられる。等方点では $H_{\varepsilon s}=\mu_q^2\,I_2$ かつ $J^{-T}J^{-1}\propto I_2$ であるため、変換後も Schur 縮退 $H_q=\mu_q^2\,\delta_{AB}$ が保たれる。
