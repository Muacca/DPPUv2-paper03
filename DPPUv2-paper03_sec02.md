# 2. モデルと既約な幾何変数

## coframeと作用

空間添字を $i,j,k=0,1,2$ 、円方向を $3$ 、向きを $\epsilon^{0123}=+1$ とする。円座標の周期は1、左不変形式は $d\sigma^i=2\epsilon_{ijk}\sigma^j\wedge\sigma^k$ を満たし、空間添字対は総和する。 $r,L>0$ 、実行列 $\det M=1$ に対して次を定義する。

式(M1)：

$$
\begin{aligned}
e^i&=rM^i{}_j\sigma^j,& e^3&=L(dy+\omega_i\sigma^i),\\
e'^a&=O^a{}_b e^b,\quad O\in SO(4),&
g&=\delta_{ab}e'^a e'^b,\qquad
\frac{d\mathcal V_4}{dy}=\frac{\pi^2Lr^3}{4}.
\end{aligned}
$$

式(M1)の空間積分済み四次元体積密度は、r,Lが一定なら $\pi^2Lr^3/4$ に積分される。round球面の幾何学的半径は $r/2$ である。三つのmixing変数を $\delta_i=\tan\chi_i$ 、 $O=O_{03}(\chi_0)O_{13}(\chi_1)O_{23}(\chi_2)$ と定義する。Oはmetricを変えないcoframe回転である。幾何学的torsionは回転後のframeで次のように指定する。

式(M2)：

$$
\begin{aligned}
T_{ijk}&=\frac{2\eta}{r}\epsilon_{ijk},\\
T^{(V)}_{abc}&=\frac13(\delta_{ac}V_b-\delta_{ab}V_c),
\qquad V_b=V\delta_{b3},\\
\mathcal S&=-\int\mathrm{vol}_4
\left(\frac{R}{2\kappa^2}+\theta N\right),\qquad
N\thinspace \mathrm{vol}_4=d(e'^a\wedge T_a).
\end{aligned}
$$

axial項は添字に3を含むとき零とする。AXは $V=0$ 、VTは $\eta=0$ 、MXは両者を許す。 $\theta$ は定数とし、曲率二乗項は加えない。計量整合接続はcoframeとこの幾何学的torsionから定まる。成分式は付録A、Paper01/02の付属コードとの符号対応は付録Iに示す。回転後のframeでtorsion成分を指定する処方と、全場を同時に回すgauge変換を区別する。

本論の符号はPaper01/02の縮約ポテンシャル $U=-S$ の規約を円方向作用へ拡張する。共通のECと完全NYのモデルでは $\mathcal S=-S_{\mathrm{01,FULL}}=-S_{\mathrm{02}}|_{\alpha=0}$ であり、一定profileでは $\mathcal S=U$ となる。以下の下限に関する主張は、この符号の $U$ または $\mathcal S$ について述べる。

## 共通の静的ポテンシャル

$M,\omega,O,r,L,\eta,V$ が一定の場合、 $S=MM^{\mathsf T}$ と置く。Cartan縮約により静的作用Uを得る。shape依存性をFにまとめると、

式(M3)：

$$
\begin{aligned}
U&=\mathcal K\left[F(S)+\frac{3\eta^2}{4}
+\frac{V^2r^2}{12}+\frac{L^2}{r^2}\lVert M\omega\rVert^2\right],\\
F(S)&=2\mathrm{tr}(S^2)-[\mathrm{tr}(S)]^2,\qquad
\mathcal K=\frac{\pi^2Lr}{\kappa^2}>0.
\end{aligned}
$$

完全NYはこの一定場族では点ごとに零であり、第6節で証明する。torsion normは回転不変で、その発散項は零となるため、scalar作用はmixing回転に依存しない。coframeは $S^3\times S^1$ 上で大域的に定義されるので、滑らかな周期的円方向変分でも完全NYの積分は零となる。その場合の局所密度は非零でもよい。

## 表現と変分族

スカラー変数は $(r,L,\eta,V)$ である。twist $\omega_i$ と三つのmixing座標はround frame近傍で空間ベクトルとして変換する。体積保存shapeは対称tracelessなhによる $M=e^h$ で表せ、正規化した五成分がspin-2をなす。軸対称参照shapeでは $M=M_0e^h$ を使う。基底は付録Cに示す。round点の静的二次展開では異なる既約表現が回転対称性によって分離する一方、スカラー同士の混合は許される。

静的極値は式(M3)について、制御変数を固定し、指定した座標を変化させた結果である。円方向作用汎関数では、 $r,L>0$ 、 $\det M=1$ のもとで全座標を滑らかな周期関数とする。スカラー制限では $M=I,\omega=0$ 、ベクトル制限では $r,L,\eta,V$ を固定して $M=I$ とする。いずれも以下の一般作用から得られる。停留・最小化の性質を述べる際は、固定する変数を明示する。

## 同時円方向変分の一般作用

$J=Lr^3$ と置く。空間の不変基底での誘導三計量・lapse・shiftは

式(M4)：

$$
\gamma=r^2M^{\mathsf T}M+L^2\omega\omega^{\mathsf T},
\qquad \mathcal N^2=\frac{J^2}{\det\gamma},
\qquad \beta=L^2\gamma^{-1}\omega.
$$

$[\beta]_\times v=\beta\times v$ と定義する。空間構造定数の規格化から、

式(M5)：

$$
\begin{aligned}
\Omega_\beta&=4[\beta]_\times,&
\ell=\mathcal L_\beta\gamma&=\gamma\Omega_\beta-\Omega_\beta\gamma,\\
X&=\gamma^{-1}\gamma',&
Y&=\gamma^{-1}\ell,\qquad \mathrm{tr}Y=0.
\end{aligned}
$$

primeは $d/dy$ を表す。非等方shapeではLie微分 $\ell$ は一般に零でない。全微分を円上で積分すると、exact作用は

式(M6)：

$$
\begin{aligned}
\mathcal S_{\mathrm{joint}}
&=\int_0^1dy\thinspace [U+\mathcal L_{\partial}],\\
\mathcal L_{\partial}
&=\frac{\pi^2J}{32\kappa^2\mathcal N^2}
\left[\mathrm{tr}(X^2)-(\mathrm{tr}X)^2-2\mathrm{tr}(XY)\right].
\end{aligned}
$$

Uは各yでの座標を代入した静的式(M3)である。shiftに依存する静的寄与は既にUに含まれるため、式(M6)へ $Y^2$ 項をさらに加えない。最初の二つのtrace項が速度二次形式、最後が速度に一次の項である。torsion微分とmixing回転は付録Jに示す全微分にまとまり、一定係数の完全NYは周期bulk作用に残らない。式(M6)は空間的に一様なこのEuclidean族の有限振幅で成立し、四次までの打切りではない。

[前の節・付録](DPPUv2-paper03_sec01.md) / [次の節・付録](DPPUv2-paper03_sec03.md) / [第1節](DPPUv2-paper03_sec01.md)
