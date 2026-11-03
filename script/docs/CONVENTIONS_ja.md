# 幾何・計算規約

[English](CONVENTIONS.md) / [再現方法](../README_ja.md) / [実装ガイド](IMPLEMENTATION_ja.md)

定義は本文第2節と付録A・B・I・Jに示す。添付原稿は
[第2節](../../LaTeX/sections/sec2.tex)、[付録A](../../LaTeX/appendices/appA.tex)、
[付録B](../../LaTeX/appendices/appB.tex)、[付録I](../../LaTeX/appendices/appI.tex)、
[付録J](../../LaTeX/appendices/appJ.tex) にある。

## 1. モデルと大域幾何

空間的一様Euclidean EC、一定係数の完全NY、曲率二乗係数alpha=0を用いる。
空間添字は0,1,2、円方向は3、向きはepsilon0123=+1である。
円座標yの周期は1。一定Lで周長はL、変化するprofileでは固定空間点に沿う長さはL(y)の積分となる。
round球面の半径はr/2、空間積分済み四次元体積密度は $\pi^2Lr^3/4$ 。
Mは行列式1の実3×3行列である。

## 2. 成分配列と静的曲率

`C[a,b,c]`、`Gamma[a,b,c]`、`R[a,b,c,d]` は順に
$`C^a{}_{bc}`$ 、 $`\Gamma^a{}_{bc}`$ 、 $`R^a{}_{bcd}`$ を保存する。
接続の最後の添字は微分方向を表し、frame添字の上下にはEuclidean単位計量を用いる。

$$
\begin{aligned}
de^a&=\frac12 C^a{}_{bc}e^b\wedge e^c,\qquad
\Gamma^a{}_{b}=\Gamma^a{}_{bc}e^c,\\
T^a{}_{bc}&=C^a{}_{bc}+\Gamma^a{}_{cb}-\Gamma^a{}_{bc},\\
\Gamma^{\mathrm{LC}}_{abc}&=\frac12(C_{abc}+C_{cba}-C_{bac}),\\
K_{abc}[T]&=-\frac12(T_{abc}+T_{bca}-T_{cab}),\\
R^a{}_{bcd}&=\Gamma^a{}_{ec}\Gamma^e{}_{bd}
-\Gamma^a{}_{ed}\Gamma^e{}_{bc}+\Gamma^a{}_{be}C^e{}_{cd}.
\end{aligned}
$$

最後の式はframe成分が一定の場合に用いる。円方向に依存する接続の微分は、
別のbulk式(M6)で扱う。添字を下げたEC曲率は各対の中で反対称だが、対交換対称性は仮定しない。
一定coframe変換e=A ebarでは、構造定数の出力添字にA、二つの入力添字に実際の逆行列を使う：
`A[a,d] inv[e,b] inv[f,c] Cbar[d,e,f]`。

## 3. Torsion入力

本文と `paper_results` は幾何学的Tを用いる：
$`T_{ijk}=(2\eta/r)\varepsilon_{ijk}`$ 、
$`T^{(V)}_{abc}=(\delta_{ac}V_b-\delta_{ab}V_c)/3`$ 、
$`V_b=V\delta_{b3}`$ 。
AXはV=0、VTはeta=0、MXは両者を許し、いずれも零振幅を含む。
比較対象のPaper01/02付属コードでは入力配列がtau=-Tなので、`eta_code=-eta`、`V_code=-V` と両係数の符号を変える。
`cartan` は幾何学的Tを直接返す。`torsion_residual` はCを含めたCartan再構成の残差を表す。

## 4. NYと二つの曲率縮約

$$
\begin{aligned}
N&=TT-REE
=\frac14\varepsilon^{abcd}T^e{}_{ab}T_{ecd}
-\frac12\varepsilon^{abcd}R_{abcd},\\
P&=\frac14\varepsilon^{cdef}R_{abcd}R_{abef},\\
B&=\sum_{c<d}(R_{01cd}R_{23cd}-R_{02cd}R_{13cd}+R_{03cd}R_{12cd}).
\end{aligned}
$$

Nはスカラー密度で、体積形式を掛けたものがNY四形式となる。
標準のform添字Pontryagin密度Pと内部縮約Bは異なる。全成分Hodge内積は2Pである。
`cartan` は上記の規格化で `TT, REE, NY, P, B` を返す。
一定unimodular族ではN=P=0。滑らかな周期profileでは局所密度が非零でも完全NYの積分は零となる。

## 5. 作用の符号と変分の定義域

本論では $\mathcal S=-\int\mathrm{vol}_4(R/(2\kappa^2)+\theta N)$ とする。
一定profileでは $\mathcal S=U$ であり、Paper01/02のU=-Sに合わせた符号である。
両論文の正の積分で定義する作用Sとは全体符号が逆になる。

スカラー固定でshapeとtwistを静的に最小化するとround・twistゼロが選ばれる。
非零twistを固定すると、静的な最小shapeは移動し得る。
滑らかな周期profileの変分は別の問題であり、周波数制限がないと連成shape–twist系列は下に有界でない。
`joint_bulk` は全微分を除いた周期bulk密度を返し、非等方shiftと速度一次項を含む。
入力は `tr(inv(M) @ M_prime)=0` を満たす。

## 6. Shape経路と制限関数

shear基底は $`\mathrm{tr}(T_AT_B)=\delta_{AB}`$ で正規化する。
`shear_curvatures` の既定経路は `M=M0 exp(q T_A)`。
`linear_normalized` は別の経路を指定する。非停留点の二階微分は経路に依存する。
静的・twistゼロでのPaper02との座標対応は $x=a^{-2}=(1+\epsilon)^{-2/3}$ である。

`scalar_bulk` はround・twistゼロで周期スカラーprofileを扱う。
`twist_bulk` はround・スカラー固定である。四次項は振幅の四次で、円方向微分については二次である。
掲載Hessianは、これら指定された幾何変分についての量である。
